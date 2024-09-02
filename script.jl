using Oscar;

chems=readdir("odebase/src",join=true)
# these filters are mostly hacky workarounds
chems=filter(filename->occursin(".jl",filename)&&(!occursin("odebase.jl",filename))&&(!occursin("rejects.jl",filename))&&(!occursin("matrix",filename))&&!occursin("#",filename),chems)

struct OdebaseNode
    ID::String
    rational::Bool
    massAction::Bool
    #? redundant
    species::Int
    deficit::Int
    numSpecies::Int
    # TODO types for the following
    # \dot{x}_i is set to 0
    param_polynomial_system::Vector
    # should replace this with a function soon
    generic_polynomial_system::Vector
    constraints::Vector
    paramsRing
    # for now, rational (remove when not in script)
    polRing::QQMPolyRing
end

# The initial values for rejects are defined by those systems that have a parameter to the power of another parameter (eg k1^k2)
# We do not save these as .jl file to begin with as of right now
rejects=Dict(vcat([id=>"Fails to load in Julia" for id in ["BIOMD0000000060","BIOMD0000000637"]],["BIOMD0000000205"=>"Contains monomial equation"])...)
odebaseSystems=OdebaseNode[]

function rand_nonzero(len::Int)
    ints=Int[]
    for x in 1:len
        num=rand(Int8)
        while num==0
            num=rand(Int8)
        end
        push!(ints,num)
    end
    return ints
end

for file in chems
    include(file);
    randCoeff=rand_nonzero(length(gens(paramsRing)));
    QQpolRing,tup=polynomial_ring(QQ,["$x" for x in gens(polRing)]);
    phi=hom(polRing,QQpolRing,c->evaluate(c,randCoeff),gens(QQpolRing));
    # we redefine polRing to be of rational type after the map
    push!(odebaseSystems, OdebaseNode(name,true,true,1,1,length(gens(polRing)),chemSystem,[phi(x) for x in chemSystem],[],paramsRing,QQpolRing));
end

unfiltered_systems=[OdebaseNode(sys.ID,true,true,1,1,sys.numSpecies,sys.param_polynomial_system,filter(l->!iszero(l),unique(sys.generic_polynomial_system)),[],sys.paramsRing,sys.polRing) for sys in odebaseSystems];

for sys in odebaseSystems
    # Note that for f=2*x1, even though this is monomial, and has no toric solutions, is_monomial returns false
 #    so we look at the length of the list of monomials, check if its 1
    if sum([length(collect(monomials(f)))==1 for f in sys.generic_polynomial_system])>0
        rejects[sys.ID]="Contains monomial equation";
        filter!(s->s!=sys,unfiltered_systems);
    end

    if length(unique(collect(Iterators.flatten([collect(monomials(f)) for f in sys.generic_polynomial_system]))))<length(sys.generic_polynomial_system)
        rejects[sys.ID]="Macaulay matrix has more rows than columns";
        filter!(s->s!=sys,unfiltered_systems);
    end
end

# quick approximation for complexity. better would be to compute the upper bound on no. of fully supported minors
## uncomment the two lines below to filter out any systems with number of species>=upperBound
#upperbound=89
#unfiltered_systems=filter(s->s.numSpecies<upperbound,unfiltered_systems);
unfiltered_systems=sort(unfiltered_systems,by= x->x.numSpecies);
const systems=copy(unfiltered_systems);

# We return the 2n perturbations of degree 1 as well as the system itself
function perturbSystem(system::OdebaseNode)
    trans=[]
    for m in gens(system.polRing)
        for k in 1:length(system.generic_polynomial_system)
            perturb=[m^Int(j==k) for j in 1:length(system.generic_polynomial_system)]
            minperturb=[m^Int(j!=k) for j in 1:length(system.generic_polynomial_system)]
            push!(trans,perturb,minperturb)
        end
    end 
    push!(trans,[1 for j in system.generic_polynomial_system])
    explodedSystems=[[i.*system.generic_polynomial_system,i] for i in trans]
    return explodedSystems
end

function IDToODE(ID)
    return filter(m->m.ID==ID,systems)[1]
end


# We work globally with QQMatrix since our generic_polynomial_system has rational coefficients
function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[QQ(coeff(f,m)) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    return M
end


function is_det_zero(mat::QQMatrix)
    # returns "no method matching AbstractFloat"??
    if number_of_columns(mat)==number_of_rows(mat)
        return !(rank(mat)==number_of_rows(mat))
    end
    error("not square")
end

function niceprod_nonrecur(arrcols::Vector{Vector{Int}})
    # We start from the biggest set because this is most likely* to contain the most duplicates
    # * for ODEBASE systems with a lot of overlap
    # worst case is being of the same complexity as Iterators.product()
    arrcols=sort(arrcols,by= x->length(x), rev=true)
    paths=Set{Set{Int}}([Set{Int}([element]) for element in arrcols[1]])
    for set in arrcols[2:end]
        newpaths=Set{Set{Int}}()
        for path in paths
            for element in set
                if !(element in path)
                    push!(newpaths,union(path,element))
                end
            end
        end
        #println(n)
        paths=newpaths
        #println(length(paths))
    end
    return paths
end

# get strict w types
function fully_supported_minors(mat::QQMatrix)
    # We tag each column to deal with two that may have identical entries
    cols=[[i,mat[:,i]] for i in 1:number_of_columns(mat)];
    cols_per_sys=[filter(m->!iszero(m[2][i]),cols) for i in 1:number_of_rows(mat)];
    cols_per_sys_num=[[c[1] for c in sys] for sys in cols_per_sys]

    minors=niceprod_nonrecur(cols_per_sys_num)

    # Convert minors from lists of integers to lists of columns
    minors=Vector{Vector{QQFieldElem}}[[cols[i][2] for i in minor] for minor in minors]
    return minors
end

# returns number of zero minors that we consider
function number_of_zero_minors(mat::QQMatrix)
    minors=fully_supported_minors(mat)
    niceminors=filter(m->is_det_zero(matrix(QQ,hcat(m...))),minors)
    return length(niceminors)
end

# returns number of zero minors that we consider
function number_of_zero_minors_provided_minors(minors::Vector{Vector{Vector{QQFieldElem}}})
    niceminors=filter(m->is_det_zero(matrix(QQ,hcat(m...))),minors)
    return length(niceminors)
end

function number_of_fully_supported_minors(mat::QQMatrix)
    return length(fully_supported_minors(mat))
end

function data_dump_matrix(sys::OdebaseNode)
    matrix=[]
    perturbations=perturbSystem(sys)
    i=1
    len=length(perturbations)
    for per in perturbations
        println("Perturbation $i/$len")
        mat=matrix_from_system(per[1])
        fullySupportedMinors=fully_supported_minors(mat)
        println("Minors computed. Now filtering for zero determinants")
        numRelevantMinors=length(fullySupportedMinors)
        numZeroMinors=number_of_zero_minors_provided_minors(fullySupportedMinors)
        numColumns=number_of_columns(mat)
        numMinors=binomial(max(numColumns,number_of_rows(mat)),min(numColumns,number_of_rows(mat)))
        row=[per[2],numRelevantMinors,numZeroMinors,numMinors,numColumns]
        push!(matrix,row)
        i=i+1
    end
    matrix=Matrix(transpose(hcat(matrix...)))
    return matrix
end

open("odebase/src/rejects.jl","w") do io
    println(io,rejects)
end

count=1
total=length(systems)
for sys in systems
    global count
    name=sys.ID
    print(name)
    println(", system $count/$total:")
    matrix=data_dump_matrix(sys)
    num=number_of_columns(matrix)
    i=1
    csv=[]
    for i in 1:number_of_rows(matrix)
        for j in 1:num
            push!(csv,"$(matrix[i,j])")
            if j+1<=num
                push!(csv,",")
            end
        end
        push!(csv,"
")
    end
    file=string(csv...)
    open("odebase/src/$name-matrix.jl", "w") do io
        write(io, file)
    end
    count=count+1
end
