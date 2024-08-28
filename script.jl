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
    constraints
    paramsRing
    polRing
end

# The initial values for rejects are defined by those systems that have a parameter to the power of another parameter (eg k1^k2)
# We do not save these as .jl file to begin with as of right now
rejects=Dict(vcat([id=>"Fails to load in Julia" for id in ["BIOMD0000000060","BIOMD0000000637"]],["BIOMD0000000205"=>"Contains monomial equation"])...)
odebaseSystems=[]

for file in chems
    include(file)
    randCoeff=rand(Int8,length(gens(paramsRing)))
    phi=hom(polRing,polRing,c->evaluate(c,randCoeff),gens(polRing))
    push!(odebaseSystems, OdebaseNode(name,true,true,1,1,length(gens(polRing)),chemSystem,[phi(x) for x in chemSystem],[],paramsRing,polRing))
end

unfiltered_systems=[OdebaseNode(sys.ID,true,true,1,1,sys.numSpecies,sys.param_polynomial_system,filter(l->!iszero(l),unique(sys.generic_polynomial_system)),[],sys.paramsRing,sys.polRing) for sys in odebaseSystems]

for sys in odebaseSystems
    # Note that for f=2*x1, even though this is monomial, and has no toric solutions, is_monomial returns false
 #    so we look at the length of the list of monomials, check if its 1
    if sum([length(collect(monomials(f)))==1 for f in sys.generic_polynomial_system])>0
        rejects[sys.ID]="Contains monomial equation"
        filter!(s->s!=sys,unfiltered_systems)
    end

    if length(unique(collect(Iterators.flatten([collect(monomials(f)) for f in sys.generic_polynomial_system]))))<length(sys.generic_polynomial_system)
        rejects[sys.ID]="Macaulay matrix has more rows than columns"
        filter!(s->s!=sys,unfiltered_systems)
    end
end

# unfiltered_systems=filter(s->s.numSpecies<6,unfiltered_systems)
const systems=copy(unfiltered_systems)

# We return the 2n pertubations of degree 1 as well as the system itself
function perturbSystem(system)
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


function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[coeff(f,m) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    return M
end


function is_det_zero(mat)
    # returns "no method matching AbstractFloat"??
    if number_of_columns(mat)==number_of_rows(mat)
        return !(rank(mat)==number_of_rows(mat))
    end
    error("not square")
end

function niceprod(arrcols,illegal)
    paths=[]
    if length(arrcols)>1
        for element in arrcols[1]
           push!(paths,vcat.(Ref(element),niceprod(arrcols[2:end],vcat(element,illegal)))...)
        end
        return paths
    else
        paths=filter(m->!(m in illegal),arrcols[1])
    end
    return paths
end

function fully_supported_minors(mat)
    # We tag each column to deal with two that may have identical entries
     cols=[[i,mat[:,i]] for i in 1:number_of_columns(mat)]
    cols_per_sys=[filter(m->!iszero(m[2][i]),cols) for i in 1:number_of_rows(mat)];
    # OLD: very slow, but doesn't crash julia for larger systems
    #degenminors=Iterators.product(cols_per_sys...)
    ## delete those with repeat columns since they will always be det 0
    ## used to collect first and then filter, this gives OOM
    #dupeminors=Iterators.filter(m->length(m)==length(unique(m)),degenminors)
    #minors=collect(dupeminors)
    # NEW: faster (?) but does crash julia on large systems. make nonrecursive TODO
    dupeminors=niceprod(cols_per_sys,[])

    # select minors unique up to pertubation of columns
    minors=unique(Set,dupeminors)
    minors=[[c[2] for c in m] for m in minors]
    return minors
end

# returns number of zero minors that we consider
function number_of_zero_minors(mat)
    minors=fully_supported_minors(mat)
    niceminors=filter(m->is_det_zero(matrix(QQ,hcat(m...))),minors)
    return length(niceminors)
end

function number_of_fully_supported_minors(mat)
    return length(fully_supported_minors(mat))
end

function data_dump_matrix(sys)
    matrix=[]
    for per in perturbSystem(sys)
        mat=matrix_from_system(per[1])
        numRelevantMinors=number_of_fully_supported_minors(mat)
        numZeroMinors=number_of_zero_minors(mat)
        numColumns=number_of_columns(mat)
        numMinors=binomial(max(numColumns,number_of_rows(mat)),min(numColumns,number_of_rows(mat)))
        row=[per[2],numRelevantMinors,numZeroMinors,numMinors,numColumns]
        push!(matrix,row)
    end
    matrix=Matrix(transpose(hcat(matrix...)))
    return matrix
end

for sys in systems
    name=sys.ID
    print(name)
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
end

open("odebase/src/rejects.jl","w") do io
    println(io,rejects)
end
