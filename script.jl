using Oscar;
using Combinatorics;
# try to find better implementation
using BipartiteMatching;

if length(ARGS)==2
    global bound::Int=parse(Int,ARGS[1])
    global slow::Bool=parse(Bool,ARGS[2])
end
include("imports.jl")

# filter out systems that we have already computed
unfiltered_systems=filter(sys->!isfile("odebase/out/$(sys.ID)-matrix.csv"),unfiltered_systems)

# quick approximation for complexity. better would be to compute the upper bound on no. of fully supported minors
## uncomment the two lines below to filter out any systems with number of species>=upperBound
if @isdefined bound
unfiltered_systems=filter(s->s.numSpecies==bound,unfiltered_systems);
else
unfiltered_systems=sort(unfiltered_systems,by= x->x.numSpecies);
end

unfiltered_systems=filter(s->s.numSpecies<=16,unfiltered_systems);
function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[QQ(coeff(f,m)) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    if number_of_columns(M)<number_of_rows(M)
        rk,M=rref(M)
        M=matrix(QQ,[M[i,:] for i in 1:rk])
    end
    return M
end

unfiltered_systems=sort(unfiltered_systems,by= x->number_of_columns(matrix_from_system(x.generic_polynomial_system)));

systems=copy(unfiltered_systems);

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
function is_det_zero(mat::QQMatrix)
    # returns "no method matching AbstractFloat"??
    if number_of_columns(mat)==number_of_rows(mat)
        return !(rank(mat)==number_of_rows(mat))
    end
    error("not square")
end

# uses too much memory
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
    niceminors=filter(m->is_det_zero(matrix(QQ,hcat(m...))),minors)
    return [length(minors),length(niceminors)]
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

function is_minor_fully_supported(submat::Vector{Vector{QQFieldElem}})
    dimen=length(submat)
    adj_matrix=BitArray{2}(hcat([[!is_zero(submat[i][j]) for i in 1:dimen] for j in 1:dimen]...))
    a,b=findmaxcardinalitybipartitematching(adj_matrix)
    # b is a bitvector where the i-th component is true (resp. false) if the i-th row is (resp. is not) matched
    return prod(b)
end

function fully_supported_minors_nonsparse(mat::QQMatrix)
    cols=[mat[:,i] for i in 1:number_of_columns(mat)];
    allminors=Combinatorics.combinations(cols,number_of_rows(mat));
    allminors=Iterators.filter(x->is_minor_fully_supported(x),allminors);
    zerominors=Iterators.filter(m->is_det_zero(matrix(QQ,hcat(m...))),allminors);
    # try to get around having to collect
    numallminors=0
    numzerominors=0
    for x in allminors
        numallminors=numallminors+1
    end
    for x in zerominors
        numzerominors=numzerominors+1
    end
    return [numallminors,numzerominors]
end

function fully_supported_minors_nonsparse_collection(mat::QQMatrix)
    cols=[mat[:,i] for i in 1:number_of_columns(mat)];
    allminors=Combinatorics.combinations(cols,number_of_rows(mat));
    allminors=Iterators.filter(x->is_minor_fully_supported(x),allminors);
    zerominors=Iterators.filter(m->is_det_zero(matrix(QQ,hcat(m...))),allminors);
    # try to get around having to collect
    numallminors=length(collect(allminors))
    allminors=nothing
    numzerominors=length(collect(zerominors))
    zerominors=nothing
    return [numallminors,numzerominors]
end


function data_dump(sys::OdebaseNode)
    matrix=[]
    perturbations=perturbSystem(sys)
    name=sys.ID
    i=1
    len=length(perturbations)
    if @isdefined slow
        if slow
            fulsup=fully_supported_minors_nonsparse
        else
            fulsup=fully_supported_minors
        end
    else
        fulsup=fully_supported_minors
    end
    
    for per in perturbations
        #println("Perturbation $i/$len")
        mat=matrix_from_system(per[1])
        #println("Minors computed. Now filtering for zero determinants")
        time=@elapsed begin
        numRelevantMinors,numZeroRelevantMinors=fulsup(mat)
        end
        numColumns=number_of_columns(mat)
        numMinors=binomial(numColumns,number_of_rows(mat))
        numZeroMinors=numMinors-numRelevantMinors+numZeroRelevantMinors
        row=["["*join(string.(per[2]), " ")*"]",numRelevantMinors,numZeroRelevantMinors,numZeroMinors,numMinors,numColumns]
        push!(matrix,row)
        log=open("output.log","a")
        println(log,"[$name:$i/$len]: $time")
        close(log)
        i=i+1
    end
    #matrix=Matrix(transpose(hcat(matrix...)))
    return matrix
end

open("odebase/src/rejects.jl","w") do io
    println(io,rejects)
end

count=1
total=length(systems)
for sys in systems
    global count
    local name=sys.ID
    print(name)
    println(", system $count/$total:")
    time=@elapsed begin
    data=data_dump(sys)
    end
    num=length(data[1])
    i=1
    csv=[]
    for result in data
        for j in 1:num
            push!(csv,"$(string(result[j]))")
            if j+1<=num
                push!(csv,",")
            end
        end
        push!(csv,"
")
    end
    file=string(csv...)
    open("odebase/out/$name-matrix.csv", "w") do io
        write(io, file)
    end
    log=open("output.log","a")
    println(log,"[$name:TOTAL]: $time")
    close(log)
    count=count+1
end
