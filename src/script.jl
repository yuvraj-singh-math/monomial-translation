using Oscar;
using OscarODEbase;

function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[QQ(coeff(f,m)) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    rk,M=rref(M)
    M=matrix(QQ,[M[i,:] for i in 1:rk])
    return M
end

# We return the 2n perturbations of degree 1 as well as the system itself
function perturbSystem(system::ODEbaseModel)
    trans=[]
    sys,specializationHomomorphism=get_polynomials_random_specialization(system,reduce=true)
    # remove 0 polynomials
    filter!(x->!iszero(x),sys)
    for m in gens(codomain(specializationHomomorphism))
        for k in 1:length(sys)
            perturb=[m^Int(j==k) for j in 1:length(sys)]
            minperturb=[m^Int(j!=k) for j in 1:length(sys)]
            push!(trans,perturb,minperturb)
        end
    end 
    push!(trans,[1 for j in sys])
    explodedSystems=[[i.*sys,i] for i in trans]
    return explodedSystems
end

function IDToODE(ID)
    return filter(m->m.ID==ID,systems)[1]
end

# We work globally with QQMatrix since our get_polynomials_random_specialization has rational coefficients
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
    numallminors=length(collect(allminors))
    allminors=nothing
    numzerominors=length(collect(zerominors))
    zerominors=nothing
    return [numallminors,numzerominors]
end


function data_dump(sys::ODEbaseModel;beginning=1,slow=false)
    matrix=[]
    perturbations=perturbSystem(sys)
    name=sys.ID
    len=length(perturbations)
    #if @isdefined slow
        if slow
            fulsup=fully_supported_minors_nonsparse_collection
        else
            fulsup=fully_supported_minors
        end
    #else
    fulsup=fully_supported_minors
    #end
    i=beginning
    for per in perturbations[beginning:end]
        #println("Perturbation $i/$len")
        mat=matrix_from_system(per[1])
        #println("Minors computed. Now filtering for zero determinants")
        working=true
        try
            global time=@elapsed begin
                global numRelevantMinors,numZeroRelevantMinors=fulsup(mat)
            end
        catch e
            println("$name takes too much memory.")
            working=false
        end
        if working
            numColumns=number_of_columns(mat)
            numMinors=binomial(numColumns,number_of_rows(mat))
            numZeroMinors=numMinors-numRelevantMinors+numZeroRelevantMinors
            row=["["*join(string.(per[2]), " ")*"]",numRelevantMinors,numZeroRelevantMinors,numZeroMinors,numMinors,numColumns]
            push!(matrix,row)
            log=open("out/output.log","a")
            println(log,"[$name:$i/$len]: $time")
            close(log)
            println("[$name:$i/$len]: $time")
            i=i+1
        end
    end
    return matrix
end
