using Oscar;

include("imports.jl")

function greedy_vertex_alignment(system::Vector,scoring_function,polRing=polRing,grad=true,ord=false,shuff=false)
    # Heuristic sorting by most to least vertices
    vertices=[collect(Nemo.exponent_vectors(f)) for f in system]
    perm=sortperm(vertices,rev=ord)
    if !shuff
        # add lex tiebreaker
        vertices=vertices[perm]
        system=system[perm]
    end
    for i in 1:(length(system)-1)
        # Here we deal *only* with the exponent vectors
        max_monomial=[0 for i in 1:length(system)]
        max_matrix=matrix_from_system(system)
        max_score=scoring_function(max_matrix)
        for v1 in union(vertices[1:i]...)
            for v2 in vertices[i+1]
                lowest_term=abs(sort(v2-v1)[1])
                monomial=prod(gens(polRing).^(v2-v1+[lowest_term for j in v2]))
                mod_system=system[:]
                deleteat!(mod_system,1:i)
                global_translation=prod(gens(polRing).^[lowest_term for j in v2])
                mod_system=[f*global_translation for f in mod_system]
                push!(mod_system,([monomial for j in 1:i].*system[1:i])...)
                mat=matrix_from_system(mod_system)
                score=scoring_function(mat)
                if score>max_score
                    max_score=score
                    max_monomial=[[lowest_term for j in v2]+Int(i>=j)*(v2-v1) for j in 1:length(system)]
                    max_matrix=matrix_from_system(system.*[prod(gens(polRing).^m) for m in max_monomial])
                elseif score==max_score
                    bound=100
                    if random_minors(mat,bound)>random_minors(max_matrix,bound)
                        max_matrix=mat
                        max_score=score
                        max_monomial=[[lowest_term for j in v2]+Int(i>=j)*(v2-v1) for j in 1:length(system)]
                    end
                end
            end
        end
        system=system.*[prod(gens(polRing).^m) for m in max_monomial]
        vertices=[collect(Nemo.exponent_vectors(f)) for f in system]
        perm=sortperm(vertices,rev=ord)
        vertices=vertices[perm]
    end
    # do a local maxima check TEST
    if grad
        perturbs=[j[1] for j in perturbSystem(system,polRing)]
        sort!(perturbs,by = x->scoring_function(matrix_from_system(x)))
        max_score=scoring_function(matrix_from_system(system))
        while scoring_function(matrix_from_system(perturbs[end]))>max_score
            system=perturbs[end]
            unsorted_perturbs=perturbs[j[1] for j in perturbSystem(system,polRing)]
            perturbs=sort(unsorted_perturbs,by = x->scoring_function(matrix_from_system(x)))
            max_score=scoring_function(matrix_from_system(system))
        end
    end
    return system
end

function testing(syss::Vector,lazy=true)
    fails=[]
    for sys in syss
        println(sys.ID)
        lowerbound=score8(matrix_from_system(sys.generic_polynomial_system))
        for x in 1:10
            println(x)
            per=bigperturb(sys.generic_polynomial_system,sys.polRing)
            worstbound=score8(matrix_from_system(per))
            modsys=greedy_vertex_alignment(per,score8,sys.polRing,lazy)
            midbound=score8(matrix_from_system(modsys))
            mid1=copy(midbound)
            println([-lowerbound,-midbound,-worstbound])
            modsys=greedy_vertex_alignment(per,score8,sys.polRing,lazy,true)
            midbound=score8(matrix_from_system(modsys))
            mid2=copy(midbound)
            println([-lowerbound,-midbound,-worstbound])
            modsys=greedy_vertex_alignment(per,score8,sys.polRing,lazy,true,true)
            println([-lowerbound,-midbound,-worstbound])
            push!(fails,mid2>mid1)
            push!(fails,-Int(mid2<mid1))
        end
    end
    println([sum(fails),10*length(syss)])
end

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

function perturbSystem(system::Vector,polring)
    trans=[]
    for m in gens(polring)
        for k in 1:length(system)
            perturb=[m^Int(j==k) for j in 1:length(system)]
            minperturb=[m^Int(j!=k) for j in 1:length(system)]
            push!(trans,perturb,minperturb)
        end
    end 
    push!(trans,[1 for j in system])
    explodedSystems=[[i.*system,i] for i in trans]
    return explodedSystems
end

function bigperturb(system::Vector,polring)
    trans=[prod(gens(polring).^rand(UInt8,length(gens(polring)))) for i in system]
    explodedSystems=[trans[i].*system[i] for i in 1:length(system)]
    return explodedSystems
end

# Could also work over finite field to make det computations faster (as all we care about is whether a determinant is zero) TODO
function randbinom(arr::Vector,num::Int)
    selection=Set()
    while length(selection)<num
        x=rand(arr)
        if !(x in selection)
            push!(selection,x)
        end
    end
    return selection
end

function random_minors(mat::QQMatrix,num::Int)
    num=min(num,binomial(number_of_columns(mat),number_of_rows(mat))//100)
    submats=Set()
    while length(submats) < num
        push!(submats, randbinom([mat[:,i] for i in 1:number_of_columns(mat)],number_of_rows(mat)))
    end
    minorszero=sum(is_det_zero.([matrix(QQ,hcat(m...)) for m in submats]))
end

# We work globally with QQMatrix since our generic_polynomial_system has rational coefficients
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

function is_det_zero(mat::QQMatrix)
    # returns "no method matching AbstractFloat"??
    if number_of_columns(mat)==number_of_rows(mat)
        return !(rank(mat)==number_of_rows(mat))
    end
    error("not square")
end
