module MonomialTranslations
export greedy_vertex_alignment
export produce_data
using Oscar;
using OscarODEbase;
const dir = Base.pkgdir(MonomialTranslations)

score(mat)=-number_of_columns(mat)

unfiltered_systems=get_odebase_model.(ODEbaseModels)

unfiltered_systems=sort(unfiltered_systems,by= x->x.numSpecies);
#unfiltered_systems=filter(s->s.numSpecies<=bound,unfiltered_systems);

#unfiltered_systems=filter(s->s.numSpecies<=16,unfiltered_systems);
unfiltered_systems=filter(s->s.massAction,unfiltered_systems);
rejects=Dict()
for sys in unfiltered_systems
    # Note that for f=2*x1, even though this is monomial, and has no toric solutions, is_monomial returns false
    # so we look at the length of the list of monomials, check if its 1
    if sum([length(collect(monomials(f)))==1 for f in generic_polynomial_system(sys)[1]])>0
        rejects[sys.ID]="Contains monomial equation";
        filter!(s->s.ID!=sys.ID,unfiltered_systems);
    end
end

systems=copy(unfiltered_systems);

include("script.jl")

#all_models=get_odebase_model.(ODEbaseModels)
#systems=generic_polynomial_system.(filter(x->x.massAction,get_odebase_model.(all_models)))

# Could also work over finite field to make det computations faster (as all we care about is whether a determinant is zero) TODO
function randbinom(arr::Vector,num::Int)
    selection=Set()
    while length(selection)<num
        x=rand(arr)
        push!(selection,x)
    end
    return selection
end

function random_minors(mat::QQMatrix,num::Int)
    num=min(num,number_of_rows(mat)//10)
    submats=Set()
    while length(submats) < num
        push!(submats, randbinom([mat[:,i] for i in 1:number_of_columns(mat)],number_of_rows(mat)))
    end
    minorszero=sum(is_det_zero.([matrix(QQ,hcat(m...)) for m in submats]))
end
function is_det_zero(mat::QQMatrix)
    if number_of_columns(mat)==number_of_rows(mat)
        return !(rank(mat)==number_of_rows(mat))
    end
    error("not square")
end

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

function vertices_of_function(f)
    points=collect(Oscar.vertices(newton_polytope(f)))
    points=[Int.(vertex) for vertex in points]
    return points
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

function greedy_vertex_alignment(system::Vector,scoring_function,grad=false,ord=false,shuff=false,tiebreaker=true)
    polRing=parent(system[1])
    points=[ vertices_of_function(f) for f in system]
    perm=sortperm(points,rev=ord)
    if !shuff
        # add lex tiebreaker
        points=points[perm]
        system=system[perm]
    end
    for i in 1:(length(system)-1)
        # Here we deal *only* with the exponent vectors
        max_monomial=[0 for i in 1:length(system)]
        max_matrix=matrix_from_system(system)
        max_score=scoring_function(max_matrix)
        for v1 in union(points[1:i]...)
            for v2 in points[i+1]
                # we adjust the difference to ensure it is positive
                # and translate every other polynomial accordingly
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
                elseif score==max_score && tiebreaker
                    # computing all the minors is too intensive so sample uniformly
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
        points=[ vertices_of_function(f) for f in system]
        perm=sortperm(points,rev=ord)
        points=points[perm]
    end
    # do a local maxima check TEST
    if grad
        perturbs=[j[1] for j in perturbSystem(system,polRing)]
        sort!(perturbs,by = x->scoring_function(matrix_from_system(x)))
        max_score=scoring_function(matrix_from_system(system))
        while scoring_function(matrix_from_system(perturbs[end]))>max_score
            system=perturbs[end]
            max_score=scoring_function(matrix_from_system(system))
            unsorted_perturbs=[j[1] for j in perturbSystem(system,polRing)]
            perturbs=sort(unsorted_perturbs,by = x->scoring_function(matrix_from_system(x)))
        end
    end
    return system
end


function produce_data(bound=16,restrict=false)
    if restrict
        filter!(x->x.numSpecies==bound,systems)
    else
        filter!(x->x.numSpecies<=bound,systems)
    end
    mkpath("out")
    mkpath("out/perturb_info")
    open("out/rejects.jl","w") do io
        println(io,rejects)
    end
    count=1
    total=length(systems)
    for sys in systems
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
        open("out/perturb_info/$name-matrix.csv", "w") do io
            write(io, file)
        end
        log=open("out/output.log","a")
        println("[$name:TOTAL]: $time")
        println(log,"[$name:TOTAL]: $time")
        close(log)
        count=count+1
    end
end
produce_data()
end
