module MonomialTranslations
export greedy_vertex_alignment
export produce_data
using Oscar;
using OscarODEbase;
const dir = Base.pkgdir(MonomialTranslations)


function __init__()
	unfiltered_systems=get_odebase_model.(ODEbaseModels)
	unfiltered_systems=sort(unfiltered_systems,by= x->x.numSpecies);
	unfiltered_systems=filter(s->s.massAction,unfiltered_systems);
	global rejects=Dict()
	for sys in unfiltered_systems
    # Note that for f=2*x1, even though this is monomial, and has no toric solutions, is_monomial returns false
    # so we look at the length of the list of monomials, check if its 1
	    if sum([length(collect(monomials(f)))==1 for f in generic_polynomial_system(sys)[1]])>0
		rejects[sys.ID]="Contains monomial equation";
		filter!(s->s.ID!=sys.ID,unfiltered_systems);
	    end
	end
	global systems=copy(unfiltered_systems);
end


#unfiltered_systems=filter(s->s.numSpecies<=bound,unfiltered_systems);

#unfiltered_systems=filter(s->s.numSpecies<=16,unfiltered_systems);
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

function vertices_of_function(f)
    points=collect(Oscar.vertices(newton_polytope(f)))
    points=[Int.(vertex) for vertex in points]
    return points
end

# To deal with generic polynomial systems
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

function integer_stuff(points::Vector{Vector{Vector{Int}}},lensys::Int,ind::Int)
    # Here we deal *only* with the exponent vectors
    max_monomial=[0 for j in 1:lensys]
    max_score=-length(union(points...))
    for v1 in union(points[1:ind]...)
        for v2 in points[ind+1]
            # we adjust the difference to ensure it is positive
            # and translate every other polynomial accordingly
            lowest_term=abs(sort(v2-v1)[1])
            monomial=[[lowest_term for j in v2]+Int(ind>=j)*(v2-v1) for j in 1:lensys]
            mod_points=[[mon+monomial[j] for mon in points[j]] for j in 1:length(points)]
            # union of the points is the same as number of monomials, which is the same as number of columns of the Macaulay matrix
            score=-length(union(mod_points...))
            if score>max_score
                max_score=score
                max_monomial=[[lowest_term for j in v2]+Int(ind>=j)*(v2-v1) for j in 1:lensys]
            end
        end
    end
    return max_monomial
end

# ,grad::Bool=false,ord::Bool=false,shuff::Bool=false,tiebreaker::Bool=true
function greedy_vertex_alignment(system::Vector)
    polRing=parent(system[1])
    points=[ vertices_of_function(f) for f in system]
    perm=sortperm(points)
    points=points[perm]
    system=system[perm]
    for i in 1:(length(system)-1)
        max_monomial=integer_stuff(points,length(system),i)
        system=system.*[prod(gens(polRing).^m) for m in max_monomial]
        points=[ vertices_of_function(f) for f in system]
    end
    return system
end


function produce_data(bound=16,restrict=false;start=1)
    if restrict
        filter!(x->x.numSpecies==bound,systems)
    else
        filter!(x->x.numSpecies<=bound,systems)
    end
    systems=systems[start:end]
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
end
