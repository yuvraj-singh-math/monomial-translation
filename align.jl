using Oscar;
using OscarODEbase;
using Random;
using MonomialTranslations;
Random.seed!(13371337)
Oscar.set_seed!(13371337)

score(mat)=-number_of_columns(mat);

systems=get_odebase_model.(ODEbaseModels);

systems=filter(x->x.massAction,systems);
sort!(systems,by=x->x.numSpecies);
systems=systems[1:end-1]

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

function bigperturb(system::Vector,specializationHom)
    trans=[prod(gens(codomain(specializationHom)).^rand(UInt8,length(gens(codomain(specializationHom))))) for i in system]
    #println(first(trans))
    #println(first(system))
    explodedSystems=[trans[i].*system[i] for i in 1:length(system)]
    return explodedSystems
end

function testing(syss::Vector,lazy=true)
    stats=[]
    results=[]
    for sys in syss
        println(sys.ID)
        gen_sys,specializationHom=get_polynomials_random_specialization(sys,reduce=true);
        filter!(x->!iszero(x),gen_sys)
        lowerbound=score(matrix_from_system(gen_sys))
	results_sys=[]
        for x in 1:10
            per=bigperturb(gen_sys,specializationHom)
            worstbound=score(matrix_from_system(per))
            time=@elapsed begin
            modsys=greedy_vertex_alignment(per);
            end
            println(time)
            midbound=score(matrix_from_system(modsys))
            println([-lowerbound,-midbound,-worstbound])
            push!(results,[sys.ID,-lowerbound,-midbound,-worstbound])
            push!(results_sys,[sys.ID,-lowerbound,-midbound,-worstbound])
        end
	average=sum([result[3] for result in results_sys])/10
	highest=max([x[3] for x in results_sys]...)
	lowest=min([x[3] for x in results_sys]...)
	push!(stats,[sys.ID,-lowerbound,lowest,average,highest,results_sys[1][4]])
    end
    return results,stats
end

outputdata,stats=testing(systems);

global csv_raw=[]
global csv=[]
for result in outputdata
    for j in 1:length(result)
        push!(csv_raw,"$(string(result[j]))")
        if j+1<=length(result)
            push!(csv_raw,",")
        end
    end
    push!(csv_raw,"
")
end
file_raw=string(csv_raw...);
open("alignment-results_raw.csv", "w") do io
    write(io, file_raw)
end

for result in stats
    for j in 1:length(result)
        push!(csv,"$(string(result[j]))")
        if j+1<=length(result)
            push!(csv,",")
        end
    end
    push!(csv,"
")
end
file=string(csv...);
open("alignment-results.csv", "w") do io
    write(io, file)
end
