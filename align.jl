using Oscar;
using OscarODEbase;
using Random;
using MonomialTranslation;
Random.seed!(13371337)
Oscar.set_seed!(13371337)

score(mat)=-number_of_columns(mat);

systems=get_odebase_model.(ODEbaseModels);
systems=filter(x->x.massAction,systems);
sort!(systems,by=x->x.numSpecies);
systems=filter(x->x.numSpecies<=16,systems);


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

function bigperturb(system::Vector,polring)
    trans=[prod(gens(polring).^rand(UInt8,length(gens(polring)))) for i in system]
    explodedSystems=[trans[i].*system[i] for i in 1:length(system)]
    return explodedSystems
end

function testing(syss::Vector,lazy=true)
    results=[]
    for sys in syss
        println(sys.ID)
        gen_sys,newring=generic_polynomial_system(sys)
        filter!(x->!iszero(x),gen_sys)
        gen_sys=unique(gen_sys)
        lowerbound=score(matrix_from_system(gen_sys))
        for x in 1:10
            per=bigperturb(gen_sys,newring)
            worstbound=score(matrix_from_system(per))
            modsys=greedy_vertex_alignment(per,score,lazy,false,false,false)
            midbound=score(matrix_from_system(modsys))
            println([-lowerbound,-midbound,-worstbound])
            push!(results,[sys.ID,-lowerbound,-midbound,-worstbound])
        end
    end
    return results
end

outputdata=testing(systems);

global csv=[]
for result in outputdata
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
