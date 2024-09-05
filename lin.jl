using Oscar;
# Making it in a separate module causes a segfault. why??? TODO
# using odebase;

# map systems to their newton polytope
# identify common edges
# greed: prefer edges that appear in more vs appear in less
# identify translations (iso to monomials) required to attain these
# greed, prefer monomials which appear more vs less
# identify circuits of associated matrix by columns, try to maximise


# form matrix out of systems
# form n-minors out of this matrix

# 092.jl is an example of a system where the number of rows exceed the number of columns 
# 40.jl has 0 fully supported zero minors, as with the degree 1 pertubations
# 233 has 0
# 267, rows>columns
# 283, rows>columns
# 363, rows>columns
# 854, wins on all counts. every fully supported minor is of det zero (it is square)
# 609 has 0
# 629, square with det 0, so it wins score4, ties score2, and loses the rest
# 868, wins score4

function score1(mat::QQMatrix)
    return number_of_zero_minors(mat)
end

function score2(mat::QQMatrix)
    return number_of_zero_minors(mat)/number_of_fully_supported_minors(mat)
end

function score3(mat::QQMatrix)
    return number_of_zero_minors(mat)/number_of_columns(mat)
end

# Ratio of number of fully supported zero minors to all minors is the winner for 413.jl
function score4(mat::QQMatrix)
    # TODO could get imprecise as systems get larger, fix
    return number_of_zero_minors(mat)/binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

function score5(mat::QQMatrix)
    return 1/binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

# prioritises the one with the least number of minors which are not fully supported zero minors
function score6(mat::QQMatrix)
    return number_of_zero_minors(mat)-binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

function score7(mat::QQMatrix)
    return number_of_zero_minors(mat)-number_of_fully_supported_minors(mat)
end


scoring_functions=[score1,score2,score3,score4,score5,score6,score7]

#system=[x^2+y^2+2x+3y+7x^2*y^3,x^2+y^7+2x^2+3y^3+7x+2y+x^7*y^7,x^2+y^2+1]
#system=[x^5*(1+x+y+x*y),x+y+x*y,x^3+y]
# system=[58*x^4*y^6 + 97*x^3*y^9, 12*x^8*y^9 + 92*x^5*y^9 + 29*x^5*y^2 + 99*x^4*y^5 + 92*x^4*y + 76*x^2*y^3 + 44*x*y^2, 27*x^10*y^5 + 87*x^9*y^8 + 70*x^8*y^5 + 100*x^3*y^8 + 36*x^2*y^8, 74*x^9*y^9 + x^8*y^2 + 22*x^5*y^4 + 7*x^4*y^8 + 2*x^4*y^3 + 71*x^4*y^2, 13*x^10*y + 95*x^7*y^2 + 18*x^7*y + 35*x^6*y^9 + 3*x^2*y^7, 79*x^9*y^9 + 98*x^8*y + 87*x^6*y^8 + 15*x^5*y^5, 99*x^10*y^7 + 49*x^8*y^7 + 100*x^7*y^9 + 7*x^6*y^5 + 83*x^6*y + 15*x^5*y^5 + 61*x^4*y^9 + 64*x^4*y^3, 73*x^9*y^2 + 79*x^7*y^10 + 12*x^7*y^9 + 67*x^7*y^8 + 15*x^4*y^6, 58*x^5*y^6, 61*x^10*y^8 + 18*x^8*y^9 + 85*x^8*y^4 + 34*x^5*y^9 + 48*x^3*y^7 + 9*x^3*y^3 + 80*x*y^4]

function greedy_vertex_alignment(system::Vector,scoring_function)
    # sorting matters for greedy algo, consider the best one TODO
    vertices=[collect(Nemo.exponent_vectors(f)) for f in system]
    # sort vertices and then the system by ascending number of monomials
    perm=sortperm(vertices)
    vertices=vertices[perm]
    system=system[perm]
    # iterate through the monomials of system i compared to i+1
    # find optimal path between these, greedy choice is to combine these all
    # if m is max monomial length, and n the number of polynomials, max paths computed is (n-1)*m^2
    # going through all of them would yield m^n paths
    # leave last fixed, we transform the rest
    # our score is the norm of a vector of number of relevant zero minors
    # and the reciprocal of the number of columns
    base_matrix=matrix_from_system(system)
    base_score=scoring_function(base_matrix)
    for i in 1:(length(system)-1)
        max_score=base_score
        # Here we deal *only* with the exponent vectors
        max_monomial=[0 for i in 1:length(system)]
        for v1 in vertices[i]
            for v2 in vertices[i+1]
                lowest_term=sort(v2-v1)[1]
                monomial=prod(gens(polRing).^(v2-v1+[abs(lowest_term) for j in v2]))
                mod_system=system[:]
                deleteat!(mod_system,i)
                global_translation=prod(gens(polRing).^[abs(lowest_term) for j in v2])
                mod_system=[f*global_translation for f in mod_system]
                mod_system=push!(mod_system,monomial*system[i])
                mat=matrix_from_system(mod_system)
                score=scoring_function(mat)
                if score>max_score
                    max_score=score
                    max_monomial=[[abs(lowest_term) for j in v2]+Int(i==j)*(v2-v1) for j in 1:length(system)]
                end
            end
        end
        system=system.*[prod(gens(polRing).^m) for m in max_monomial]
        # do a local maxima check TODO
        perturbs=[j[1] for j in perturbSystem(system,polRing)]
        sort!(perturbs,by = x->scoring_function(matrix_from_system(x)))
        while scoring_function(matrix_from_system(perturbs[end]))>max_score
            system=perturbs[end]
            unsorted_perturbs=perturbs=[j[1] for j in perturbSystem(system,polRing)]
            perturbs=sort(unsorted_perturbs,by = x->scoring_function(matrix_from_system(x)))
            max_score=scoring_function(matrix_from_system(system))
        end
    end
    return system
end

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

# we do not want to remove any coefficients entirely
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

