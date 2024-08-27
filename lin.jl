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

# TODO: better det verification
# TODO: filter out row>columns
# deg 1 pertubations

# paramsRing
# polRing

# map parameters to random values
# some systems have repeated equations, or zero equations, as they arise from ODEs

# right now structs are *mutable*. this is bad practice (?)
# TODO: come up with better implementation

chems=readdir("odebase/src",join=true)
chems=filter(filename->occursin(".jl",filename)&&(!occursin("odebase.jl",filename))&&!occursin("#",filename),chems)

# figure out way to select by no. of species, &c.
# compute degenericity of matrices

mutable struct OdebaseNode
    ID::String
    rational::Bool
    massAction::Bool
    species::Int
    deficit::Int
    numSpecies::Int
    # TODO types for the following
    # \dot{x}_i is set to 0
    param_polynomial_system::Vector
    generic_polynomial_system::Vector
    constraints
    paramsRing
    polRing
end

# fix 60, 637
odebaseSystems=[]
for file in chems
    include(file)
    randCoeff=rand(Int8,length(gens(paramsRing)))
    phi=hom(polRing,polRing,c->evaluate(c,randCoeff),gens(polRing))
    push!(odebaseSystems, OdebaseNode(name,true,true,1,1,length(gens(polRing)),chemSystem,[phi(x) for x in chemSystem],[],paramsRing,polRing))
end
randCoeff=nothing
phi=nothing
chems=nothing

systems=copy(odebaseSystems);
#systems=filter(m->m.numSpecies<10,systems)
# OOM errors when numParams is too large

# we want to filter out systems where rows>columns
systems=filter(m->length(unique(collect(Iterators.flatten([collect(monomials(f)) for f in m.generic_polynomial_system]))))>=length(m.generic_polynomial_system),systems);
for system in systems
    system.generic_polynomial_system=unique(system.generic_polynomial_system);
    system.generic_polynomial_system=filter(l->!iszero(l),system.generic_polynomial_system);
end

#expsystem=[prod(gens(polRing).^abs.(rand(Int8,length(gens(polRing)))))*f for f in system]

# we consider 2n transformations of the system, moving one component by gens(S)

function perturbSystem(system)
    trans=[]
    for m in gens(system.polRing)
        for k in 1:length(system.generic_polynomial_system)
            perturb=[m^Int(j==k) for j in 1:length(system.generic_polynomial_system)]
            minperturb=[m^Int(j!=k) for j in 1:length(system.generic_polynomial_system)]
            push!(trans,perturb,minperturb)
        end
    end 
    explodedSystems=[i.*system.generic_polynomial_system for i in trans]
    return explodedSystems
end

function IDToODE(ID)
    return filter(m->m.ID==ID,systems)[1]
end

# we change base ring for the code to work
#explodedSystems=[[change_base_ring(QQ,f) for f in sys] for sys in explodedSystems]
#system=[change_base_ring(QQ,f) for f in system]
#expsystem=[change_base_ring(QQ,f) for f in expsystem]


function score1(mat)
    return number_of_zero_minors(mat)
end

function score2(mat)
    return number_of_zero_minors(mat)/number_of_fully_supported_minors(mat)
end

function score3(mat)
    return number_of_zero_minors(mat)/number_of_columns(mat)
end

# Ratio of number of fully supported zero minors to all minors is the winner for 413.jl
function score4(mat)
    # TODO could get imprecise as systems get larger, fix
    return number_of_zero_minors(mat)/binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

function score5(mat)
    return 1/binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

# prioritises the one with the least number of minors which are not fully supported zero minors
function score6(mat)
    return number_of_zero_minors(mat)-binomial(max(number_of_columns(mat),number_of_rows(mat)),min(number_of_columns(mat),number_of_rows(mat)))
end

function score7(mat)
    return number_of_zero_minors(mat)-number_of_fully_supported_minors(mat)
end


scoring_functions=[score1,score2,score3,score4,score5,score6,score7]

function testing(systems,scores)
    total=length(systems)
    totalsuccess=[]
    for system in systems
        println(system.ID)
        original_scores=[s(matrix_from_system(system.generic_polynomial_system)) for s in scores]
        perturbations=perturbSystem(system)
        perturbed_scores=[[score(matrix_from_system(sys)) for score in scores] for sys in perturbations]
        # check whether the original system beats *every* pertubation
        # 1004 does not beat the pertubations, but is matched
        success=sum([Int.(original_scores.>m) for m in perturbed_scores])
        success=[Int(m==length(perturbations)) for m in success]
        push!(totalsuccess,success)
        println(success)
    end
    print("TOTAL SUCCESSES (max is ")
    print(total)
    println("):")
    totalsuccess=sum(totalsuccess)
    println(totalsuccess)
    return totalsuccess
end



#system=[x^2+y^2+2x+3y+7x^2*y^3,x^2+y^7+2x^2+3y^3+7x+2y+x^7*y^7,x^2+y^2+1]
#system=[x^5*(1+x+y+x*y),x+y+x*y,x^3+y]
# system=[58*x^4*y^6 + 97*x^3*y^9, 12*x^8*y^9 + 92*x^5*y^9 + 29*x^5*y^2 + 99*x^4*y^5 + 92*x^4*y + 76*x^2*y^3 + 44*x*y^2, 27*x^10*y^5 + 87*x^9*y^8 + 70*x^8*y^5 + 100*x^3*y^8 + 36*x^2*y^8, 74*x^9*y^9 + x^8*y^2 + 22*x^5*y^4 + 7*x^4*y^8 + 2*x^4*y^3 + 71*x^4*y^2, 13*x^10*y + 95*x^7*y^2 + 18*x^7*y + 35*x^6*y^9 + 3*x^2*y^7, 79*x^9*y^9 + 98*x^8*y + 87*x^6*y^8 + 15*x^5*y^5, 99*x^10*y^7 + 49*x^8*y^7 + 100*x^7*y^9 + 7*x^6*y^5 + 83*x^6*y + 15*x^5*y^5 + 61*x^4*y^9 + 64*x^4*y^3, 73*x^9*y^2 + 79*x^7*y^10 + 12*x^7*y^9 + 67*x^7*y^8 + 15*x^4*y^6, 58*x^5*y^6, 61*x^10*y^8 + 18*x^8*y^9 + 85*x^8*y^4 + 34*x^5*y^9 + 48*x^3*y^7 + 9*x^3*y^3 + 80*x*y^4]

function greedy_vertex_alignment(system,scoring_function)
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
    monomialTranslations=[]
    # our score is the norm of a vector of number of relevant zero minors
    # and the reciprocal of the number of columns
    base_matrix=matrix_from_system(system)
    base_score=scoring_function(base_matrix)
    for i in 1:(length(system)-1)
        max_score=base_score
        # Here we deal *only* with the exponent vectors
        max_monomial=[0 for i in system]
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
                    max_monomial=v2-v1
                end
            end
        end
        push!(monomialTranslations,max_monomial)
    end
    push!(monomialTranslations,[0 for i in system])
    montrans=[sum(monomialTranslations[i:end]) for i in 1:length(system)]
    lowest_term=abs(min([sort(m)[1] for m in montrans]...))
    var=gens(polRing)
    trans=[(Int(i!=j)*((var.^[lowest_term for j in var])))+(Int(i==j)*(var.^([lowest_term for j in var]+montrans[i]))) for j in 1:length(system) for i in 1:length(system)]
    for m in trans
        system=system.*m
    end
    return system
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

# returns number of zero minors that we consider
function number_of_zero_minors(mat)
    # we tag each column so that we can consider columns that appear more than once
    cols=[[i,mat[:,i]] for i in 1:number_of_columns(mat)]
    cols_per_sys=[filter(m->!iszero(m[2][i]),cols) for i in 1:number_of_rows(mat)]
    square_dim=min(number_of_columns(mat),number_of_rows(mat))
    # We assume the number of columns is larger than the number of rows for now, TODO
    minors=collect(Iterators.product(cols_per_sys...))
    # delete those with repeat columns since they will always be det 0
    minors=filter(m->length(m)==length(unique(m)),minors)
    minors=[Set(m) for m in minors]
    minors=unique(minors)
    minors=[[c[2] for c in m] for m in minors]
    # replace det with G-J to find zero entry in 
    niceminors=filter(m->is_det_zero(matrix(QQ,hcat(m...))),minors)
    # return the number of relevant zero minors
    return length(niceminors)
end

function number_of_fully_supported_minors(mat)
    # we tag each column so that we can consider columns that appear more than once
    cols=[[i,mat[:,i]] for i in 1:number_of_columns(mat)]
    cols_per_sys=[filter(m->!iszero(m[2][i]),cols) for i in 1:number_of_rows(mat)]
    square_dim=min(number_of_columns(mat),number_of_rows(mat))
    # We assume the number of columns is larger than the number of rows for now, TODO
    minors=collect(Iterators.product(cols_per_sys...))
    # delete those with repeat columns since they will always be det 0
    minors=filter(m->length(m)==length(unique(m)),minors)
    minors=[Set(m) for m in minors]
    minors=unique(minors)
    minors=[[c[2] for c in m] for m in minors]
    return length(minors)
end
