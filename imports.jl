using Oscar;
chems=readdir("odebase/src/odes",join=true)
# these filters are mostly hacky workarounds
chems=filter(filename->occursin(".jl",filename)&&(!occursin("00553",filename))&&(!occursin("00552",filename))&&(!occursin("00346",filename))&&(!occursin("0091",filename))&&(!occursin("00108",filename))&&(!occursin("odebase.jl",filename))&&(!occursin("rejects.jl",filename))&&(!occursin("matrix",filename))&&!occursin("#",filename),chems)

struct OdebaseNode
    ID::String
    description::String
    numSpecies::Int
    numParams::Int
    speciesNames::Vector{String}
    paramNames::Vector{String}
    numIrr::Int
    numRev::Int
    deficiency::Int
    rational::Bool
    polynomial::Bool
    massAction::Bool
    paramsRing
    polRing
    ODEs::Vector
    constraints::Vector
    paramValues::Vector
    stoichMatrix::QQMatrix
    reconStoichMatrix::QQMatrix
    kineticMatrix::QQMatrix
    generic_polynomial_system::Vector
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
    println(ID)
    randCoeff=rand_nonzero(length(gens(paramsRing)));
    QQpolRing,tup=polynomial_ring(QQ,["$x" for x in gens(polRing)]);
    phi=hom(polRing,QQpolRing,c->evaluate(c,randCoeff),gens(QQpolRing));
    # we redefine polRing to be of rational type after the map
     push!(odebaseSystems,OdebaseNode(ID,desc,length(Oscar.gens(polRing)),length(gens(paramsRing)),speciesNames,paramNames,irr,rev,def,rat,pol,mass_bool,paramsRing,QQpolRing,chemSystem,constraints,[],matrix(QQ,[[]]),matrix(QQ,[[]]),matrix(QQ,[[]]),[phi(x) for x in union(chemSystem,constraints)]));
end

unfiltered_systems=[OdebaseNode(sys.ID,sys.description,sys.numSpecies,sys.numParams,sys.speciesNames,sys.paramNames,sys.numIrr,sys.numRev,sys.deficiency,sys.rational,sys.polynomial,sys.massAction,sys.paramsRing,sys.polRing,sys.ODEs,sys.constraints,sys.paramValues,sys.stoichMatrix,sys.reconStoichMatrix,sys.kineticMatrix,filter(x->!is_zero(x),unique(sys.generic_polynomial_system))) for sys in odebaseSystems];

for sys in odebaseSystems
    # Note that for f=2*x1, even though this is monomial, and has no toric solutions, is_monomial returns false
 #    so we look at the length of the list of monomials, check if its 1
    if sum([length(collect(monomials(f)))==1 for f in sys.generic_polynomial_system])>0
        rejects[sys.ID]="Contains monomial equation";
        filter!(s->s.ID!==sys.ID,unfiltered_systems);
    end

    #if length(unique(collect(Iterators.flatten([collect(monomials(f)) for f in sys.generic_polynomial_system]))))<length(sys.generic_polynomial_system)
        #rejects[sys.ID]="Macaulay matrix has more rows than columns";
        #filter!(s->s.ID!=sys.ID,unfiltered_systems);
    #end
end

unfiltered_systems=sort(unfiltered_systems,by= x->x.numSpecies);
systems=copy(unfiltered_systems);

