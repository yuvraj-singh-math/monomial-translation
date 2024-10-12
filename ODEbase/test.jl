using Oscar;

chems=readdir("/home/main/notes/odebase")
chems=filter(filename->occursin(".jl",filename)&&(filename!="test.jl"),chems)

# figure out way to select by no. of species, &c.
# compute degenericity of matrices

struct OdebaseNode
    ID::String
    rational::Bool
    massAction::Bool
    species::Int
    defecit::Int
    numParams::Int
    # TODO types for the following
    # \dot{x}_i is set to 0
    param_polynomial_system
    constraints
    paramsRing
end

for file in chems

end

