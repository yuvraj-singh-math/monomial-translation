# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10"])
polRing,(x1,x2)=polynomial_ring(paramsRing,["x1","x2"])
chemSystem=[(1*k10*k1 + (-1)*k10*k4*x1 + 1*k10*k6*(1 - k2*x2)*x2*x1)/k10, ((-1)*k10*k3*x2*x1 + 1*k10*k5*x2*(1 - (x2/k7)^k9))/k10 ]
constraints=[ ]
name="BIOMD0000000919"
pol=true
mass=false
rev=0
irr=5
def=0
rat=true
desc="Ledzewicz2013 - On optimal chemotherapy with a strongly targeted agent for a model of tumor immune system interactions with generalized logistic growth"
stoichMatrix=[ [ 1,-1,1, 0,0], [ 0, 0,0,-1,1 ] ]
reconStoichMatrix=[ [ 1,-1,1, 0,0], [ 0, 0,0,-1,1 ] ]
kineticMatrix=[ [ 0,1,0,0,0], [ 0,0,0,1,0 ] ]
paramValues=[ 1181/10000 33/12500 1 37451/100000 5599/10000 121/25000 780 1 1 1 ]
paramNames=[ "alpha" "beta" "gamma" "delta" "mu_C" "mu_I" "x_inf" "kappa" "v" "compartment" ]
speciesNames=[ "y" "x" ]
