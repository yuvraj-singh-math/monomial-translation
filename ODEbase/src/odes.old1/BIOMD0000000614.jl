# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4)=rational_function_field(QQ,["k1","k2","k3","k4"])
polRing,(x1)=polynomial_ring(paramsRing,["x1"])
chemSystem=[(1*k4*k1 + 1*k4*k2*k3*x1 + (-1)*k4*k1*x1 + (-1)*k4*k2*k3*x1*x1)/k4 ]
constraints=[ ]
name="BIOMD0000000614"
pol=true
mass=false
rev=4
irr=0
def=0
rat=true
desc="Kamihira2000"
stoichMatrix=[ [ 1,1,-1,-1 ] ]
reconStoichMatrix=[ [ 1,-1,1,-1,-1,1,-1,1 ] ]
kineticMatrix=[ [ 0,1,0,1,1,0,1,0 ] ]
paramValues=[ 279/100000000 229/100 117/2000000 1 ]
paramNames=[ "k1" "k2" "a" "compartment_" ]
speciesNames=[ "f" ]
