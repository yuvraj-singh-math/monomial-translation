# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15"])
polRing,(x1,x2,x3,x4,x5,x6,x7)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7"])
chemSystem=[(1*k1*x2 + (-1)*k2*x1*x3)/k15, (1*k13*x7 + (-1)*k14*x2)/k15, (1*k3*x4 + (-1)*k4*x3)/k15, (1*k9*x1 + (-1)*k10*x6*x4)/k15, (1*k5*x2 + (-1)*k6*x3*x5)/k15, (1*k7 + (-1)*k8*x5*x6)/k15, (1*k11*x1 + (-1)*k12*x7)/k15 ]
constraints=[ ]
ID="BIOMD0000000229"
pol=true
mass_bool=false
rev=14
irr=0
def=0
rat=true
desc="Ma2002_cAMP_oscillations"
stoichMatrix=[ [ 1,-1,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0], [ 0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,1,-1], [ 0, 0,1,-1,0, 0,0, 0,0, 0,0, 0,0, 0], [ 0, 0,0, 0,0, 0,0, 0,1,-1,0, 0,0, 0], [ 0, 0,0, 0,1,-1,0, 0,0, 0,0, 0,0, 0], [ 0, 0,0, 0,0, 0,1,-1,0, 0,0, 0,0, 0], [ 0, 0,0, 0,0, 0,0, 0,0, 0,1,-1,0, 0 ] ]
reconStoichMatrix=[ [ 1,-1,-1,1,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0], [ 0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,1,-1,-1,1], [ 0, 0, 0,0,1,-1,-1,1,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0], [ 0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,1,-1,-1,1,0, 0, 0,0,0, 0, 0,0], [ 0, 0, 0,0,0, 0, 0,0,1,-1,-1,1,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0], [ 0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,1,-1,-1,1,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0], [ 0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,0, 0, 0,0,1,-1,-1,1,0, 0, 0,0 ] ]
kineticMatrix=[ [ 0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0], [ 0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0 ] ]
paramValues=[ 2,9//10,5//2,3//2,3//5,4//5,1,13//10,3//10,4//5,7//10,49//10,23,9//2,1 ]
paramNames=[ "k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","compartment" ]
speciesNames=[ "ACA","CAR1","PKA","incAMP","ERK2","REGA","excAMP" ]
