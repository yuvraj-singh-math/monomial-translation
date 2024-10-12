# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8"])
chemSystem=[(1*k9*k3*x3 + (-1)*k9*k1*x1 + (-2)*k9*k6*x1*(x1 - 1)/2*x4 + 2*k9*k8*x6)/k9, (1*k9*k4*x4 + (-1)*k9*k2*x2 + (-2)*k9*k5*x2*(x2 - 1)/2*x3 + 2*k9*k7*x5)/k9, ((-1)*k9*k5*x2*(x2 - 1)/2*x3 + 1*k9*k7*x5)/k9, ((-1)*k9*k6*x1*(x1 - 1)/2*x4 + 1*k9*k8*x6)/k9, (1*k9*k5*x2*(x2 - 1)/2*x3 + (-1)*k9*k7*x5)/k9, (1*k9*k6*x1*(x1 - 1)/2*x4 + (-1)*k9*k8*x6)/k9, ((-1)*k9*k3*x3 + 1*k9*k1*x1)/k9, ((-1)*k9*k4*x4 + 1*k9*k2*x2)/k9 ]
constraints=[ 1*x1 + 2*x6 + 1*x7 - k10, 1*x2 + 2*x5 + 1*x8 - k11, 1*x3 + 1*x5 - k12, 1*x4 + 1*x6 - k13 ]
name="BIOMD0000000483"
pol=true
mass=false
rev=0
irr=8
def=0
rat=true
desc="Cao2008 - Network of a toggle switch"
stoichMatrix=[ [ 1, 0,-1, 0,-2, 0, 0, 2], [ 0, 1, 0,-1, 0,-2, 2, 0], [ 0, 0, 0, 0, 0,-1, 1, 0], [ 0, 0, 0, 0,-1, 0, 0, 1], [ 0, 0, 0, 0, 0, 1,-1, 0], [ 0, 0, 0, 0, 1, 0, 0,-1], [ -1, 0, 1, 0, 0, 0, 0, 0], [ 0,-1, 0, 1, 0, 0, 0, 0 ] ]
reconStoichMatrix=[ [ 1, 0,-1, 0,-2, 0, 0, 2], [ 0, 1, 0,-1, 0,-2, 2, 0], [ 0, 0, 0, 0, 0,-1, 1, 0], [ 0, 0, 0, 0,-1, 0, 0, 1], [ 0, 0, 0, 0, 0, 1,-1, 0], [ 0, 0, 0, 0, 1, 0, 0,-1], [ -1, 0, 1, 0, 0, 0, 0, 0], [ 0,-1, 0, 1, 0, 0, 0, 0 ] ]
kineticMatrix=[ [ 0,0,1,0,2,0,0,0], [ 0,0,0,1,0,2,0,0], [ 0,0,0,0,0,1,0,0], [ 0,0,0,0,1,0,0,0], [ 0,0,0,0,0,0,1,0], [ 0,0,0,0,0,0,0,1], [ 1,0,0,0,0,0,0,0], [ 0,1,0,0,0,0,0,0 ] ]
paramValues=[ 1 1 100 100 1/100000 1/100000 1/10 1/10 1 0 0 0 0 ]
paramNames=[ "da" "db" "sa" "sb" "ba" "bb" "ua" "ub" "default" "__cm_k10" "__cm_k11" "__cm_k12" "__cm_k13" ]
speciesNames=[ "Pa" "Pb" "Da" "Db" "BDa" "BDb" "ESA" "ESB" ]
