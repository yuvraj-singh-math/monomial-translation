# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"])
chemSystem=[0/k14, 0/k14, 0/k14, k11*k3*x5*x6 - k12*k4*x1*x6 - k3*x1 - k9*x1*x7^(k13 + 1) + k4*x4 + k10*x4*x7, k6*x7*x6 - k8*x2 + k9*x1*x7^(k13 + 1) + k10*x4*x7, k5*x7^k13*x5 - k7*x3 + k9*x1*x7^(k13 + 1), k12*k4*x1*x6 - k4*x4 - k10*x4*x7 ]
constraints=[ x5 - x8 - x1 - x3 - x4, x6 - x9 - x1 - x2 - 2*x4, x7 - x10 - x2 - k13*x3 ]
name="BIOMD0000000291"
pol=true
mass=false
rev=0
irr=0
def=0
rat=true
desc="Nikolaev2005_AlbuminBilirubinAdsorption"
stoichMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
reconStoichMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
kineticMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
paramValues=[ 121/250 3979/50000 1019/200000000 83/3125000 5489/1000000 1613/5000000000 301/100000 1011/10000000000 337/20000 53/400 95000 3000 1 1 ]
paramNames=[ "k1" "k2" "k3" "k4" "k5" "k6" "k7" "k8" "k9" "k10" "K_AlB" "K_AlB2" "n" "compartment" ]
speciesNames=[ "x1" "x2" "x3" "x4" "x5" "x6" "x7" "A0" "B0" "C0" ]
