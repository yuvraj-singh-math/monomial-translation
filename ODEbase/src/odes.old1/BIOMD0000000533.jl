# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11"])
chemSystem=[(-k1)*(k7 + x8)*x1 - k2*x7*x1, k2*x7*x1 - k3*x7^k6*x2, k3*x7^k6*x2 - 4*k4*x3, 4*k4*x3 - 4*k4*x4, 4*k4*x4 - 4*k4*x5, 4*k4*x5 - 4*k4*x6, 4*k4*x6 - k5*x7, k5*x7, k1*(k7 + x8)*x1 ]
constraints=[ x10 - x1 + x9, x11 - x8 + x3 + x4 + x5 + x6 + x7 ]
name="BIOMD0000000533"
pol=true
mass=false
rev=0
irr=0
def=0
rat=true
desc="Steckmann2012 - Amyloid beta-protein fibrillogenesis (kinetics of secondary structure conversion)"
stoichMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
reconStoichMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
kineticMatrix=[ [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ], [ ] ]
paramValues=[ 59/100 84/125 339/500 49/1250 277/500 2 0 1 ]
paramNames=[ "k0" "k1" "k2" "k3" "k4" "q" "epsilon" "cell" ]
speciesNames=[ "RCT0" "alpha" "BN1" "BN2" "BN3" "BN4" "BTX" "BM" "RCT1" "RC" "beta" ]
