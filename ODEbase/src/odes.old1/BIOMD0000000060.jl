# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13"])
polRing,(x1,x2,x3,x4)=polynomial_ring(paramsRing,["x1","x2","x3","x4"])
chemSystem=[1*(k3*x3 - k4*k5^k6*x1)/k2, 1*(k7*k8^k9*x3 - k10*x2)/k2, ((-1)*(k3*x3 - k4*k5^k6*x1) + (-1)*(k7*k8^k9*x3 - k10*x2) + (-1)*(k11*x3 - k12*x4))/k2, 1*(k11*x3 - k12*x4)/k2 ]
constraints=[ 1*x1 + 1*x2 + 1*x3 + 1*x4 - k13 ]
name="BIOMD0000000060"
pol=true
mass=true
rev=3
irr=0
def=0
rat=true
desc="Keizer1996_Ryanodine_receptor_adaptation"
stoichMatrix=[ [ 1, 0, 0], [ 0, 1, 0], [ -1,-1,-1], [ 0, 0, 1 ] ]
reconStoichMatrix=[ [ 1,-1, 0, 0, 0, 0], [ 0, 0, 1,-1, 0, 0], [ -1, 1,-1, 1,-1, 1], [ 0, 0, 0, 0, 1,-1 ] ]
kineticMatrix=[ [ 0,1,0,0,0,0], [ 0,0,0,1,0,0], [ 1,0,1,0,1,0], [ 0,0,0,0,0,1 ] ]
paramValues=[ 0 1 144/5 1500 9/10 4 1500 9/10 3 3859/10 7/4 1/10 1 ]
paramNames=[ "Open_probability" "compartment" "__lp_r2_ka_minus" "__lp_r2_ka_plus" "__lp_r2_Ca" "__lp_r2_n" "__lp_r3_kb_plus" "__lp_r3_Ca" "__lp_r3_m" "__lp_r3_kb_minus" "__lp_r4_kc_plus" "__lp_r4_kc_minus" "__cm_k13" ]
speciesNames=[ "Pc1" "Po2" "Po1" "Pc2" ]
