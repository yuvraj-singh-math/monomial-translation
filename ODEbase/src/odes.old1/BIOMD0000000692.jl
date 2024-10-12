# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8"])
chemSystem=[1*k3*(k7*x5 - k6*x1*x6)/k3/k3, (-1)*k3*(k13*x2*x8 - k12*x3)/k3/k3, ((-1)*k3*(k9*x3 - k8*x4)/k3 + 1*k3*(k13*x2*x8 - k12*x3)/k3)/k3, ((-1)*k3*(k5*x4 - k4*x5)/k3 + 1*k3*(k9*x3 - k8*x4)/k3)/k3, (1*k3*(k5*x4 - k4*x5)/k3 + (-1)*k3*(k7*x5 - k6*x1*x6)/k3)/k3, (1*k3*(k7*x5 - k6*x1*x6)/k3 + (-1)*k3*(k11*x6 - k10*x7*x8)/k3)/k3, 1*k3*(k11*x6 - k10*x7*x8)/k3/k3, (1*k3*(k11*x6 - k10*x7*x8)/k3 + (-1)*k3*(k13*x2*x8 - k12*x3)/k3)/k3 ]
constraints=[ 1*k3*x1 + 1*k3*x2 + 1*k3*x3 + 1*k3*x4 + 1*k3*x5 - k3*k14, 1*k3*x2 + 1*k3*x3 + 1*k3*x4 + 1*k3*x5 + 1*k3*x6 + 1*k3*x7 - k3*k15, 1*k3*x3 + 1*k3*x4 + 1*k3*x5 + 1*k3*x6 + 1*k3*x8 - k3*k16 ]
name="BIOMD0000000692"
pol=true
mass=true
rev=5
irr=0
def=0
rat=true
desc="Phillips2003 - The Mechanism of Ras GTPase Activation by Neurofibromin"
stoichMatrix=[ [ 0, 1, 0, 0, 0], [ 0, 0, 0, 0,-1], [ 0, 0,-1, 0, 1], [ -1, 0, 1, 0, 0], [ 1,-1, 0, 0, 0], [ 0, 1, 0,-1, 0], [ 0, 0, 0, 1, 0], [ 0, 0, 0, 1,-1 ] ]
reconStoichMatrix=[ [ 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0,-1, 1], [ 0, 0, 0, 0,-1, 1, 0, 0, 1,-1], [ -1, 1, 0, 0, 1,-1, 0, 0, 0, 0], [ 1,-1,-1, 1, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,-1, 0, 0,-1, 1, 0, 0], [ 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], [ 0, 0, 0, 0, 0, 0, 1,-1,-1, 1 ] ]
kineticMatrix=[ [ 0,0,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,0], [ 0,0,0,0,1,0,0,0,0,1], [ 1,0,0,0,0,1,0,0,0,0], [ 0,1,1,0,0,0,0,0,0,0], [ 0,0,0,1,0,0,1,0,0,0], [ 0,0,0,0,0,0,0,1,0,0], [ 0,0,0,0,0,0,0,1,1,0 ] ]
paramValues=[ 0 100 1309/2500000000000000 14399/125000000000000000 51051/5000000000000000 282741/5000000000000000000000 1309/62500000000000 14399/5000000000000000 43773/200000000000000 314159/500000000000000000 121737/5000000000000000 33301/10000000000000000 314159/500000000000000000 2 1 62496021135727/6250000000000 ]
paramNames=[ "Pi_curve" "hplc_curve" "geometry" "__lp_r2_kb" "__lp_r2_kf" "__lp_r3_kb" "__lp_r3_kf" "__lp_r4_kb" "__lp_r4_kf" "__lp_r5_kb" "__lp_r5_kf" "__lp_r6_kb" "__lp_r6_kf" "__cm_k14" "__cm_k15" "__cm_k16" ]
speciesNames=[ "Pi" "RasGTP" "RasGTP_minus_NF1" "RasGTP_minus_NF1_star_" "RasGDP_minus_NF1_Pi" "RasGDP_NF1" "RasGDP" "NF1" ]