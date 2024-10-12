# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15"])
polRing,(x1,x2,x3,x4,x5,x6)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6"])
chemSystem=[((-1)*k4*k5*x3*x1 + (-1)*k4*k7*x4*x1 + (-1)*k4*k2*x1 + 1*k4*k12)/k4, (1*k4*k5*x3*x1 + 1*k4*k7*x4*x1 + (-1)*k4*k1*x2 + (-1)*k4*(k8*x2*k14 - k9*k15) + 1*k4*k10*k15)/k4, ((-1)*k4*k6*x2^k3*x3 + (-1)*k4*k1*x3 + 1*k4*k11)/k4, (1*k4*k6*x2^k3*x3 + (-1)*k4*k1*x4)/k4, 0, 0 ]
constraints=[ ]
name="BIOMD0000000630"
pol=true
mass=false
rev=1
irr=13
def=4
rat=true
desc="Venkatraman2011 - PLS-UPA behaviour in the presence of substrate competition_1_1_1_1"
stoichMatrix=[ [ -1, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0,0,1,0], [ 1, 0, 1, 0, 0,-1, 0,-1, 1, 0, 0,0,0,0], [ 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0,1,0,0], [ 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,0,0,0], [ 0, 0, 0, 0, 0, 0, 0,-1, 0,-1, 0,0,0,1], [ 0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1,0,0,0 ] ]
reconStoichMatrix=[ [ -1, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0,0,1,0], [ 1, 0, 1, 0, 0,-1, 0,-1, 1, 1, 0, 0,0,0,0], [ 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0,1,0,0], [ 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,0,0,0], [ 0, 0, 0, 0, 0, 0, 0,-1, 1, 0,-1, 0,0,0,1], [ 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0,-1,0,0,0 ] ]
kineticMatrix=[ [ 1,0,1,0,1,0,0,0,0,0,0,0,0,0,0], [ 0,1,0,0,0,1,0,1,0,0,0,0,0,0,0], [ 1,1,0,1,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,1,0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,1,0,0,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,1,0,1,0,0,0 ] ]
paramValues=[ 21/250 4/125 2 1 7/200 40 9/10 0 2/125 1/50 2/625 1/100 1/100 0 0 ]
paramNames=[ "parameter_1" "parameter_2" "parameter_13" "compartment_1" "__lp_r2_k1" "__lp_r3_parameter_8" "__lp_r4_k1" "__lp_r9_k1" "__lp_r9_k2" "__lp_r10_k1" "__lp_r13_v" "__lp_r14_v" "__lp_r15_v" "species_5" "species_6" ]
speciesNames=[ "species_1" "species_2" "species_3" "species_4" "species_5" "species_6" ]
