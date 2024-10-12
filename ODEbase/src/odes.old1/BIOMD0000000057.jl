# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35"])
polRing,(x1,x2,x3,x4,x5,x6)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6"])
chemSystem=[((-1)*k28*((k10*k8 + k11*k7)/(k8 + k7*(1 + k8/k4))*k29*x1 - (k14 + k16*k7)/(1 + k7/k17)*x2) + (-1)*k28*((k3*k4 + k5)*k7/(k4 + k7*(1 + k4/k8))*x1 - (k30 + k31)*x3))/k28, (1*k28*((k10*k8 + k11*k7)/(k8 + k7*(1 + k8/k4))*k29*x1 - (k14 + k16*k7)/(1 + k7/k17)*x2) + (-1)*k28*(k19*k17/(k17 + k7)*x2 - k32*x4) + (-1)*k28*((k21*k17 + k22)*k7/(k17 + k7)*x2 - k4*(k24 + k25)/(k4 + k7)*x5))/k28, 1*k28*((k3*k4 + k5)*k7/(k4 + k7*(1 + k4/k8))*x1 - (k30 + k31)*x3)/k28, 1*k28*(k19*k17/(k17 + k7)*x2 - k32*x4)/k28, (1*k28*((k21*k17 + k22)*k7/(k17 + k7)*x2 - k4*(k24 + k25)/(k4 + k7)*x5) + (-1)*k28*((k3*k4 + k5)*k7/(k4 + k7)*x5 - (k33 + k34)*x6))/k28, 1*k28*((k3*k4 + k5)*k7/(k4 + k7)*x5 - (k33 + k34)*x6)/k28 ]
constraints=[ 1*k28*x1 + 1*k28*x2 + 1*k28*x3 + 1*k28*x4 + 1*k28*x5 + 1*k28*x6 - k28*k35 ]
name="BIOMD0000000057"
pol=true
mass=true
rev=5
irr=0
def=0
rat=true
desc="Sneyd2002_IP3_Receptor"
stoichMatrix=[ [ -1,-1, 0, 0, 0], [ 1, 0,-1,-1, 0], [ 0, 1, 0, 0, 0], [ 0, 0, 1, 0, 0], [ 0, 0, 0, 1,-1], [ 0, 0, 0, 0, 1 ] ]
reconStoichMatrix=[ [ -1, 1,-1, 1, 0, 0, 0, 0, 0, 0], [ 1,-1, 0, 0,-1, 1,-1, 1, 0, 0], [ 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 1,-1,-1, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1 ] ]
kineticMatrix=[ [ 1,0,1,0,0,0,0,0,0,0], [ 0,1,0,0,1,0,1,0,0,0], [ 0,0,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,1,0,0,0,0], [ 0,0,0,0,0,0,0,1,1,0], [ 0,0,0,0,0,0,0,0,0,1 ] ]
paramValues=[ 0 2221/7265 16/25 3/25 17/10 4/5 10 1/40 10761/7265 187/5 17/10 72204/3235 1/25 7/5 149/5 5/2 547/10 6017/64700 11/100 492580/647 4 4707 1791/12650 27/50 57/5 2221/1265 10 1 10 1/25 4/5 149/5 1/25 4/5 1 ]
paramNames=[ "open_probability" "Phi1" "k1" "L1" "l2" "lminus2" "c" "L3" "Phi2" "k2" "l4" "Phi_minus2" "kminus1" "kminus2" "kminus3" "lminus4" "L5" "Phi3" "k3" "Phi4" "k4" "l6" "Phi_minus4" "kminus4" "lminus6" "Phi5" "IP3" "compartment" "__lp_r2_IP3" "__lp_r3_kminus1" "__lp_r3_lminus2" "__lp_r4_kminus3" "__lp_r6_kminus1" "__lp_r6_lminus2" "__cm_k35" ]
speciesNames=[ "R" "O" "I1" "S" "A" "I2" ]
