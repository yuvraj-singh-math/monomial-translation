# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53,k54,k55,k56)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k42","k43","k44","k45","k46","k47","k48","k49","k50","k51","k52","k53","k54","k55","k56"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12"])
chemSystem=[((-1)*k53*(k29/(k6/(83/50))*x1*x3 - k29*x6) + (-1)*k53*(k1*k21/k2*x1 - k21*x2) + (-1)*k53*(k29/(k7*k6/(83/50))*x1*x4 - k29*x11) + (-1)*k53*(k29/(k8*k6/(83/50))*x1*x5 - k29*x10))/k53, (1*k53*(k1*k21/k2*x1 - k21*x2) + (-1)*k53*(k29/(k3*k6/(83/50))*x2*x3 - k29*x8) + (-1)*k53*(k29/(k4*k7*k6/(83/50))*x2*x4 - k29*x9) + (-1)*k53*(k29/(k5*k8*k6/(83/50))*x2*x5 - k29*x12))/k53, ((-1)*k53*(k29/(k6/(83/50))*x1*x3 - k29*x6) + (-1)*k53*(k29/(k3*k6/(83/50))*x2*x3 - k29*x8) + (-1)*k53*k52*x3 + 1*k53*(83/50*13/100*x5*x7 - 83/50*13/100/(83/50)*k51*x3))/k53, ((-1)*k53*(k29/(k7*k6/(83/50))*x1*x4 - k29*x11) + (-1)*k53*(k29/(k4*k7*k6/(83/50))*x2*x4 - k29*x9) + 1*k53*k52*x3 + (-1)*k53*k10*x4 + 1*k53*k48*x5)/k53, ((-1)*k53*(k29/(k8*k6/(83/50))*x1*x5 - k29*x10) + (-1)*k53*(k29/(k5*k8*k6/(83/50))*x2*x5 - k29*x12) + 1*k53*k10*x4 + (-1)*k53*(83/50*13/100*x5*x7 - 83/50*13/100/(83/50)*k51*x3) + (-1)*k53*k48*x5)/k53, (1*k53*(k29/(k6/(83/50))*x1*x3 - k29*x6) + (-1)*k53*(k1*k21/(k3*k2)*x6 - k21*x8) + (-1)*k53*k12*x6 + 1*k53*(83/50*13/100*x10*x7 - 83/50*13/100/(83/50)*k51/k8*x6))/k53, (1*k53*k52*x3 + (-1)*k53*(83/50*13/100*x5*x7 - 83/50*13/100/(83/50)*k51*x3) + 1*k53*k12*x6 + 1*k53*k13*x8 + (-1)*k53*(83/50*13/100*x10*x7 - 83/50*13/100/(83/50)*k51/k8*x6) + (-1)*k53*(83/50*13/100*x12*x7 - k3*83/50*13/100/(83/50)*k51/(k5*k8)*x8))/k53, (1*k53*(k29/(k3*k6/(83/50))*x2*x3 - k29*x8) + 1*k53*(k1*k21/(k3*k2)*x6 - k21*x8) + (-1)*k53*k13*x8 + 1*k53*(83/50*13/100*x12*x7 - k3*83/50*13/100/(83/50)*k51/(k5*k8)*x8))/k53, (1*k53*(k29/(k4*k7*k6/(83/50))*x2*x4 - k29*x9) + 1*k53*(k1*k21/(k4*k2)*x11 - k21*x9) + 1*k53*k13*x8 + (-1)*k53*k16*x9 + 1*k53*k13*x12)/k53, (1*k53*(k29/(k8*k6/(83/50))*x1*x5 - k29*x10) + 1*k53*k16*x11 + (-1)*k53*(83/50*13/100*x10*x7 - 83/50*13/100/(83/50)*k51/k8*x6) + (-1)*k53*(k1*k21/(k5*k2)*x10 - k21*x12) + (-1)*k53*k12*x10)/k53, (1*k53*(k29/(k7*k6/(83/50))*x1*x4 - k29*x11) + (-1)*k53*(k1*k21/(k4*k2)*x11 - k21*x9) + 1*k53*k12*x6 + (-1)*k53*k16*x11 + 1*k53*k12*x10)/k53, (1*k53*(k29/(k5*k8*k6/(83/50))*x2*x5 - k29*x12) + 1*k53*k16*x9 + (-1)*k53*(83/50*13/100*x12*x7 - k3*83/50*13/100/(83/50)*k51/(k5*k8)*x8) + 1*k53*(k1*k21/(k5*k2)*x10 - k21*x12) + (-1)*k53*k13*x12)/k53 ]
constraints=[ 1*k53*x1 + 1*k53*x2 + 1*k53*x6 + 1*k53*x8 + 1*k53*x9 + 1*k53*x10 + 1*k53*x11 + 1*k53*x12 - k53*k54, 1*k53*x3 + 1*k53*x4 + 1*k53*x5 + 1*k53*x6 + 1*k53*x8 + 1*k53*x9 + 1*k53*x10 + 1*k53*x11 + 1*k53*x12 - k53*k55, 1*k53*x3 + 1*k53*x6 + 1*k53*x7 + 1*k53*x8 - k53*k56 ]
name="BIOMD0000000637"
pol=true
mass=true
rev=13
irr=9
def=9
rat=true
desc="Bush2016 - Simplified Carrousel model of GPCR"
stoichMatrix=[ [ -1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 1, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ -1, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0,-1, 0, 0,-1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0,-1, 0, 0,-1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0], [ 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 1, 0, 0,-1,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 1], [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0,-1, 0,-1, 0], [ 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 0, 0,-1 ] ]
reconStoichMatrix=[ [ -1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,-1, 0, 0, 0, 0,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ -1, 1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0], [ 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 1, 1, 0, 0,-1, 1,-1, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1], [ 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 1, 0, 0,-1, 1, 0,-1, 0], [ 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 1, 1,-1, 0, 0,-1 ] ]
kineticMatrix=[ [ 1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0], [ 0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0], [ 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1 ] ]
paramValues=[ 0 28/5 1 1 1 36 1 1 31/50000 1/500 1079/5000 31/50000 3/2 11/100 1079/5000 11/100 1079/5000 3308 2042 0 1/1000 0 1/1000 0 1/1000 0 1/1000 83/18000 1/10 83/18000 1/10 83/18000 1/10 83/18000 1/10 83/18000 1/10 83/18000 1/10 13/10000 13/10000 13/10000 0 0 25525/49 0 0 31/50000 31/50000 3/2 1/100 31/50000 98/25 41350/49 25525/49 25525/49 ]
paramNames=[ "L" "K_d_L_R" "lambda" "lambda_t" "lambda_d" "K_d_R_G" "eta" "rho" "k_Ef_G" "k_Hf_Gt" "k_Af_Gd" "k_Ef_RG" "k_Ef_LRG" "k_Hf_LRGt" "k_Af_LRGd" "k_Hf_RGt" "k_Af_RGd" "Rtot" "Gtot" "k_on_L_R" "k_off_L_R" "k_on_L_RG" "k_off_L_RG" "k_on_L_RGt" "k_off_L_RGt" "k_on_L_RGd" "k_off_L_RGd" "k_on_R_G" "k_off_R_G" "k_on_LR_G" "k_off_LR_G" "k_on_R_Gt" "k_off_R_Gt" "k_on_LR_Gt" "k_off_LR_Gt" "k_on_R_Gd" "k_off_R_Gd" "k_on_LR_Gd" "k_off_LR_Gd" "k_Ar_Gd" "k_Ar_RGd" "k_Ar_LRGd" "tot_LR" "tot_RG" "tot_G" "tot_Gd" "tot_Gt" "k_Ef_Gd" "k_Ef_RGd" "k_Ef_LRGd" "K_d_Gd_Gbg" "ModelValue_47" "PM" "__cm_k54" "__cm_k55" "__cm_k56" ]
speciesNames=[ "R" "LR" "G" "Gt" "Gd" "RG" "Gbg" "LRG" "LRGt" "RGd" "RGt" "LRGd" ]
