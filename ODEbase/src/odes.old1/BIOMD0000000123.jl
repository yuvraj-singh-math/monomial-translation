# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k26,k25,k28,k27,k30,k29)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k26","k25","k28","k27","k30","k29"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14"])
chemSystem=[((-3)*k24*(k19*x8*x1^3 - k20*x4) + 1*(k23*k21*x2 - k24*k22*x1))/k24, ((-3)*k23*(k19*x14*x2^3 - k20*x10) + (-1)*(k23*k21*x2 - k24*k22*x1))/k23, (1*k24*(k1*x5 - k2*x3) + (-1)*k24*(k3*x3*x4 - k4*x6) + (-1)*(k24*k5*x3 - k23*k6*x9))/k24, ((-1)*k24*(k3*x3*x4 - k4*x6) + (-1)*(k24*k7*x4 - k23*k8*x10) + 1*k24*(k11*x7 - k12*x5*x4) + 1*k24*(k19*x8*x1^3 - k20*x4))/k24, ((-1)*k24*(k1*x5 - k2*x3) + 1*k24*(k11*x7 - k12*x5*x4) + 1*(k23*k15*x11 - k24*k16*x5))/k24, (1*k24*(k3*x3*x4 - k4*x6) + (-1)*k24*(k9*x6 - k10*x7) + (-1)*(k24*k13*x6 - k23*k14*x12))/k24, (1*k24*(k9*x6 - k10*x7) + (-1)*k24*(k11*x7 - k12*x5*x4) + 1*(k23*k17*x13 - k24*k18*x7))/k24, ((-1)*k24*(k19*x8*x1^3 - k20*x4) + 1*(k23*k8*x14 - k24*k7*x8))/k24, (1*(k24*k5*x3 - k23*k6*x9) + 1*k23*(k4*x12 - k3*x9*x10) + 1*k23*(k1*x11 - k2*x9))/k23, (1*(k24*k7*x4 - k23*k8*x10) + 1*k23*(k11*x13 - k12*x11*x10) + 1*k23*(k4*x12 - k3*x9*x10) + 1*k23*(k19*x14*x2^3 - k20*x10))/k23, (1*k23*(k11*x13 - k12*x11*x10) + (-1)*(k23*k15*x11 - k24*k16*x5) + (-1)*k23*(k1*x11 - k2*x9))/k23, (1*(k24*k13*x6 - k23*k14*x12) + (-1)*k23*(k9*x12 - k10*x13) + (-1)*k23*(k4*x12 - k3*x9*x10))/k23, (1*k23*(k9*x12 - k10*x13) + (-1)*k23*(k11*x13 - k12*x11*x10) + (-1)*(k23*k17*x13 - k24*k18*x7))/k23, ((-1)*k23*(k19*x14*x2^3 - k20*x10) + (-1)*(k23*k8*x14 - k24*k7*x8))/k23 ]
constraints=[ 3*k24*x1 + 3*k23*x2 + 9*k24*x4 + 9*k24*x6 + 9*k24*x7 + 9*k23*x10 + 9*k23*x12 + 9*k23*x13 - k23*k26 + k24*k25, 1*k24*x3 + 1*k24*x5 + 1*k24*x6 + 1*k24*x7 + 1*k23*x9 + 1*k23*x11 + 1*k23*x12 + 1*k23*x13 - k23*k28 + k24*k27, 1*k24*x4 + 1*k24*x6 + 1*k24*x7 + 1*k24*x8 + 1*k23*x10 + 1*k23*x12 + 1*k23*x13 + 1*k23*x14 - k23*k30 + k24*k29 ]
name="BIOMD0000000123"
pol=true
mass=false
rev=17
irr=0
def=3
rat=true
desc="Fisher2006_NFAT_Activation"
stoichMatrix=[ [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0,-1], [ 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0,-1, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [ -1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [ 0, 1, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0], [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0], [ 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0,-1, 0 ] ]
reconStoichMatrix=[ [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 1,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0,-1, 1], [ 1,-1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0,-1, 1, 0, 0,-1, 1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], [ -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 0, 0], [ 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0,-1, 1, 0, 0 ] ]
kineticMatrix=[ [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,1], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,1,0], [ 0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 1,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0], [ 1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0], [ 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0 ] ]
paramValues=[ 2/78125 8/3125 663/100 21/12500 3/3125 3/2000 23/25000 19/10000 8/3125 1/2 21/12500 663/100 1/200 1/2 1/200 1/2 1/200 1/2 1 1 21/100 1/2 269/1000000000000000 113/1000000000000000 15000783/5000000 30090063/10000000 95581/10000000 16993/10000000 48641/5000000 501987/10000000 ]
paramNames=[ "k1" "k2" "k16" "k15" "k18" "k17" "k6" "k5" "k14" "k13" "k12" "k11" "k10" "k9" "k3" "k4" "k7" "k8" "k19" "k20" "k21" "k22" "cytosol" "nucleus" "__cm_k26" "__cm_k25" "__cm_k28" "__cm_k27" "__cm_k30" "__cm_k29" ]
speciesNames=[ "Ca_Nuc" "Ca_Cyt" "NFAT_Nuc" "Act_C_Nuc" "NFAT_Pi_Nuc" "NFAT_Act_C_Nuc" "NFAT_Pi_Act_C_Nuc" "Inact_C_Nuc" "NFAT_Cyt" "Act_C_Cyt" "NFAT_Pi_Cyt" "NFAT_Act_C_Cyt" "NFAT_Pi_Act_C_Cyt" "Inact_C_Cyt" ]