# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17"])
polRing,(x1,x2,x3,x4)=polynomial_ring(paramsRing,["x1","x2","x3","x4"])
chemSystem=[(1*k17*k1*x3 + 1*k17*k7*x1*x3*(1 - (x1 + x2)/k5) + (-1)*k17*k13*x1)/k17, (1*k17*k2*x4 + 1*k17*k8*x2*x4*(1 - (x2 + x1)/k5) + (-1)*k17*k14*x2)/k17, (1*k17*k3*x1 + 1*k17*k9*x3*(1 - (x3 + x4)/k6) + (-1)*k17*k15*x3 + (-1)*k17*k12*x3 + 1*k17*k11*x4)/k17, (1*k17*k12*x3 + (-1)*k17*k11*x4 + 1*k17*k4*x2 + 1*k17*k10*x4*x2*(1 - (x4 + x3)/k6) + (-1)*k17*k16*x4)/k17 ]
constraints=[ ]
name="BIOMD0000000770"
pol=true
mass=false
rev=0
irr=14
def=0
rat=true
desc="Eftimie2017/1 - interaction of Th and macrophage"
stoichMatrix=[ [ 1,1,-1,0,0, 0,0,0, 0, 0, 0,0,0, 0], [ 0,0, 0,1,1,-1,0,0, 0, 0, 0,0,0, 0], [ 0,0, 0,0,0, 0,1,1,-1,-1, 1,0,0, 0], [ 0,0, 0,0,0, 0,0,0, 0, 1,-1,1,1,-1 ] ]
reconStoichMatrix=[ [ 1,1,-1,0,0, 0,0,0, 0, 0, 0,0,0, 0], [ 0,0, 0,1,1,-1,0,0, 0, 0, 0,0,0, 0], [ 0,0, 0,0,0, 0,1,1,-1,-1, 1,0,0, 0], [ 0,0, 0,0,0, 0,0,0, 0, 1,-1,1,1,-1 ] ]
kineticMatrix=[ [ 0,0,1,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,1,0,0,1 ] ]
paramValues=[ 1/125 1/1000 1/1000 1/1000 100000000 1000000000 9/100 9/100 1/50 1/50 9/100 1/20 3/100 3/100 1/50 1/50 1 ]
paramNames=[ "ah1" "ah2" "am1" "am2" "m1" "m2" "ph1" "ph2" "pm1" "pm2" "rm1" "rm2" "eh1" "eh2" "em1" "em2" "tme" ]
speciesNames=[ "H1" "H2" "M1" "M2" ]
