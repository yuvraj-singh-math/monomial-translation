# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k42","k43","k44","k45"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22","x23","x24"])
chemSystem=[0, 0, 0, (1*k31*x23 + (-1)*k28*x4)/k32, 0, (1*(k23 + k22*x23) + (-1)*k20*x6)/k32, (1*k19*x6 + (-1)*k21*x7)/k32, 0, ((-1)*k11*x9 + 1*(k5*x11 + k4*x11*x7))/k32, 0, (1*k1*x12 + (-1)*k6*x11 + 1*k8*x15 + (-1)*k9*x11*x14 + 1*k10*x17 + (-1)*(k5*x11 + k4*x11*x7) + (-1)*k7*x11*x21)/k32, ((-1)*k1*x12 + (-1)*k2*x12 + 1*k3)/k32, 0, ((-1)*k12*x14 + (-1)*k9*x11*x14 + 1*(k13*x22*x21 - k14*x14) + 1*(k29*x20 - k30*x14))/k32, ((-1)*k8*x15 + 1*k7*x11*x21)/k32, 0, (1*k9*x11*x14 + (-1)*k10*x17)/k32, 0, 0, ((-5)*(k29*x20 - k30*x14) + 1*(k15*x23*x24 - k16*x20))/k33, ((-1)*k24*x21 + 1*k25*x4 + (-1)*(k13*x22*x21 - k14*x14) + (-1)*k7*x11*x21 + (-1)*(k26*x21 - k27*x24))/k32, (1*k12*x14 + 1*k10*x17 + (-1)*(k13*x22*x21 - k14*x14) + (-1)*(k17*x22 - k18*x23))/k32, (5*(k17*x22 - k18*x23) + (-1)*(k15*x23*x24 - k16*x20))/k33, ((-1)*(k15*x23*x24 - k16*x20) + 5*(k26*x21 - k27*x24))/k33 ]
constraints=[ 25*x14 + 25*x17 + 5*x20 + 25*x22 + 5*x23 - k44 + k45 ]
name="BIOMD0000000226"
pol=true
mass=false
rev=0
irr=24
def=4
rat=true
desc="Radulescu2008_NFkB_hierarchy_M_14_25_28_Lipniacky"
stoichMatrix=[ [ -1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 1,-1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0,-1, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0,-1,0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 1, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 1, 0, 0,-1, 1,-1, 1, 0,0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0,-1, 0,-1, 0, 0, 0, 0, 0,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 1, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0,-1, 0, 0, 0,-1, 0, 0,0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0,-1, 0, 0, 0,0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 1,-1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,-5, 1, 0, 0, 0], [ 0, 0,-1,1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0,-1,-1, 0, 0, 0, 0, 0,-1], [ 0, 0, 0,0, 0, 0, 1, 0, 0, 0, 0, 1, 0,0, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 5, 0,-1, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 5 ] ]
reconStoichMatrix=[ [ -1, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 1,-1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0,-1, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0,-1,0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 1, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 1, 0, 0,-1, 1,-1, 1, 0,0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0,-1, 0,-1, 0, 0, 0, 0, 0,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 1, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0,-1, 0, 0, 0,-1, 0, 0,0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0,-1, 0, 0, 0,0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 1,-1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,-5, 1, 0, 0, 0], [ 0, 0,-1,1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0,-1,-1, 0, 0, 0, 0, 0,-1], [ 0, 0, 0,0, 0, 0, 1, 0, 0, 0, 0, 1, 0,0, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 5, 0,-1, 0, 0, 0], [ 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 5 ] ]
kineticMatrix=[ [ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0], [ 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0 ] ]
paramValues=[ 1/400 1/8000 1/400000 1/10 3/2000 1/8000 1/5 1/10 1 1/10 1/8000 1/50000 92/5 0 92/5 0 1/400 0 1/2 1/2500 3/10000 1/2000000 0 1/10000 1/2 1/1000 1/2000 1/2500 1/100 0 1/2000000 1 1 0 0 0 0 0 0 0 0 0 0 3/2 0 ]
paramNames=[ "k1" "k2" "k3" "k4" "k5" "k6" "k7" "k8" "k9" "k10" "k11" "k12" "kf13" "kr13" "kf14" "kr14" "kf15" "kr15" "k16" "k17" "k18" "k20" "k19" "k21" "k22" "kf23" "kr23" "k27" "kf28" "kr28" "k26" "default" "c2" "s121" "s122" "s124" "s126" "s129" "s131" "s134" "s150" "s153" "s154" "__cm_k44" "__cm_k45" ]
speciesNames=[ "s121" "s122" "s124" "s125" "s126" "s127" "s128" "s129" "s130" "s131" "s132" "s133" "s134" "s135" "s139" "s150" "s152" "s153" "s154" "s159" "s160" "s161" "s164" "s167" ]