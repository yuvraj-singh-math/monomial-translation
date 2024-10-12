# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53,k54)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k42","k43","k44","k45","k46","k47","k48","k49","k50","k51","k52","k53","k54"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22"])
chemSystem=[(1*k3*(k6*x4 - k7*x3*x1) + 1*k3*(k12*x15 - k13*x10*x1) + 1*k3*(k16*x13 - k17*x8*x1) + 1*k3*(k18*x11 - k19*x9*x1) + (-1)*k3*k29*x1 + 1*k3*k38*x20*x2 + 1*k3*k44*x17*x2)/k3, (1*k3*k29*x1 + (-1)*k3*k38*x20*x2 + (-1)*k3*k44*x17*x2)/k3, (1*k3*(k4*x9 - k5*x8*x3) + 1*k3*(k6*x4 - k7*x3*x1) + 1*k3*(k20*x15 - k21*x11*x3) + 1*k3*(k22*x10 - k23*x9*x3) + 1*k3*(k24*x6 - k25*x3*x4) + 1*k3*(k26*x11 - k27*x13*x3))/k3, ((-1)*k3*(k6*x4 - k7*x3*x1) + 1*k3*(k8*x11 - k9*x8*x4) + 1*k3*(k10*x15 - k11*x9*x4) + 1*k3*(k24*x6 - k25*x3*x4) + (-1)*k3*k30*x4 + 1*k3*k39*x20*x5 + 1*k3*k45*x17*x5)/k3, (1*k3*k30*x4 + (-1)*k3*k39*x20*x5 + (-1)*k3*k45*x17*x5)/k3, (1*k3*(k14*x15 - k15*x8*x6) + (-1)*k3*(k24*x6 - k25*x3*x4) + (-1)*k3*k31*x6 + 1*k3*k40*x20*x7 + 1*k3*k46*x17*x7)/k3, (1*k3*k31*x6 + (-1)*k3*k40*x20*x7 + (-1)*k3*k46*x17*x7)/k3, (1*k3*(k4*x9 - k5*x8*x3) + 1*k3*(k8*x11 - k9*x8*x4) + 1*k3*(k14*x15 - k15*x8*x6) + 1*k3*(k16*x13 - k17*x8*x1))/k3, ((-1)*k3*(k4*x9 - k5*x8*x3) + 1*k3*(k10*x15 - k11*x9*x4) + 1*k3*(k18*x11 - k19*x9*x1) + 1*k3*(k22*x10 - k23*x9*x3))/k3, (1*k3*(k12*x15 - k13*x10*x1) + (-1)*k3*(k22*x10 - k23*x9*x3))/k3, ((-1)*k3*(k8*x11 - k9*x8*x4) + (-1)*k3*(k18*x11 - k19*x9*x1) + 1*k3*(k20*x15 - k21*x11*x3) + (-1)*k3*(k26*x11 - k27*x13*x3) + (-1)*k3*k33*x11 + 1*k3*k42*x20*x12 + 1*k3*k48*x17*x12)/k3, (1*k3*k33*x11 + (-1)*k3*k42*x20*x12 + (-1)*k3*k48*x17*x12)/k3, ((-1)*k3*(k16*x13 - k17*x8*x1) + 1*k3*(k26*x11 - k27*x13*x3) + (-1)*k3*k32*x13 + 1*k3*k41*x20*x14 + 1*k3*k47*x17*x14)/k3, (1*k3*k32*x13 + (-1)*k3*k41*x20*x14 + (-1)*k3*k47*x17*x14)/k3, ((-1)*k3*(k10*x15 - k11*x9*x4) + (-1)*k3*(k12*x15 - k13*x10*x1) + (-1)*k3*(k14*x15 - k15*x8*x6) + (-1)*k3*(k20*x15 - k21*x11*x3) + (-1)*k3*k28*x15 + 1*k3*k43*x20*x16 + 1*k3*k49*x17*x16)/k3, (1*k3*k28*x15 + (-1)*k3*k43*x20*x16 + (-1)*k3*k49*x17*x16)/k3, ((-1)*k3*k34*x17 + 1*k3*k35*x18 + 1*k3*k36*x18*x19 + (-1)*k3*k44*x17*x2 + (-1)*k3*k45*x17*x5 + (-1)*k3*k46*x17*x7 + (-1)*k3*k47*x17*x14 + (-1)*k3*k48*x17*x12 + (-1)*k3*k49*x17*x16)/k3, (1*k3*k34*x17 + (-1)*k3*k35*x18 + (-1)*k3*k36*x18*x19 + 1*k3*k44*x17*x2 + 1*k3*k45*x17*x5 + 1*k3*k46*x17*x7 + 1*k3*k47*x17*x14 + 1*k3*k48*x17*x12 + 1*k3*k49*x17*x16)/k3, 0/k3, (1*k3*k37*x21 + (-1)*k3*k38*x20*x2 + (-1)*k3*k39*x20*x5 + (-1)*k3*k40*x20*x7 + (-1)*k3*k41*x20*x14 + (-1)*k3*k42*x20*x12 + (-1)*k3*k43*x20*x16)/k3, ((-1)*k3*k37*x21 + 1*k3*k38*x20*x2 + 1*k3*k39*x20*x5 + 1*k3*k40*x20*x7 + 1*k3*k41*x20*x14 + 1*k3*k42*x20*x12 + 1*k3*k43*x20*x16)/k3, 0/k3 ]
constraints=[ 1*k3*x1 + 1*k3*x2 + 1*k3*x4 + 1*k3*x5 + 1*k3*x6 + 1*k3*x7 + 1*k3*x11 + 1*k3*x12 + 1*k3*x13 + 1*k3*x14 + 1*k3*x15 + 1*k3*x16 - k3*k50, 1*k3*x3 + 1*k3*x4 + 1*k3*x5 + 2*k3*x6 + 2*k3*x7 + 1*k3*x9 + 2*k3*x10 + 1*k3*x11 + 1*k3*x12 + 2*k3*x15 + 2*k3*x16 - k3*k51, 1*k3*x8 + 1*k3*x9 + 1*k3*x10 + 1*k3*x11 + 1*k3*x12 + 1*k3*x13 + 1*k3*x14 + 1*k3*x15 + 1*k3*x16 - k3*k52, 1*k3*x17 + 1*k3*x18 - k3*k53, 1*k3*x20 + 1*k3*x21 - k3*k54 ]
name="BIOMD0000000200"
pol=true
mass=true
rev=12
irr=22
def=18
rat=true
desc="Bray1995_chemotaxis_receptorlinkedcomplex"
stoichMatrix=[ [ 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0], [ 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0,-1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0], [ 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ -1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0,-1, 0, 0, 0, 0,-1, 1, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0], [ 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0], [ 0, 0, 0,-1,-1,-1, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
reconStoichMatrix=[ [ 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0], [ 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0], [ 1,-1, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ -1, 1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0], [ 0, 0, 0, 0, 0, 0,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
kineticMatrix=[ [ 0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0], [ 0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0], [ 0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0], [ 0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ] ]
paramValues=[ 4 1 141/100 73/20000 1000000 447/50000 1000000 297 1000000 16/25 1000000 14/125 1000000 229/10000 1000000 393/10 1000000 727 1000000 787/100000000 1000000 511/10000 1000000 51/500 1000000 169/2500 1000000 31/2 227/10000 227/10000 227/10000 227/10000 227/10000 31/25000 37/1000 500000 7/20 6000000 6000000 6000000 6000000 6000000 6000000 30000000 30000000 30000000 30000000 30000000 30000000 1/400000 1/200000 1/400000 1/100000 1/500000 ]
paramNames=[ "Hill" "Bias" "cell" "__lp_r2_k1" "__lp_r2_k2" "__lp_r3_k1" "__lp_r3_k2" "__lp_r4_k1" "__lp_r4_k2" "__lp_r5_k1" "__lp_r5_k2" "__lp_r6_k1" "__lp_r6_k2" "__lp_r7_k1" "__lp_r7_k2" "__lp_r8_k1" "__lp_r8_k2" "__lp_r9_k1" "__lp_r9_k2" "__lp_r10_k1" "__lp_r10_k2" "__lp_r11_k1" "__lp_r11_k2" "__lp_r12_k1" "__lp_r12_k2" "__lp_r13_k1" "__lp_r13_k2" "__lp_r14_k1" "__lp_r15_k1" "__lp_r16_k1" "__lp_r17_k1" "__lp_r18_k1" "__lp_r19_k1" "__lp_r20_k1" "__lp_r21_k1" "__lp_r22_k1" "__lp_r23_k1" "__lp_r24_k1" "__lp_r25_k1" "__lp_r26_k1" "__lp_r27_k1" "__lp_r28_k1" "__lp_r29_k1" "__lp_r30_k1" "__lp_r31_k1" "__lp_r32_k1" "__lp_r33_k1" "__lp_r34_k1" "__lp_r35_k1" "__cm_k50" "__cm_k51" "__cm_k52" "__cm_k53" "__cm_k54" ]
speciesNames=[ "AA" "AAp" "W" "WAA" "WAAp" "WWAA" "WWAAp" "TT" "TTW" "TTWW" "TTWAA" "TTWAAp" "TTAA" "TTAAp" "TTWWAA" "TTWWAAp" "Y" "Yp" "Z" "B" "Bp" "SetYp" ]
