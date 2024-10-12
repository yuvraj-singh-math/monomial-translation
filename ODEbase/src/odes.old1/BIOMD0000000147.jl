# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,k73,k74,k75,k76,k77)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k42","k43","k44","k45","k46","k47","k48","k49","k50","k51","k52","k53","k54","k55","k56","k57","k58","k59","k60","k61","k62","k63","k64","k65","k66","k67","k68","k69","k70","k71","k72","k73","k74","k75","k76","k77"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22","x23","x24"])
chemSystem=[(1*k5*k6 + (-1)*k5*k7*x1 + 1*k5*k75*x10^2)/k5, (1*k5*k8*x1 + (-1)*k4*(k9*x2*x9 - k10*x4) + (-1)*k4*(k11*x2*x8 - k12*x5) + (-1)*k4*k19*x2 + (-1)*(k4*k25*x2 - k5*k26*x3))/k4, ((-1)*k5*(k13*x3*x10 - k14*x6) + (-1)*k5*k20*x3 + 1*(k4*k25*x2 - k5*k26*x3))/k5, (1*k4*(k9*x2*x9 - k10*x4) + (-1)*k4*(k17*x4*x8 - k18*x7) + (-1)*k4*k23*x4)/k4, (1*k4*(k11*x2*x8 - k12*x5) + (-1)*k4*(k15*x5*x9 - k16*x7) + (-1)*k4*k22*x5 + 1*k5*k27*x6)/k4, (1*k5*(k13*x3*x10 - k14*x6) + (-1)*k5*k21*x6 + (-1)*k5*k27*x6)/k5, (1*k4*(k15*x5*x9 - k16*x7) + 1*k4*(k17*x4*x8 - k18*x7) + (-1)*k4*k24*x7)/k4, ((-1)*k4*(k11*x2*x8 - k12*x5) + (-1)*k4*(k17*x4*x8 - k18*x7) + 1*k4*k22*x5 + 1*k4*k24*x7 + 1*k4*k31*x12 + 1*k4*k33*x14 + (-1)*k4*(k37*x11*x8 - k38*x12) + (-1)*k4*(k43*x16*x8 - k44*x14) + (-1)*(k4*k50*x8 - k5*k51*x10) + (-1)*k4*(k57*x19*x8 - k58*x21) + (-1)*k4*(k63*x24*x8 - k64*x23) + 1*k4*k68*x21 + 1*k4*k70*x23)/k4, ((-1)*k4*(k9*x2*x9 - k10*x4) + (-1)*k4*(k15*x5*x9 - k16*x7) + 1*k4*k23*x4 + 1*k4*k24*x7 + 1*k4*k31*x12 + 1*k4*k32*x11 + (-1)*k4*(k39*x14*x9 - k40*x12) + (-1)*k4*(k45*x16*x9 - k46*x11) + (-1)*k4*(k55*x19*x9 - k56*x24) + (-1)*k4*(k61*x21*x9 - k62*x23) + 1*k4*k69*x24 + 1*k4*k70*x23 + (-1)*k4*k74*x9)/k4, ((-1)*k5*(k13*x3*x10 - k14*x6) + 1*k5*k21*x6 + 1*k5*k34*x13 + (-1)*k5*(k41*x15*x10 - k42*x13) + 1*(k4*k50*x8 - k5*k51*x10) + (-1)*k5*(k59*x20*x10 - k60*x22) + 1*k5*k67*x22)/k5, ((-1)*k4*k32*x11 + (-1)*k4*(k37*x11*x8 - k38*x12) + 1*k4*(k45*x16*x9 - k46*x11))/k4, ((-1)*k4*k31*x12 + 1*k4*(k37*x11*x8 - k38*x12) + 1*k4*(k39*x14*x9 - k40*x12))/k4, ((-1)*k5*k28*x13 + (-1)*k5*k34*x13 + 1*k5*(k41*x15*x10 - k42*x13))/k5, (1*k5*k28*x13 + (-1)*k4*k33*x14 + (-1)*k4*(k39*x14*x9 - k40*x12) + 1*k4*(k43*x16*x8 - k44*x14))/k4, (1*(k4*k29*x16 - k5*k30*x15) + (-1)*k5*k35*x15 + (-1)*k5*(k41*x15*x10 - k42*x13))/k5, ((-1)*(k4*k29*x16 - k5*k30*x15) + (-1)*k4*k36*x16 + (-1)*k4*(k43*x16*x8 - k44*x14) + (-1)*k4*(k45*x16*x9 - k46*x11) + 1*k5*k47*x17)/k4, ((-1)*k5*k48*x17 + 1*k5*k49)/k5, (1*k5*k52 + (-1)*k5*k53*x18)/k5, (1*k5*k54*x18 + (-1)*k4*(k55*x19*x9 - k56*x24) + (-1)*k4*(k57*x19*x8 - k58*x21) + (-1)*k4*k65*x19 + (-1)*(k4*k71*x19 - k5*k72*x20))/k4, ((-1)*k5*(k59*x20*x10 - k60*x22) + (-1)*k5*k66*x20 + 1*(k4*k71*x19 - k5*k72*x20))/k5, (1*k4*(k57*x19*x8 - k58*x21) + (-1)*k4*(k61*x21*x9 - k62*x23) + (-1)*k4*k68*x21 + 1*k5*k73*x22)/k4, (1*k5*(k59*x20*x10 - k60*x22) + (-1)*k5*k67*x22 + (-1)*k5*k73*x22)/k5, (1*k4*(k61*x21*x9 - k62*x23) + 1*k4*(k63*x24*x8 - k64*x23) + (-1)*k4*k70*x23)/k4, (1*k4*(k55*x19*x9 - k56*x24) + (-1)*k4*(k63*x24*x8 - k64*x23) + (-1)*k4*k69*x24)/k4 ]
constraints=[ 1*k4*x5 + 1*k5*x6 + 1*k4*x7 + 1*k4*x8 + 1*k5*x10 + 1*k4*x12 + 1*k5*x13 + 1*k4*x14 + 1*k4*x21 + 1*k5*x22 + 1*k4*x23 - k4*k76 + k5*k77 ]
name="BIOMD0000000147"
pol=true
mass=false
rev=19
irr=32
def=18
rat=true
desc="ODea2007_IkappaB"
stoichMatrix=[ [ 1,-1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1], [ 0, 0,1,-1,-1, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 1, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 1, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0,-1, 0, 0,-1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,-1, 0, 0,-1, 0,0, 0,0,-1,0, 0,0, 0,-1, 0, 0,-1, 0, 0, 0, 1, 0, 1, 0, 0, 0,0], [ 0, 0,0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,-1, 0, 0,-1,0, 0,0, 0,0, 0,0,-1, 0, 0,-1, 0, 0, 0, 0, 0, 1, 1, 0, 0,-1,0], [ 0, 0,0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0,0, 0,0, 1,0, 0,0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 1,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 1, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0,-1, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0,-1,-1,1, 0,0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,-1,1, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,1,-1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,1,-1,-1, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 1, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 1, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0,0, 0,0, 1, 0, 0, 0,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0,0 ] ]
reconStoichMatrix=[ [ 1,-1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1], [ 0, 0,1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 1,-1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 1,-1, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0,0, 0,0,-1, 1,0, 0,0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,0], [ 0, 0,0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1,0, 0,0, 0, 0,0, 0,0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,-1,0], [ 0, 0,0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,0, 0,0, 1,-1,0, 0,0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 1,-1,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1,-1, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 1,-1, 1,1, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,-1,1, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,1,-1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,1,-1, 1,-1, 1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 1, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1,-1, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 1,-1, 0, 0,-1, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,0], [ 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0, 0, 0,0, 0,0, 1,-1, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,0 ] ]
kineticMatrix=[ [ 0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0], [ 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0 ] ]
paramValues=[ 3/50 3/100 1/100 1 1 231/1250000 21/1250 153/625 27/20 3/40 30 3/50000 30 3/50000 111/10 3/40 30 3/50000 3/25 3/25 3/50000 3/50000 9/125 9/25 9/500 3/250 207/250 207/500 9/500 3/250 3/25 3/125 3/50000 3/50000 9/50 9/50 30 3/50000 72/25 21/200 30 3/50000 30 3/50000 9/25 21/200 153/625 21/1250 267/6250000 27/5 3/625 381/12500000 21/1250 153/625 27/50 21/200 30 3/50000 30 3/50000 21/5 21/200 30 3/50000 9/50 9/50 3/50000 3/50000 9/250 9/50 9/500 3/250 207/500 0 99/50 1/10 0 ]
paramNames=[ "Total_IkBalpha" "Total_IkBbeta" "Total_IkBeps" "cytoplasm" "nucleus" "__lp_r2_tr2a" "__lp_r3_tr3a" "__lp_r4_tr1a" "__lp_r5_a1" "__lp_r5_d1_1" "__lp_r6_a4_1" "__lp_r6_d4_1" "__lp_r7_a4_2" "__lp_r7_d4_2" "__lp_r8_a7" "__lp_r8_d1_2" "__lp_r9_a4_3" "__lp_r9_d4_3" "__lp_r10_deg1_c" "__lp_r11_deg1_n" "__lp_r12_deg4_n" "__lp_r13_deg4_c" "__lp_r14_r1" "__lp_r15_r4" "__lp_r16_tp1a" "__lp_r16_tp2a" "__lp_r17_k2_a" "__lp_r18_k2_b" "__lp_r19_tp1b" "__lp_r19_tp2b" "__lp_r20_r5" "__lp_r21_r2" "__lp_r22_deg5_c" "__lp_r23_deg5_n" "__lp_r24_deg2_n" "__lp_r25_deg2_c" "__lp_r26_a5_3" "__lp_r26_d5_3" "__lp_r27_a8" "__lp_r27_d2_2" "__lp_r28_a5_2" "__lp_r28_d5_2" "__lp_r29_a5_1" "__lp_r29_d5_1" "__lp_r30_a2" "__lp_r30_d2_1" "__lp_r31_tr1b" "__lp_r32_tr3b" "__lp_r33_tr2b" "__lp_r34_k1_2" "__lp_r34_k1_1" "__lp_r35_tr2e" "__lp_r36_tr3e" "__lp_r37_tr1e" "__lp_r38_a3" "__lp_r38_d3_1" "__lp_r39_a6_1" "__lp_r39_d6_1" "__lp_r40_a6_2" "__lp_r40_d6_2" "__lp_r41_a9" "__lp_r41_d3_2" "__lp_r42_a6_3" "__lp_r42_d6_3" "__lp_r43_deg3_c" "__lp_r44_deg3_n" "__lp_r45_deg6_n" "__lp_r46_deg6_c" "__lp_r47_r3" "__lp_r48_r6" "__lp_r49_tp1e" "__lp_r49_tp2e" "__lp_r50_k2_e" "__lp_r51_k_IKK_deg" "__lp_r52_tr2a_i" "__cm_k76" "__cm_k77" ]
speciesNames=[ "IkBa_mRNA" "IkBa_cytoplasm" "IkBa_nucleus" "IkBaIKK" "IkBaNFkB_cytoplasm" "IkBaNFkB_nucleus" "IkBaIKKNFkB" "NFkB_cytoplasm" "IKK" "NFkB_nucleus" "IkBbIKK" "IkBbIKKNFkB" "IkBbNFkB_nucleus" "IkBbNFkB_cytoplasm" "IkBb_nucleus" "IkBb_cytoplasm" "IkBb_mRNA" "IkBe_mRNA" "IkBe_cytoplasm" "IkBe_nucleus" "IkBeNFkB_cytoplasm" "IkBeNFkB_nucleus" "IkBeIKKNFkB" "IkBeIKK" ]
