# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20","k21","k22","k23","k24","k25","k26","k27","k28","k29","k30","k31","k32","k33","k34","k35","k36","k37","k38","k39","k40","k41","k42","k43","k44","k45","k46","k47"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","x16","x17","x18","x19","x20","x21","x22","x23","x24","x25","x26","x27","x28"])
chemSystem=[((-1)*k1*k12*x1*k35 + 1*k1*k22*x3)/k1, (1*k1*k12*x1*k35 + (-1)*k1*k19*x2)/k1, (1*k1*k19*x2 + (-1)*k1*k22*x3)/k1, 0, 0, (1*k2*k14*x10 + (-1)*k2*k21*x6)/k2, (1*k2*k18*x12 + (-1)*k2*k29*x7*x26)/k2, (1*k2*k11*x9*k36 + (-1)*k2*k24*x8*x16)/k2, ((-1)*k2*k11*x9*k36 + 1*k2*k24*x8*x16)/k2, (1*k2*k6*x14*x13 + 1*k2*k13*x11*x13 + (-1)*k2*k14*x10)/k2, (1*k7*x14*x2 + (-1)*k2*k13*x11*x13 + (-1)*k2*k15*x11)/k2, ((-1)*k2*k18*x12 + 1*k2*k29*x7*x26)/k2, (1*k2*k4*x22*x8 + (-1)*k2*k5*x13)/k2, ((-1)*k2*k6*x14*x13 + (-1)*k7*x14*x2 + 1*k2*k15*x11 + 1*k2*k21*x6)/k2, (1*k2*k8*x20*x11 + (-1)*k2*k9*x15*x12 + 1*k2*k16*x21 + (-1)*k2*k17*x15 + (-1)*k2*k30*x15*x25)/k2, (1*k2*k3*x23 + (-1)*k2*k10*x16*x15 + (-1)*k2*k20*x16*x21)/k2, ((-1)*k2*k23*x17*x21 + (-1)*k2*k25*x17*x15 + 1*k2*k26*x19 + 1*k2*k27*x18 + (-1)*k2*k28*x17*x8)/k2, ((-1)*k2*k27*x18 + 1*k2*k28*x17*x8)/k2, (1*k2*k23*x17*x21 + 1*k2*k25*x17*x15 + (-1)*k2*k26*x19)/k2, ((-1)*k2*k8*x20*x11 + 1*k2*k17*x15)/k2, (1*k2*k9*x15*x12 + (-1)*k2*k16*x21 + 1*k2*k30*x15*x25)/k2, ((-1)*k2*k4*x22*x8 + 1*k2*k5*x13)/k2, ((-1)*k2*k3*x23 + 1*k2*k10*x16*x15 + 1*k2*k20*x16*x21)/k2, (1*k2*k31*x25 + (-1)*k32*x24*x2)/k2, ((-1)*k2*k31*x25 + 1*k32*x24*x2)/k2, (1*k33*x27*x2 + (-1)*k2*k34*x26)/k2, ((-1)*k33*x27*x2 + 1*k2*k34*x26)/k2, 0 ]
name="BIOMD0000000581"

