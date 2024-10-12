# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"])
chemSystem=[(-1)*k1*(k2*x1*x10 - k3*x2)/k1, (1*k1*(k2*x1*x10 - k3*x2) + (-1)*k1*(k4*x2*x4 - k5*x3) + 1*k1*k10*x7)/k1, (1*k1*(k4*x2*x4 - k5*x3) + (-1)*k1*(k6*x3 - k7*x6*x5))/k1, ((-1)*k1*(k4*x2*x4 - k5*x3) + 1*k1*k11*x9)/k1, (1*k1*(k6*x3 - k7*x6*x5) + (-1)*k1*(k8*x5*x8 - k9*x7))/k1, 1*k1*(k6*x3 - k7*x6*x5)/k1, (1*k1*(k8*x5*x8 - k9*x7) + (-1)*k1*k10*x7)/k1, (-1)*k1*(k8*x5*x8 - k9*x7)/k1, (1*k1*k10*x7 + (-1)*k1*k11*x9)/k1, (-1)*k1*(k2*x1*x10 - k3*x2)/k1 ]
name="BIOMD0000000080"

