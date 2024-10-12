# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15"])
chemSystem=[(1*k17*k6 + (-1)*k17*k7*x1 + (-1)*k17*k8*x5*x1 + (-1)*k17*k9*x5*x1)/k17, (1*k17*k8*x5*x1 + (-1)*k17*k2*x2 + (-1)*k17*k12*x2 + 1*k17*k11*x4 + (-1)*k17*k13*x6*x2)/k17, (1*k17*k12*x2 + (-1)*k17*k3*x3 + (-1)*k17*k13*x6*x3)/k17, (1*k17*k9*x5*x1 + (-1)*k17*k11*x4 + (-1)*k17*k7*x4)/k17, (1*k17*k4*x3 + (-1)*k17*k5*x5)/k17, (1*k17*k10*x15 + 1*k17*k14*(x2 + x3)*x6 + (-1)*k17*k16*x6)/k17, (-1)*k17*k15*x7/k17, (2*k17*k15*x7 + (-1)*k17*k15*x8)/k17, (2*k17*k15*x8 + (-1)*k17*k15*x9)/k17, (2*k17*k15*x9 + (-1)*k17*k15*x10)/k17, (2*k17*k15*x10 + (-1)*k17*k15*x11)/k17, (2*k17*k15*x11 + (-1)*k17*k15*x12)/k17, (2*k17*k15*x12 + (-1)*k17*k15*x13)/k17, (2*k17*k15*x13 + (-1)*k17*k15*x14)/k17, (2*k17*k15*x14 + (-1)*k17*k15*x15)/k17 ]
constraints=[ ]
name="BIOMD0000000687"
pol=true
mass=false
rev=0
irr=25
def=0
rat=true
desc="Wodarz2007 - Cytomegalovirus infection model with cytotoxic T lymphocyte response"
stoichMatrix=[ [ 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1, 0,-1,-1, 1,-1, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 1, 0, 0,-1,-1, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0,-1,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,-1,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,1,1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 2,-1, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 2,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 2,-1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 2,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 2,-1, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 2,-1, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 2,-1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-1 ] ]
reconStoichMatrix=[ [ 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1, 0,-1,-1, 1,-1, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 1, 0, 0,-1,-1, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0,-1,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,-1,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,1,1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 2,-1, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 2,-1, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 2,-1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 2,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 2,-1, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 2,-1, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 2,-1, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-1 ] ]
kineticMatrix=[ [ 0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0], [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 ] ]
paramValues=[ 175/11 1/10 1/5 1 1 10 1/10 1/10 1/2 1/5 1/10 1/100 1/1000000 31/2 1 1/10 1 ]
paramNames=[ "R0" "a0" "a1" "k" "u" "lambda" "d" "beta" "gamma" "alpha" "phi" "eta" "pa" "ca" "r" "ba" "COMpartment" ]
speciesNames=[ "x_0" "y_0" "y_1" "L_0" "v_0" "z_a" "m_0_0" "m_1_0" "m_2_0" "m_3_0" "m_4_0" "m_5_0" "m_6_0" "m_7_0" "m_8_0" ]
