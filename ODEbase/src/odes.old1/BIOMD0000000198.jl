# scraped from ODEBASE, all credit goes to them
using Oscar;
paramsRing,(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20)=rational_function_field(QQ,["k1","k2","k3","k4","k5","k6","k7","k8","k9","k10","k11","k12","k13","k14","k15","k16","k17","k18","k19","k20"])
polRing,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)=polynomial_ring(paramsRing,["x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12"])
chemSystem=[0, (-1)*k17*(k1*k18*x2 - k2*x3)/k17, (1*k17*(k1*k18*x2 - k2*x3) + (-1)*k3*k17*x3)/k17, (1*k3*k17*x3 + (-1)*k17*(k4*x4 - k5*x5))/k17, 1*k17*(k4*x4 - k5*x5)/k17, (-1)*k17*(k6*k18*x6 - k7*x7)/k17, (1*k17*(k6*k18*x6 - k7*x7) + (-1)*k8*k17*x7)/k17, (1*k8*k17*x7 + (-1)*k17*(k9*k18*x8 - k10*x9))/k17, (1*k17*(k9*k18*x8 - k10*x9) + (-1)*k17*(k11*x9 - k12*x10))/k17, 1*k17*(k11*x9 - k12*x10)/k17 ]
constraints=[ x11 - x5 + x10, x12 - x2 + x3 + x4 + x6 + x7 + x8 + x9, 1*k17*x2 + 1*k17*x3 + 1*k17*x4 + 1*k17*x5 - k17*k19, 1*k17*x6 + 1*k17*x7 + 1*k17*x8 + 1*k17*x9 + 1*k17*x10 - k17*k20 ]
name="BIOMD0000000198"
pol=true
mass=false
rev=5
irr=2
def=0
rat=true
desc="Stone1996 - activation of soluble guanylate cyclase by nitric oxide"
stoichMatrix=[ [ -1, 0, 0,-1, 0,-1, 0], [ -1, 0, 0, 0, 0, 0, 0], [ 1,-1, 0, 0, 0, 0, 0], [ 0, 1,-1, 0, 0, 0, 0], [ 0, 0, 1, 0, 0, 0, 0], [ 0, 0, 0,-1, 0, 0, 0], [ 0, 0, 0, 1,-1, 0, 0], [ 0, 0, 0, 0, 1,-1, 0], [ 0, 0, 0, 0, 0, 1,-1], [ 0, 0, 0, 0, 0, 0, 1], [ 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0 ] ]
reconStoichMatrix=[ [ -1, 1, 0, 0, 0,-1, 1, 0,-1, 1, 0, 0], [ -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 1,-1, 1, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 1,-1, 1, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ] ]
kineticMatrix=[ [ 1,0,0,0,0,1,0,0,1,0,0,0], [ 1,0,0,0,0,0,0,0,0,0,0,0], [ 0,1,1,0,0,0,0,0,0,0,0,0], [ 0,0,0,1,0,0,0,0,0,0,0,0], [ 0,0,0,0,1,0,0,0,0,0,0,0], [ 0,0,0,0,0,1,0,0,0,0,0,0], [ 0,0,0,0,0,0,1,1,0,0,0,0], [ 0,0,0,0,0,0,0,0,1,0,0,0], [ 0,0,0,0,0,0,0,0,0,1,1,0], [ 0,0,0,0,0,0,0,0,0,0,0,1], [ 0,0,0,0,0,0,0,0,0,0,0,0], [ 0,0,0,0,0,0,0,0,0,0,0,0 ] ]
paramValues=[ 700 800 850 20 1/5 700 800 850 5 25 8/5 1/50 11/100 1/40 8/125 11/250 1 1/2 14/125 36/125 ]
paramNames=[ "k1" "k2" "k3" "k4" "k5" "k6" "k7" "k8" "k9" "k10" "k11" "k12" "e5c" "e5c_NO" "e6c_NO" "ext" "cytosol" "NO" "__cm_k19" "__cm_k20" ]
speciesNames=[ "NO" "sGCfast" "NO_sGCfast" "NO_sGCfast_6coord" "NO_sGCfast_5coord" "sGCslow" "NO_sGCslow" "NO_sGCslow_6coord" "NO_sGCslow_6coord_NO_int" "NO_sGCslow_5coord" "NO_sGC_5coord_tot" "sGC_inact_tot" ]
