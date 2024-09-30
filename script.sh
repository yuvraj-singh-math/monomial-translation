#!/usr/bin/env bash
species=("2" "4" "5" "6" "7" "8" "9" "10" "11" "13" "14" "16" "17" "18" "20" "22" "27" "28" "30" "35" "52" "57" "86")

slow=false
for num in ${species[@]};do
julia script.jl $num $slow
if [ $? != 0 ]; then
if [ $slow == false ]; then
slow=true
julia script.jl $num $slow
if [ $? != 0 ]; then
echo fails at $num species
exit 1
fi
fi
fi
done