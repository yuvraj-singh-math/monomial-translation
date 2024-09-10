
#!/usr/bin/env bash
slow=false
julia script.jl 35 $slow
if [ $? != 0 ]; then
if [ $slow == false ]; then
slow=true
julia script.jl 35 $slow
if [ $? != 0 ]; then
echo fails at 35 species
exit 1
fi
fi
fi

