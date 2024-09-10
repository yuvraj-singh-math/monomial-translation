
#!/usr/bin/env bash
slow=false
julia script.jl 18 $slow
if [ $? != 0 ]; then
if [ $slow == false ]; then
slow=true
julia script.jl 18 $slow
if [ $? != 0 ]; then
echo fails at 18 species
exit 1
fi
fi
fi

