
#!/usr/bin/env bash
slow=false
julia script.jl 8 $slow
if [ $? != 0 ]; then
if [ $slow == false ]; then
slow=true
julia script.jl 8 $slow
if [ $? != 0 ]; then
echo fails at 8 species
exit 1
fi
fi
fi

