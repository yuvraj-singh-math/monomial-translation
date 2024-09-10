
#!/usr/bin/env bash
slow=false
julia script.jl 17 $slow
if [ $? != 0 ]; then
if [ $slow == false ]; then
slow=true
julia script.jl 17 $slow
if [ $? != 0 ]; then
echo fails at 17 species
exit 1
fi
fi
fi

