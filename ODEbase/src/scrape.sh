#!/usr/bin/env bash
# https://www.odebase.org/detail/1325
n=1325
#<1986
while [ $n -le 1987 ]; do
url="https://www.odebase.org/detail/$n"
html=$(curl -s $url)

id=$(echo "$html" | sed  '7q;d' | sed  's/^.*; //' | sed  's/:.*//')
dir="/home/main/Downloads/odebase/$id"
echo $dir
if [ -d $dir ]; then
	echo "yay!!!"
fi
desc=$(echo "$html" | sed  '7q;d' | sed  's/^.*: //')
irr=$(echo "$html" | sed  '91q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')
rev=$(echo "$html" | sed  '92q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')
def=$(echo "$html" | sed  '93q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')

if [ "$(echo "$html" |sed  '96q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')" = "Yes" ]; then
	rat=true
else
	rat=false
fi

if [ "$(echo "$html" |sed  '97q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')" = "Yes" ]; then
	pol=true
else
	pol=false
fi


if [ "$(echo "$html" |sed  '98q;d' | sed  's/^.*<\/strong> //' | sed  's/<\/p>*//')" = "Yes" ]; then
	mass=true
else
	mass=false
fi

n=$(( $n + 1))
done
