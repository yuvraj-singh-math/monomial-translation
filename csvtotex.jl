list=[]
file="alignment-results.csv"
open(file, "r") do io
lines=readlines(io)
for line in lines
    data=split(line,",")
    #data=vcat(data[1],parse.(ref(int),data[2:end]))
    push!(list,data)
end
end

stats=[]
for data in list
	name=data[1][end-3:end]
	mean=data[4]
	worst=data[5]
	best=data[3]
	original=data[2]
	upperbound=data[6]
	row=[name,original,best,mean,worst,upperbound]
	println(row)
	push!(stats,row)
end

csv=[]
global i=0
for result in stats
    push!(csv,
          "
    \\coordinate (original$(result[1])) at ($(i),$(result[2]));
    \\coordinate (best$(result[1])) at ($(i),$(result[3]));
    \\coordinate (average$(result[1])) at ($(i),$(result[4]));
    \\coordinate (worst$(result[1])) at ($(i),$(result[5]));

    \\draw
    (best$(result[1]))++(-0.35,0) -- ++(0.7,0)  % horizontal line at best
    (worst$(result[1]))++(-0.35,0) -- ++(0.7,0) % horizontal line at worst
    (best$(result[1])) -- (worst$(result[1]));       % vertical line in between
    \\fill[black] (original$(result[1])) circle (2.5pt);
    \\fill[white,draw=black] (average$(result[1])) circle (1.5pt);

"
          )
    global i=i+1
end

file=string(csv...);
open("graph.tex", "w") do io
    write(io, file)
end
