using Oscar;
using OscarODEbase;

if length(ARGS)==1
    global directory=ARGS[1]
end

function perturbInfo(sys::ODEbaseModel)
    location=sys.ID
    if @isdefined directory
        file="$directory/$location-matrix.csv"
    else
        file="../out/perturb_info/$location-matrix.csv"
    end
    
    if !isfile(file)
        error("System has not been analysed")
    end
    list=[]
    open(file, "r") do io
        lines=readlines(io)
        for line in lines
            data=split(line,"],")
            data=vcat([String(data[1])*"]"],split(data[2],","))
            data=vcat([String(data[1])],parse.(Ref(Int),data[2:end]))
            push!(list,data)
        end
    end
    return list
end

unfiltered_systems=get_odebase_model.(ODEbaseModels);

systems=filter(x->isfile("../out/perturb_info/$(x.ID)-matrix.csv"),unfiltered_systems)
systems=filter(x->x.massAction,systems);

# 1        2                                          3                                             4                              5                       6
# name, num nontrivial minors, num zero nontrivial minors, num zero minors, num minors, num columns
global score1(list::Vector)=-list[5]
global score2(list::Vector)=list[4]
global score3(list::Vector)=list[3]
global score4(list::Vector)=list[4]/list[5]
global score5(list::Vector)=list[3]/list[2]

global scores=[score1,score2,score3,score4,score5]

global result=[]

for sys in systems
    data=perturbInfo(sys)
    wins=[Int(score(data[end])>=max([score(datum) for datum in data]...)) for score in scores]
    push!(result,vcat([sys.ID],wins))
    println(wins)
end

csv=[]
for datum in result
    for j in 1:length(datum)
        push!(csv,"$(string(datum[j]))")
        if j+1<=length(datum)
            push!(csv,",")
        end
    end
    push!(csv,"
")
end
file=string(csv...);
open("score-results.csv", "w") do io
    write(io, file)
end
