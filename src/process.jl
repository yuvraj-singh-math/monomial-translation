using Oscar;
using OscarOdebase;
Random.seed!(13371337)
Oscar.set_seed!(13371337)

# limitations: returns the perturbation as a string
function perturbInfo(sys::OdebaseModel)
    location=sys.ID
    file="../out/perturb_info/$location-matrix.csv"
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

unfiltered_systems=get_odebase_system.(odebaseSystems);

systems=filter(x->isfile("../out/perturb_info/$(x.ID)-matrix.csv"),unfiltered_systems);
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
    up=max([score1(datum) for datum in data[1:end]]...)
    #original=data[end]
    filtered=filter(x->score1(x)==up,data)
    wins=[score(data[end])>max([score(datum) for datum in filtered]...) for score in scores]
    if !(data[end] in filtered)
        println("$(sys.ID) bad")
    end
    #temp=[length(filter(x->score(x)>=score(data[end]),filtered)) for score in scores]
    #println(temp)
    if length(filtered)>1
        println(sys.ID)
        for datum in data
            if score1(datum)==up
                println(datum)
            end
        end
        println(wins)
        push!(result,wins)
    end
end

