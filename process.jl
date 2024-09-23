using Odebase;

# limitations: returns the perturbation as a string
function perturbInfo(sys::Odebase.OdebaseNode)
    location=sys.ID
    file="odebase/out/$location-matrix.csv"
    if !isfile(file)
        error("System has not been analysed")
    end
    list=[]
    open(file, "r") do io
        lines=readlines(io)
        for line in lines
            data=replace(line,r".*Elem"=>"")
            data=split(data,"],")
            data=vcat([String(data[1])*"]"],split(data[2],","))
            data=vcat([String(data[1])],parse.(Ref(Int),data[2:end]))
            push!(list,data)
        end
    end
    return list
end

unfiltered_systems=Odebase.get_odebase_system.(odebaseSystems);

systems=filter(x->(x.numSpecies==4)&&(x.massAction),unfiltered_systems);
