#module utility_funcs

using Printf

#export create_file, write_output, write_text

function create_file(path,evod,mmax,Nsites)
    f=open(path,"w")
    if evod == "dvr"
        println(f,"# Ngrid= ",mmax)
        println(f,"# DVR basis")
    else
        println(f,"# mmax= ",mmax)
        println(f,"# m-states: ",evod," states (FBR)")
    end
    println(f,"# Nr. of sites: ",Nsites)
    println(f)
    close(f)
end

function write_output(path,g,observable)
    text=' '
    for b=1:length(observable)
        text*=string(observable[b],' ')
    end
    f=open(path,"a")
    println(f,round(g,digits=4)," ",text)
    close(f)
end

function write_text(path,g,txt)
    f=open(path,"a")
    println(f,round(g,digits=4)," ",txt)
    close(f)
end

###
###save data from Observers, code by Wladislaw Krinitsin
# util function to save data in obs
function savedata(name::String, obs)
    name == "" && return
    h5open(name*".h5", "w") do file    
      # iterate through the fields of obs and append the data to the dataframe
      for n in names(obs)
        #@show n
        create_group(file, n)
        for (i,data) in enumerate(obs[:,n])
           #@show i,data
          file[n][string(i)] = [data[i] for i in 1:length(data)]
        end
      end 
    end
  end
  
  # util function to save data in obs
  function saveparams(name::String, params)
    name == "" && return
    h5open(name*".h5", "cw") do file    
      create_group(file, "params")
      for (name,val) in params
        file["params"][name] = val
      end 
    end
  end# util function to save data in obs

#end