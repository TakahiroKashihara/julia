using Distributed
rmprocs(workers())
addprocs(1)
@everywhere function plus(x)
    return(x + 1,x + 2)
end

function main()
    @distributed for i in 1:3
        open("./ditest.dat","a+") do f
            fi,gi = plus(i)
            
            println(f,"$(i) $(fi) $(gi)")            
        end
    end
end

@time main()


rmprocs(workers())
