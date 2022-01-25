using Distributed
addprocs(1)
@everywhere function plus(x)
    return(x + 1)
end

function main()
    @distributed for i in 1:3
        open("./ditest.dat","a+") do f
            fi = plus(i)
            println(f,"$(i) $(fi)")            
        end
    end
end

@time main()