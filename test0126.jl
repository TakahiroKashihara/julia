using Distributed
addprocs(1)
@everywhere using SharedArrays

@everywhere function plus(x)
    return(x+1,x+2)
end
function main()
    a = SharedArray{Float64}(10)
    b = SharedArray{Float64}(10)

    @sync @distributed for i = 1:10
    p,q = plus(i)
    a[i] = p
    b[i] = q
    end
    for i in 1:length(a)

        open("./ditest.dat","a+") do f
            fi = a[i]
            gi = b[i]
            println(f,"$(i) $(fi) $(gi)")            
        end
    end
end
@time main()
rmprocs(workers())