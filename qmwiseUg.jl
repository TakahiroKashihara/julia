using Distributed
rmprocs(workers())
addprocs(20)
@everywhere using SharedArrays#@everywhereがいらない？
@everywhere using LinearAlgebra
@everywhere using ProgressBars
@everywhere using Plots
@everywhere using LaTeXStrings

@everywhere struct Parm
    R::Float64
    Nk::Int
    Nw::Int
    d::Float64
    beta::Int
    mu::Float64
    M::Float64
    U::Float64
end
#自己エネルギー
@everywhere struct SE
    omegah::Vector{Float64}
    omega::Vector{Float64}
    imsA::Vector{Float64}
    resA::Vector{Float64}
    imsB::Vector{Float64}
    resB::Vector{Float64}
end
@everywhere struct Va
    BZ::Array{Any}
    BZ1::Vector{Any}
    hsbz::Vector{Any}
    
end

# これまでの定義
@everywhere function H(k::Vector{Float64},M::Float64)
    kx = k[1]
    ky = k[2]
    h11 = M + cos(kx) + cos(ky)
    h12 = complex(sin(kx),-sin(ky))
    h21 = conj(h12)
    h22 = -(M + cos(kx) + cos(ky))
    H = zeros(ComplexF64,2,2)
    H[1,1] = h11
    H[1,2] = h12
    H[2,1] = h21
    H[2,2] = h22

    return(Hermitian(H)/100)#1/100にした
end
# 新しい定義
@everywhere function H(k::Vector{Float64},p::Parm)
    kx = k[1]
    ky = k[2]
    h11 = p.M + cos(kx) + cos(ky)
    h12 = complex(sin(kx),-sin(ky))
    h21 = conj(h12)
    h22 = -(p.M + cos(kx) + cos(ky))
    H = zeros(ComplexF64,2,2)
    H[1,1] = h11
    H[1,2] = h12
    H[2,1] = h21
    H[2,2] = h22

    return(Hermitian(H)/100)#1/100にした
end
# 速度演算子は同じ
@everywhere function H_dx(k::Vector{Float64})
    kx = k[1]
    ky = k[2]
    h11 = -sin(kx)
    h12 = cos(kx)
    h21 = conj(h12)
    h22 = sin(kx)
    H = zeros(ComplexF64,2,2)
    H[1,1] = h11
    H[1,2] = h12
    H[2,1] = h21
    H[2,2] = h22

    return(Hermitian(H)/100)#1/100にした
end

@everywhere function H_dy(k::Vector{Float64})
    kx = k[1]
    ky = k[2]
    h11 = -sin(ky)
    h12 = complex(0,-cos(ky))
    h21 = conj(h12)
    h22 = sin(ky)
    H = zeros(ComplexF64,2,2)
    H[1,1] = h11
    H[1,2] = h12
    H[2,1] = h21
    H[2,2] = h22

    return(Hermitian(H)/100)#1/100にした
end


@everywhere function qmgk(k::Vector{Float64},p::Parm,se::SE)
    gxxk = 0.
    gyyk = 0.
    gxyk = 0.
    for NE in 1:length(se.omega)
        for NE1 in 1:length(se.omegah)
            if 1 <= (NE-NE1) && (NE + NE1) <= length(se.omega)#w-w1もw+w1も有限
                Sw = [complex(se.resA[NE],se.imsA[NE]) 0.;0. complex(se.resB[NE],se.imsB[NE])]
                GRw = inv(((se.omega[NE] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2)-(H(k,p)+Sw)))
                GAw = GRw'
                #GR(w-w1)\neq 0,GR(w+w1)\neq 0

                Swmw1 = [complex(se.resA[NE - NE1],se.imsA[NE - NE1]) 0.;0. complex(se.resB[NE - NE1],se.imsB[NE - NE1])]
                GRwmw1 = inv(((se.omega[NE - NE1] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2) - (H(k,p) + Swmw1)))
                GAwmw1 = GRwmw1' 
                Swpw1 = [complex(se.resA[NE + NE1],se.imsA[NE + NE1]) 0.;0. complex(se.resB[NE + NE1],se.imsB[NE + NE1])]
                GRwpw1 = inv(((se.omega[NE + NE1] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2)-(H(k,p)+Swpw1)))
                GAwpw1 = GRwpw1'
                
                #sigam_xx
                Fxx = tr(H_dx(k)*GRwpw1*H_dx(k)*(GRw - GAw))
                Sxx = tr(H_dx(k)*(GRw - GAw)*H_dx(k)*GAwmw1)
                
                gxxw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fxx + Sxx)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gxxk += gxxw/pi

                #sigma_yy
                Fyy = tr(H_dy(k)*GRwpw1*H_dy(k)*(GRw - GAw))
                Syy = tr(H_dy(k)*(GRw - GAw)*H_dy(k)*GAwmw1)

                gyyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fyy + Syy)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gyyk += gyyw/pi

                #sigma_xy

                Fxy = (tr(H_dx(k)*GRwpw1*H_dy(k)*(GRw - GAw)) + tr(H_dy(k)*GRwpw1*H_dx(k)*(GRw - GAw)))/2
                Sxy = (tr(H_dx(k)*(GRw - GAw)*H_dy(k)*GAwmw1) + tr(H_dy(k)*(GRw - GAw)*H_dx(k)*GAwmw1))/2

                gxyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fxy + Sxy)/se.omegah[NE1]/2/pi/se.omegah[NE1] 
                gxyk += gxyw/pi
            elseif 1 > (NE-NE1) &&  (NE + NE1) <= length(se.omega)
                Sw = [complex(se.resA[NE],se.imsA[NE]) 0.;0. complex(se.resB[NE],se.imsB[NE])]
                GRw = inv(((se.omega[NE] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2)-(H(k,p)+Sw)))
                GAw = GRw'
                #GR(w-w1) = 0,GR(w+w1)\neq 0
                Swpw1 = [complex(se.resA[NE + NE1],se.imsA[NE + NE1]) 0.;0. complex(se.resB[NE + NE1],se.imsB[NE + NE1])]
                GRwpw1 = inv(((se.omega[NE + NE1] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2)-(H(k,p)+Swpw1)))
                GAwpw1 = GRwpw1'

                #sigam_xx
                Fxx = tr(H_dx(k)*GRwpw1*H_dx(k)*(GRw - GAw))
                #Sxx = tr(H_dx(k)*(GRw - GAw)*H_dx(k)*GAwmw1)
                #sigma_yy
                gxxw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fxx)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gxxk += gxxw/pi
                
                #sigma_yy
                Fyy = tr(H_dy(k)*GRwpw1*H_dy(k)*(GRw - GAw))
                #Syy = tr(H_dy(k)*(GRw - GAw)*H_dy(k)*GAwmw1)

                gyyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fyy)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gyyk += gyyw/pi

                #sigma_xy
                Fxy = (tr(H_dx(k)*GRwpw1*H_dy(k)*(GRw - GAw)) + tr(H_dy(k)*GRwpw1*H_dx(k)*(GRw - GAw)))/2
 

                gxyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Fxy)/se.omegah[NE1]/2/pi/se.omegah[NE1] 
                gxyk += gxyw/pi
            elseif 1 <= (NE-NE1) && (NE + NE1) > length(se.omega)
                Sw = [complex(se.resA[NE],se.imsA[NE]) 0.;0. complex(se.resB[NE],se.imsB[NE])]
                GRw = inv(((se.omega[NE] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2)-(H(k,p)+Sw)))
                GAw = GRw'
                #GR(w-w1)\neq 0,GR(w+w1) = 0
                Swmw1 = [complex(se.resA[NE - NE1],se.imsA[NE - NE1]) 0.;0. complex(se.resB[NE - NE1],se.imsB[NE - NE1])]
                GRwmw1 = inv(((se.omega[NE - NE1] + p.mu + complex(0,p.d))*Matrix{Float64}(I,2,2) - (H(k,p) + Swmw1)))
                GAwmw1 = GRwmw1' 
                #sigam_xx
                #Fxx = tr(H_dx(k)*GRwpw1*H_dx(k)*(GRw - GAw))
                Sxx = tr(H_dx(k)*(GRw - GAw)*H_dx(k)*GAwmw1)
               
                gxxw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Sxx)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gxxk += gxxw/pi
                
                #sigma_yy
                #Fyy = tr(H_dy(k)*GRwpw1*H_dy(k)*(GRw - GAw))
                Syy = tr(H_dy(k)*(GRw - GAw)*H_dy(k)*GAwmw1)
                gyyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Syy)/se.omegah[NE1]/2/pi/se.omegah[NE1]
                gyyk += gyyw/pi 

                #sigma_xy
                Sxy = (tr(H_dx(k)*(GRw - GAw)*H_dy(k)*GAwmw1) + tr(H_dy(k)*(GRw - GAw)*H_dx(k)*GAwmw1))/2

                gxyw = - (1/(exp(p.beta*se.omega[NE])+1))*real(Sxy)/se.omegah[NE1]/2/pi/se.omegah[NE1] 
                gxyk += gxyw/pi
            
            end
        end
    end
    return((gxxk + gyyk)*(p.R/p.Nw)^2,(gxxk*gyyk-gxyk^2)*(p.R/p.Nw)^4,sqrt(abs((gxxk*gyyk-gxyk^2)*(p.R/p.Nw)^4)))
end

function ok()
    println("ok")
    return()
end
function mainU()
    p = Parm(#=reshapeと合わせる=#0.2, #=Nk=#10, #=Nw.これもreshapeと合わせる=#1000, 1/20/100, 10000, parse.(Float64,ARGS[2])/2, parse.(Float64,ARGS[1]), parse.(Float64,ARGS[2]))
    se = SE(swise(p.M,p.U)...)
    Nkh = Int(Nk/2)
    println("NK = $(p.Nk)")
   
    hsbz= []
       

    for i in 1:Nkh
        ky = 0.
        kx = (pi/Nkh)*i
        push!(hsbz,[kx,ky])
    end 
    for i in 1:Nkh
        kx = Float64(pi) 
        ky = (pi/Nkh) * i
        push!(hsbz,[kx,ky])
    end
    for i in 1:Nkh
        kx = Float64(pi) - (Float64(pi)/Nkh)*i
        ky = kx
        push!(hsbz,[kx,ky])
    end
    
    kx = range(-pi,pi,length = p.Nk)
    ky = range(-pi,pi,length = p.Nk )
 
    BZ = Array{Any}(undef,length(kx),length(ky))
    for i in 1:length(ky)
        for j in 1:length(kx)
            BZ[j,i] = [kx[i],ky[j]]
        end
    end
    BZ1 = []
    for i in 1:length(kx)
        for j in 1:length(ky)
            push!(BZ1,[kx[i],ky[j]])
        end
    end
    
    ok()
    
 
    
    va = Va(BZ,BZ1,hsbz)
    ok()
    #高対称点
    kxhs = SharedArray{Float64}(length(va.hsbz))
    kyhs = SharedArray{Float64}(length(va.hsbz))
    trghs = SharedArray{Float64}(length(va.hsbz))
    detghs = SharedArray{Float64}(length(va.hsbz))
    sqdetghs = SharedArray{Float64}(length(va.hsbz))
    @sync @distributed for i = 1:length(va.hsbz)
        kxhs[i] = va.hsbz[i][1]
        kyhs[i] = va.hsbz[i][2]
        trgk,detgk,sqdetgk = qmgk(va.hsbz[i],p,se)
        trghs[i] = trgk
        detghs[i] = detgk
        sqdetghs[i] = sqdetgk

    end
    
    plot(trghs,ylabel = "QuantumMetric",label = L"\mathrm{Tr}(g)",xticks = 
    ([1,Int(length(va.hsbz)/3),Int(length(va.hsbz)*2/3),length(va.hsbz)],[L"\Gamma",L"X",L"M",L"\Gamma"]))
    plot!(detghs,ylabel = "QuantumMetric",label = L"\mathrm{det}(g)",xticks = 
    ([1,Int(length(va.hsbz)/3),Int(length(va.hsbz)*2/3),length(va.hsbz)],[L"\Gamma",L"X",L"M",L"\Gamma"]))
    plot!(sqdetghs,ylabel = "QuantumMetric",label = L"\sqrt{\mathrm{det}(g)}",xticks = 
    ([1,Int(length(va.hsbz)/3),Int(length(va.hsbz)*2/3),length(va.hsbz)],[L"\Gamma",L"X",L"M",L"\Gamma"]))

    for i in 1:length(va.hsbz)
        open("./metric_hs.dat","a+") do f
            x = kxhs[i]
            y = kyhs[i]
            trgk = trghs[i]
            detgk = detghs[i]
            sqdetgk = sqdetghs[i]
            println(f,"$(x) $(y) $(trgk) $(detgk) $(sqdetgk)")
        end
    end

    println("hsbz_ok")
    #BZ全て
    kxs = SharedArray{Float64}(length(kx),length(ky))
    kys = SharedArray{Float64}(length(kx),length(ky))    
    trg = SharedArray{Float64}(length(kx),length(ky))
    detg = SharedArray{Float64}(length(kx),length(ky))
    sqdetg = SharedArray{Float64}(length(kx),length(ky))
    @sync @distributed for i = 1:p.Nk
        for j in 1:p.Nk
            kxs[i,j] = va.BZ[i,j][1]
            kys[i,j] = va.BZ[i,j][2]
            trgk,detgk,sqdetgk = qmgk(va.BZ[i,j],p,se)
            trg[i,j] = trgk
            detg[i,j] = detgk
            sqdetg[i,j] = sqdetgk
        end
    end
    
    for i in 1:p.Nk
        for j in 1:p.Nk
            open("./metric.dat","a+") do f
                x = kxs[i,j]
                y = kys[i,j]
                trgk = trg[i,j]
                detgk = detg[i,j]
                sqdetgk = sqdetg[i,j]
                #println("kx = $(kx), ky = $(ky),trg =  $(trgk),detg =  $(detgk),sqdetg =  $(sqdetgk)")

                println(f,"$(x) $(y) $(trgk) $(detgk) $(sqdetgk)")
            end
        end

    end
    println("write_ok")
    
    #符号
    trsgn = SharedArray{Float64}(length(kx),length(ky))
    detsgn = SharedArray{Float64}(length(kx),length(ky))
    @sync @distributed for i = 1:Nk
        for j in 1:Nk
            if trg[i,j] >=0
                sgn[i,j] = 1
            else
                sgn[i,j] = -1
            end
        end
    end
    @sync @distributed for i = 1:Nk
        for j in 1:Nk
            if detg[i,j] >=0
                sgn[i,j] = 1
            else
                sgn[i,j] = -1
            end
        end
    end
    ENV["GKSwstype"]="nul"
    heatmap(kx,ky,trsgn,xlabel = "\$kx\$",ylabel = "\$ky\$",title = L"sgn\mathrm{tr}(g)")
    savefig("./trsgn.png")
    heatmap(kx,ky,detsgn,xlabel = "\$kx\$",ylabel = "\$ky\$",title = L"sgn\mathrm{det}(g)")
    savefig("./detsgn.png")
    println("sgn_ok")

    #積分値
    open("./metric_int.dat","a+") do f
        trg_int = sum(trg)*(2*pi/p.Nk)^2
        sqdetg_int = sum(sqdetg)*(2*pi/p.Nk)^2
        println(f,"$(trg_int) $(sqdetg_int)")
    end
    println("int_ok")
    
    
end

include("retestserver.jl")

@time mainU()

println("completed,M=$(ARGS[1]),U=$(ARGS[2])")
rmprocs(workers())
#addprocs(2),Nk = 5で223秒

#=a = SharedArray{Float64}(2,2)
@distributed for i in 1:2
    for j in 1:2
        a[i,j] = i+j
    end
end=#
    