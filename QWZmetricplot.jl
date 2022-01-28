using Distributed
rmprocs(workers())
addprocs(1)
@everywhere using SharedArrays#@everywhereÇ™Ç¢ÇÁÇ»Ç¢ÅH
@everywhere using LinearAlgebra
@everywhere using ProgressBars
@everywhere using Plots

@everywhere struct Va
    BZ::Array{Any}
    
end
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

    return(Hermitian(H)/100)#1/100Ç…ÇµÇΩ
end
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

    return(Hermitian(H)/100)#1/100Ç…ÇµÇΩ
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

    return(Hermitian(H)/100)#1/100Ç…ÇµÇΩ
end
@everywhere function trgkm(k::Vector{Float64},M::Float64)
     
    v,w = eigen(H(k,M))
    pv = v[2]
    mv = v[1]
    pw = w[:,2]
    mw = w[:,1]
    
    #g_xx
    gxx = (abs(dot(mw,H_dx(k) * pw)))^2/(pv-mv)^2
    #g_yy
    gyy = (abs(dot(mw,H_dy(k) * pw)))^2/(pv-mv)^2
    #tr_g
    g = gxx + gyy
           
    return(g)
end
@everywhere function detgkm(k::Vector{Float64},M::Float64)
    v,w = eigen(H(k,M))
    pv = v[2]
    mv = v[1]
    pw = w[:,2]
    mw = w[:,1]
   
    #g_xx
    gxx = (abs(dot(mw,H_dx(k) * pw)))^2/(pv-mv)^2
    #g_yy
    gyy = (abs(dot(mw,H_dy(k) * pw)))^2/(pv-mv)^2
    #g_xy
    gxy = real((dot(mw,H_dx(k) * pw)*dot(pw,H_dy(k) * mw))/(pv-mv)^2)

    #det_g
    g = sqrt(abs(gxx*gyy - gxy^2))
    return(g)
end
@everywhere function trgkm(k::Vector{Float64},M::Float64)
     
    v,w = eigen(H(k,M))
    pv = v[2]
    mv = v[1]
    pw = w[:,2]
    mw = w[:,1]
    
    #g_xx
    gxx = (abs(dot(mw,H_dx(k) * pw)))^2/(pv-mv)^2
    #g_yy
    gyy = (abs(dot(mw,H_dy(k) * pw)))^2/(pv-mv)^2
    #tr_g
    g = gxx + gyy
           
    return(g)
end
@everywhere function detg(Nk::Int,M::Float64)
    
    kx = range(0.,2*pi,length = Nk)
    ky = range(0.,2*pi,length = Nk)
    BZ = Array{Any}(undef,length(kx),length(ky))
    for i in 1:length(ky)
        for j in 1:length(kx)
            BZ[j,i] = [kx[i],ky[j]]
        end
    end
    va = Va(BZ)
    A = SharedArray{Float64}(length(kx),length(ky))
    @sync @distributed for i = 1:Nk
        for j in 1:Nk
            A[i,j] = detgkm(va.BZ[i,j],M)
        end
    end
    return(sum(A)*(2*pi/Nk)^2)
end

function main()
    M = 1.5:0.05:3.0
    NK = [100,200,300,400]
    plot(M,detg.(NK[1],M),label = "Nk = 100")
    plot!(M,detg.(NK[2],M),label = "Nk = 200")
    plot!(M,detg.(NK[3],M),label = "Nk = 300")
    plot!(M,detg.(NK[4],M),label = "Nk = 400")
end


