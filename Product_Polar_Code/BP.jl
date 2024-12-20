module BPs
include("AWGN.jl")
include("PolarEncode.jl")
"""
BPProblem
"""
struct BPProblem
    N::Int
    K::Int
    mi::Int
    M::Matrix{Int}
    L::Matrix{Float64}
    R::Matrix{Float64}
end

BPProblem(N,K,mi)  = begin 
    M = indexMatrix(N)
    n = Int(log2(N))
    L = zeros(N,n+1)
    R = zeros(N,n+1)
    BPProblem(N,K,mi,M,L,R)
end

"""

"""
function initialize!(t::BPProblem,ch::Channel,p::Polar)
    updatellr!(ch)
    initializeLandR!(t.L,t.R,p.F,ch.llr)
    nothing
end

"""
information bitsにエラーがあるかを調べる

```
u[i] != 0 (frozen value) 
```

ゼロでないビットが一つでもあるならfalseを返す
"""
function iserror(t::BPProblem)
    A = getinfobits(t.polar)
    @views !isnothing(findfirst(<(0.0),t.L[A,1]))
end

"""

"""
@inline function errors(t::BPProblem)
    A = getinfobits(t.polar)
    sum(ifelse(t.L[i,1] < 0.0,1,0) for i in A) 
end


@inline function f(x, y)
    return 0.9375*sign(x)*sign(y)*min(abs(x), abs(y))
end

function Gcheck(v,u)
    polar_encoder_FN!(v)
    u == v
end

harddecision(T,s...) = sum(s) < 0.0 ? one(T) : zero(T)

function BP(
    N,mi,
    L,R,M,
    checkcond)

    n = Int(log2(N))

    for iter=1:mi
        @views for m=1:n
            @simd for i in M[:,m]
                R[i,m+1] = f(R[i,m],L[i+N>>m,m+1]+R[i+N>>m,m])
                R[i+N>>m,m+1] = R[i+N>>m,m]+f(R[i,m],L[i,m+1])
            end
        end
        
        @views for m=n:-1:1
            @simd for i in M[:,m]
                L[i,m] = f(L[i,m+1],L[i+N>>m,m+1]+R[i+N>>m,m])
                L[i+N>>m,m] = L[i+N>>m,m+1]+f(L[i,m+1],R[i,m])
            end
        end

        if checkcond()
            return true
        end
    end
    
    return false
end

function initializeLandR!(L,R,F,llr)
    L .= 0.0
    R .= 0.0
    R[F,1] .= floatmax()
    L[:,end] .= llr
    nothing
end

"""
BP復号のsolver関数

sim:シミュレーション回数

t:BPProblem型の変数

ch:AWGN型
"""
function solver(sim,t::BPProblem,ch::AWGN,p::polar)
    n = Int(log2(t.N)); F = getfrozenbits(p); A = getinfobits(p)
    u = zeros(Int,t.N); v = zeros(Int,t.N)
    function checkcond()
        for i=1:t.N
            v[i] = harddecision(Int,t.L[i,end],t.R[i,end])
            u[i] = harddecision(Int,t.L[i,1],t.R[i,1])
        end 
        Gcheck(u,v)
    end

    bler = 0; ber = 0; errornum = zeros(Int,sim)
    
    for i=1:sim
        initialize!(t,ch,p)
        success = BP(t.N,t.mi,t.L,t.R,t.M,checkcond)
        if !iserror(t)
            continue
        end
        errornum[i] = errors(t)
        bler += 1
        ber += errornum[i]
    end
    return bler/sim,ber/(sim*t.K), errornum
end

"""
parameter

**input**
 
LLR:受信信号の対数尤度比
F:凍結ビットの集合 例:[1,3]

**output**

Gmatrixの終了判定を満たすかどうか
"""
function BP(llr,F)
    N = length(llr)
    K = N - length(F)
    bpproblem = BPProblem(N,K,100)
    u = zeros(Int,N)
    y = zeros(Int,N)
    function Gmatrix()
    end

    
end

end