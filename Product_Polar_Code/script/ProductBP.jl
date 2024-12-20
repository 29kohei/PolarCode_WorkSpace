include("../BP.jl")
using Parameters
using Random

struct AWGN
    N :: Int
    K :: Int
    eb::Float64
    snrdb :: Float64
    R :: Float64
    rng :: MersenneTwister
    llr ::Vector{Float64}

    AWGN(N,K,seed,eb) = begin 
        R = K/N
        snrdb = eb + 10.0*log10(R)
        rng = MersenneTwister(seed)
        llr = zeros(N)
        new(N,K,eb,snrdb,R,rng,llr)
    end
end

function updatellr!(ch::AWGN)
    llr = ch.llr; rng = ch.rng; snrdb = ch.snrdb
    for i in eachindex(llr)
        y = -10^(snrdb/20) + sqrt(1/2)*randn(rng)
        llr[i] = -4*y*(10^(snrdb/20))
    end
    nothing
end

"""
依存関係:

Imax,Lmatrix更新回数の最大値

alpha,soft値の重み

Le La Lapp Lch,復号に必要な行列
"""
struct ProdTurboProblem
    Imax::Int
    alpha :: Float64
    Le :: Matrix{Float64}
    La :: Matrix{Float64}
    Lapp :: Matrix{Float64}
    Lch :: Matrix{Float64}
end

"""
コンストラクタ 
"""
ProdTurboProblem(N,Imax,alpha) = begin
    Le = zeros(N,N)
    La = zeros(N,N)
    Lapp = zeros(N,N)
    Lch = zeros(N,N)
    ProdTurboProblem(Imax,alpha,Le,La,Lapp,Lch)
end

"""

"""
function initialize!(t::ProdBPProblem,llr)
    t.La .= zero(eltype(t.La))
    t.Le .= zero(eltype(t.Le))
    t.Lapp .= zero(eltype(t.Lapp))

    nothing
end

function initialize!(i,::Type{ROW},t::ProdBPProblem)
    x = t.bpproblem; F = x.polar.F
    @views initializeLandR!(x.L,x.R,F,t.Lch[i,:])
    @views x.L[:,end] .=  x.L[:,end] .+ t.La[i,:]
    nothing
end

function initialize!(i,::Type{COLUMN},t::ProdBPProblem)
    x = t.bpproblem; F = x.polar.F
    @views initializeLandR!(x.L,x.R,F,t.Lch[:,i])
    @views x.L[:,end] .=  x.L[:,end] .+ t.La[:,i]
    nothing
end

function update!(i,::Type{ROW},t::ProdBPProblem)
    x = t.bpproblem
    for k=1:x.N
        t.Lapp[i,k] = x.L[k,end] + x.R[k,end]
    end
    nothing
end

function update!(i,::Type{COLUMN},t::ProdBPProblem)
    x = t.bpproblem
    for k=1:x.N
        t.Lapp[k,i] = x.L[k,end] + x.R[k,end]
    end
    nothing
end

function update!(t::ProdBPProblem)
    t.Le .= t.Lapp .- t.Lch .- t.La
    t.La .= t.alpha .*t.Le
    nothing
end

function errors(::Type{INFOBITS},t::ProdBPProblem)
    A = t.bpproblem.polar.A
    sum(ifelse(t.Lapp[i,j] < 0.0,1,0) for i in A, j in A) 
end

function errors(::Type{ALLBITS},t::ProdBPProblem)
    sum(ifelse(t.Lapp[i] < 0.0,1,0) for i in eachindex(t.Lapp)) 
end

"""
依存関係

BPProblem
"""
function Gmatrix(u,v,t::bpproblem)
    for i=1:bpproblem.N
        v[i] = harddecision(Int,bpproblem.L[i,end],bpproblem.R[i,end])
        u[i] = harddecision(Int,bpproblemR[i,1],bpproblem.R[i,1])
    end
    Gcheck(u,v)
end

"""

"""
function solver(sim,prodbpproblem::ProdBPProblem, ch::AWGN)
    @unpack Imax,alpha,Le,La,Lapp,Lch = prodbpproblem
    d = Dict()
    t = bpproblem;
    n = Int(log2(bpproblem.N))

    v = zeros(Int,bpproblem.N); u = zeros(Int,bpproblem.N)

    for s = 1:sim
        #初期化 + LLRの更新 
        for _=1:Imax
            #復号処理
            #row方向から復号する
            #符号語かどうかチェックする
            #次にcolumn方向を復号する
        end
    end
    
    return bler/sim, ber/sim/ch.K, sucnt, error
end