#=
how to use
AWGNオブジェクトの準備
 ch = AWGN(N,K,EbN0,seed)
BPProblemオブジェクトの準備
 bp = BPProblem(N,K,max_iter)

#sample code
N = 256; K= 128; EbN0 =3.0; seed = 8; max_iter = 100; sim = 10000
ch = AWGN(N,K,EbN0,seed)
bp = BPProblem(N,K,max_iter)
solver(sim,ch,bp)
=#

#FrozenBitsSelect.jl
function frozenbits_select_FN(N,K)    
    designSNR = 0; n = Int64(log2(N)); S = 10^(designSNR/10); 
    z = zeros(N)
    #z[1] = exp(-S)
    z[1] = 0.32 #design snr = 0.3
    #z[1]=0.36787944117144233
    for j=1:n
        u=2^j
        for t=1:u>>1
            T=z[t]
            z[t]=2*T-T^2
            z[u>>1+t]=T^2
        end
    end
    
    frozenbits = [1:N;]
    sort!(frozenbits,by=x->z[x],rev=true)
    for _=1:K; pop!(frozenbits); end
    return frozenbits
end

function frozenbits_select_GN(N,K)
    FN=frozenbits_select_FN(N,K)
    @. FN = (bit_reversal_n(FN-1,Int64(log2(N))))+1
end

function bit_reversal_n(b,n) 
    b = (b & 0xffffffff00000000) >> 32 | (b & 0x00000000ffffffff) << 32
    b = (b & 0xffff0000ffff0000) >> 16 | (b & 0x0000ffff0000ffff) << 16
    b = (b & 0xff00ff00ff00ff00) >> 8  | (b & 0x00ff00ff00ff00ff) << 8
    b = (b & 0xf0f0f0f0f0f0f0f0) >> 4  | (b & 0x0f0f0f0f0f0f0f0f) << 4
    b = (b & 0xcccccccccccccccc) >> 2  | (b & 0x3333333333333333) << 2
    b = (b & 0xaaaaaaaaaaaaaaaa) >> 1  | (b & 0x5555555555555555) << 1
    Int64(b>>(64-n))
end

#@inline infobits_select(frozenbits,N) = [i for i=1:N if !in(i,frozenbits)]
@inline infobits_select(frozenbits,N) = filter(p->!in(p,frozenbits),1:N)

function Zvalue(N)
    designSNR = 0; n = Int(log2(N)); S = 10^(designSNR/10); 
    z = zeros(N)
    #z[1] = exp(-S)
    z[1] = 0.32 #design snr = 0.3
    #z[1]=0.36787944117144233
    for j=1:n
        u=2^j
        for t=1:u>>1
            T=z[t]
            z[t]=2*T-T^2
            z[u>>1+t]=T^2
        end
    end
    z
end

struct RM; end
function frozenbits_select_FN(::Type{RM},N,K,r)
    z = Zvalue(N); n = Int(log2(N))
    G = kron_generate_G(n)
    index = filter(i -> r <= sum(G[i,:]),1:N)
    sort!(index,by=x->z[x],rev=true)
    [i for i=1:N if in(i,index)]
end

function frozenbits_select_GN(::Type{RM},N,K,r)
    FN = frozenbits_select_FN(RM,N,K,r)
    @. FN = (bit_reversal_n(FN-1,Int64(log2(N))))+1
end

function rKRM(N,r)
    m = Int(log2(N))
    sum(binomial(m,i) for i=0:r)
end

function rRateRM(N,K)
    m = Int(log2(N))
    for i=0:m
        if rKRM(N,i) > K
            return i
        end
    end
end

#Polar.jl
struct Polar
    F::Vector{Int}
    A::Vector{Int}
end
getfrozenbits(p::Polar) = p.F
getinfobits(p::Polar) = p.A

abstract type batt end

Polar(::Type{batt},N,K) = 
begin
    F = frozenbits_select_FN(N,K); sort!(F)
    A = infobits_select(F,N); sort!(A)
    Polar(F,A)
end

Polar(::Type{RM},N,K,r) = 
begin
    F = frozenbits_select_FN(RM,N,K,r); sort!(F)
    A = infobits_select(F,N); sort!(A)
    Polar(F,A)
end

#AWGN.jl
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

#IndexMatrix.jl
function indexMatrix(N)
    x = [1:N;]
    n = Int(log2(N))
    M = zeros(eltype(x),N-1, n)
    for k = 0:n-1
        for i = 0:1<<(k+1):N
            if (1<<k+i) < N
                for j=1:2^k; M[i+j,n-k] = x[i+j]; end;
            end
        end
    end
    A = reshape(M[M.>0],N>>1,n)
end

#PolarEncode.jl
polar_encoder_FN!(seq) = polar_encoder!(seq,true)
polar_encoder_GN!(seq) = polar_encoder!(seq,false)

ope_xor(x::Integer,y::Integer)= xor(x,y)
ope_xor(x::T,y::T) where {T<:AbstractVector}= ifelse(isequal(x,y),zero(T),one(T)) 

function polar_encoder!(seq,V)
    N = length(seq); n = Int(log2(N))
    for i = (V ? (1:n) : (n:-1:1))
        kernel = N >> i
        @simd for j=1:N
            if iszero(mod(div(j-1,kernel),2))
                seq[j] = ope_xor(seq[j],seq[j+kernel])
            end
        end
    end
    nothing
end

function polar_parity_check(seq,frozenbits,cache_seq=similar(seq))
    @. cache_seq = seq
    polar_encoder_GN!(cache_seq)
    error = 0
    @simd for i in frozenbits
        error += (cache_seq[i]!=0)
    end
    (error == 0)
end

view_polar_infoseq(seq,ib) = @view seq[ib]

function kron_generate_G(n)
    i = 1
    kernel = [1 0; 1 1]
    G = [1]
    while i <= n
        G = kron(kernel,G)
       i += 1
    end
    G
end

function polar_systematic_encoder(seq,N,A)
    n = Int(log2(N))
    G = kron_generate_G(n)
    uA = mod.(G[A,A]*seq,2)
    u = zeros(eltype(uA),N)
    u[A] .= uA
    polar_encoder_FN!(u)
    u
end
@inline function f(x, y)
    return 0.9375*sign(x)*sign(y)*min(abs(x), abs(y))
end

function Gcheck(v,u)
    polar_encoder_FN!(v)
    u == v
end

harddecision(::Type{T},s...) where {T<:Real} = sum(s) < 0.0 ? one(T) : zero(T)
harddecision(s...) = sum(s) < 0.0 ? one(Int) : zero(Int)

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


struct BPProblem
    N::Int
    K::Int
    mi::Int
    M::Matrix{Int}
    L::Matrix{Float64}
    R::Matrix{Float64}
    polar::Polar
end

BPProblem(N,K,mi)  = begin 
    M = indexMatrix(N)
    n = Int(log2(N))
    L = zeros(N,n+1)
    R = zeros(N,n+1)
    polar = Polar(batt,N,K)
    BPProblem(N,K,mi,M,L,R,polar)
end

function initialize!(L,R,F,ch::AWGN)
    updatellr!(ch)
    initializeLandR!(L,R,F,ch.llr)
    nothing
end

function iserror(t::BPProblem)
    A = getinfobits(t.polar)
    @views !isnothing(findfirst(<(0.0),t.L[A,1]))
end

function errors(t::BPProblem)
    A = getinfobits(t.polar)
    sum(ifelse(t.L[i,1] <= 0.0,1,0) for i in A) 
end

function solver(sim,t::BPProblem,ch::AWGN)
    n = Int(log2(t.N)); F = getfrozenbits(t.polar); A = getinfobits(t.polar)
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
        initialize!(t.L,t.R,F,ch)
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


