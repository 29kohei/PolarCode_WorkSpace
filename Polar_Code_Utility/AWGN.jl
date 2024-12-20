using Random
#=
AWGN構造体

=#
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


#=
a::AWGN
a.LLRをSNRを使って
=#
function updatellr!(ch::AWGN)
    llr = ch.llr; rng = ch.rng; snrdb = ch.snrdb
    for i in eachindex(llr)
        y = -10^(snrdb/20) + sqrt(1/2)*randn(rng)
        llr[i] = -4*y*(10^(snrdb/20))
    end
    nothing
end