include("AWGN.jl")
include("Polar.jl")
include("SCL.jl")

struct SCLProblem
    N::Int
    K::Int
    pm::PathManager
    p::Polar
end

#=
SCLProblem(N,K,L) = begin
    p = Polar(BattGN,N,K)
    pm = PathManager(N,L)
    SCLProblem(N,K,pm,p)
end
=#

@inline function initialize!(t::SCLProblem,a::AWGN)
    updatellr!(a)
    initialize!(t.pm,a.llr)
    nothing
end

function solver(sim,t::SCLProblem,a::AWGN)
    v = [in(i,t.p.F) for i=1:t.N]
    @inline isfrozen(i)  = v[i]
    bler = 0; ber = 0
    for i=1:sim
        initialize!(t,a)
        p = SCL(t.pm,isfrozen)
        tmp = errors(t.pm,p)
        if !iszero(tmp)
            bler += 1
            ber += tmp
        end
    end
    return bler/sim,ber/sim/t.K
end


