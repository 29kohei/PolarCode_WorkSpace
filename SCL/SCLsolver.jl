include("AWGN.jl")
include("Polar.jl")
include("SCL.jl")
#=
N,符号長
K,情報ビット数
pm,パスマネージャー構造体

SCLProblemの役割

依存性とは？
依存性
=#
struct SCLProblem
    N::Int
    K::Int
    pm::PathManager
    p::Polar
end

#=
2024/12/06
なんのためのコンストラクタなのこれ？
=#
SCLProblem(N,K,L) = begin
    p = Polar(BattGN,N,K)
    pm = PathManager(N,L)
    SCLProblem(N,K,pm,p)
end

#=
2024/12/06
どういうメソッドかわかりづらすぎる
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

