#=
依存関係
../SCLsolver.jl

=#

include("../SCLsolver.jl")
using Dates
function main()
    sim = 1000; N = 1024; K = 512; L = 1; eb = 3.0; seed = 8
    t = SCLProblem(N,K,L)
    a = AWGN(N,K,seed,eb)
    res = solver(sim,t,a)
    data = """
    $(today())
    SCL復号 
    符号長:$N
    情報ビット数:$K
    エネルギー対雑音密度:$eb
    seed:$seed
    リスト:$L
    シミュレーション結果
    BLER:$(res[1])
    BER:$(res[2])
    """
    open("result.txt","a") do f
        println(f,data)
    end
end

main()