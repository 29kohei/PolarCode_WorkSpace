include("../SCLsolver.jl")
include("../FrozenBitsSelect.jl")

#=
凍結ビット集合のソート問題
　ソートするべきかどうか
=#
"""
sim,シミュレーション回数

ebno,エネルギー対雑音密度

N,128

K,64

L,リスト
"""
function main(;sim=100,ebno=3.0,N=128,K=64,L=1,seed=8)
    F = frozenbits_select_GN(N,K); sort!(F)
    A = infobits_select(F,N); sort!(A)
    pm = PathManager(N,L)
    p = Polar(F,A)
    sclproblem = SCLProblem(N,K,pm,p)
    awgn = AWGN(N,K,seed,ebno)
    result = solver(sim,sclproblem,awgn)
    fn = "result_SCL.txt"
    output = """
    手法：LLRを用いない、確率を計算するSCL復号
    符号長:$N
    情報ビット数:$K
    リストサイズ:$L
    凍結ビット集合:$F
    情報ビット集合:$A
    エネルギー対雑音密度:$(ebno)
    シミュレーション回数:$(sim)
    """
    open(fn,"a") do file
        println(file,output)
        println(file,result)
    end
end