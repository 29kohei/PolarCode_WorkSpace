#=
2024/12/12
未完 MinDistance関数, MLerror関数
=#

#=
SC復号の誤り率を計算する関数
z値を格納したベクトル Vector{FLoat64}
Aは情報ビット、エラーを訂正するとき誤りが発生する場所
=#

#=
function SCerror(z,A)
    e = 1.0
    for i in A
        e *= 1.0-z[i]
    end
    1-e
end
=#

#=
シミュレーションでMLの確率を計算する
どうやって計算するか？
    GRANDを使うか？
    GRANDで計算しよう

生成行列Gから計算する
=#
function MLerror(G)
    
end

#=
普通に計算して求める
2^K個ある
@simd マクロって何？
Gを考える

Argument
G,生成行列 k×N
=#
function MinDistance(G)
    k,N = size(G)
    #=
    2^kの系列を生成する
    =#
    u = zeros(Int,k)
    #=
    u[i] == 1 && G[?,i] == 1
    が偶数個なら0
    =#
    for i=1:2^k
    end
end

