#=
Output
多項式を返す

多項式の表現方法
例えば、
最大の次数(max(i|x^i))とx^iの係数 coefficient x[i+1] なぜなら、次数ゼロの係数がx[1]だから
struct Poly
    max_deg::Int
    coefficient::Vector{Int} #サイズ=max_deg+1
end
構造体を使っても問題ないか？
Poly構造体を作る
=#
#=
CalcA
cosetの重み分布を計算する
    cosetとはなにか？
    pathで決まる
=#
function CalcA(n,u)
    if n==1
        
    else
    end
end

#=
Argument
F:凍結ビットの集合
=#
function CEW(F,)
    #maxFを計算する
    maxF = maximum(F)
    #CalcAを計算する
    #=
    uはベクトル Vector{Int}
    length(u) = maxF
    =#
    for u in U   
        if u[maxF] == 0
        #u[maxF] == 1
        else
            
        end
    end
end

#=
凍結ビットの制約を満たす系列を生成する関数
=#
function f()
end

#=
ZはPrecoding Matrix
Zは行列
例えば、
Z = [
 0 0 1 0
 0 1 0 0
 1 0 0 1]

Zから凍結ビットの集合を求める
SC aimedの場合、行列から一意に決まる
=#

function PrecodingMatrixInfoSets(Z)
    #=
    x:行数
    y:列数
    =#
    x,_ = size(Z)
    for i=1:x
        @views a = findfirst(isequal(1),Z[i,:])
        #println(a)
    end
    a
end








