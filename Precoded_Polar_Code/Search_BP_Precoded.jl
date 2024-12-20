#=
2024/12/13
とりあえずプログラムを完成させる
凍結ビットの集合をとりあえず作る
    SCの誤り率が小さいものの組み合わせを一通り試す

2024/12/04
なんかめんどい、とりあえず variable rowの個数を減らせばいいのでは？
variable rowの個数を w以下にする
SC的にとりあえずやるか
    でなかったらどうしようか笑

どのrowをvariable rowに選択するか？

=#

#=
modeの説明
sc scoreにscの誤り率を使う
=#

include("Precoding_Search_Utility.jl")
using DataStructures

#=
Vector型ではなく、Setを使う
Argument:
F:凍結ビットの集合の集合 Vector{Vector{Int}}、SC復号の誤り率を計算してそれで選択する
動的凍結ビットはこの集合中から選ぶ
つまり動的凍結ビットの集合をAとしたら、AはFの要素である
N:符号長
G:polar符号の生成行列、2種類ある？
=#
function Search_BP_PrecodingMatrix(F,d,zvalue,N,K,G,mode=:sc)
    #=
    Fの最大値を見つける
    maxF Vector{Int} それぞれの凍結ビットの集合の最大値が格納されている
    =#
    maxF = maximum.(F)
    #Setにすることで重複が消える
    maxF_set = Set(maxF)

    #参考論文のhに対応する
    a = CalcUpperI(d,G)

    #=
    初期のprecoding matrixを作成する
    initsetは初期precoding matrixの集合
    どうするか？凍結ビットも格納するようにする
    =#
    initset = InitialPrecodingMatrix(a,maxF,N)

    #initから優先度付きキューを作成する
    Q = PriorityQueue{Matrix{Int},Float64}()
    for z in initset
        #GetFrozenSetはSet型を返す
        fr = GetFrozenSet(z)
        score = GetScore(fr,D,F,mode=mode)
        insert!(Q,(z,fr),score)
    end

    #凍結ビットの集合から計算する方法？

    while !isempty(Q)
        p = dequeue!(Q)
        k,_ = size(l)
        if k == K
            return l
        else
            #=
            BP復号の場合このzをどう作成するか
            例えば、w=2ならzはどうなるか？
            w=2は何を意味するのか？
                動的凍結ビットに関係する情報ビットの個数の最大値
            とりあえずSCで作る
                SCで作る場合のやり方
                    VR : variable row 
                    未完なZにvR
            =#
            VR = CreateVR()
            for vr in VR
                z = [l;vr]
                b = GetFrozenBitForVR(vr)
                push!(fr,b)
                if MinDistance(z) >= d
                    score = GetScore(fr,D,F,mode=mode)
                    #キューに追加する
                    insert!(Q,(z,fr))
                end
            end
            
        end
    end
end

#=
FはSet
frは凍結ビットのSet
=#

function GetScore(fr,F,D;mode)
    #=
    return 1.0 のときどうするか？
    =#
    if mode == :sc
        e = 1.0
        for f in F
            if issubset(fr,f)
                v = D[f]
                e = min(e,v)
            end
        end
        return e
    end
end


#=
B = fix可能なinfoのインデックスはmaxFより大きくてN以下の自然数
C = 最小距離を考慮したときにfix可能なインデックスはaより大きくN以下の自然数
=#
function InitialPrecodingMatrix(a,maxF,N)
    init = Matrix{Int}[]
    for k in maxF
        if a < k
            #=
            [0.....0]
            =#
            z = zeros(Int,N-k,N)
            for j=k+1:N
                z[j-k,j] = 1
            end
            push!(init,z)
        end
    end
    init
end

#=
Aとはなにか？
=#
function CalcUpperI(d,G)
    N,_= size(G)
    for i=N:-1:1
        @views w = sum(G[i,:])
        if 2^w < d
            return i
        end
    end
    return -1
end

#=
variable rowから凍結ビットを求める
=#
function GetFrozenBitForVR(z)
    N = length(z)
    for i=N:-1:1
        if z[i] == 1
            return i
        end
    end
    return -1
end

#=
Bを作る
=#
function CreateB(F,fr,U)
    B = Int[]
    for i in F
        if !issubset(i,fr)
            continue
        end

        u = setdiff(U,fr)
        m = maximum(intersect(u,i))
        push!(B,m)
    end
    B
end

#=
BからVRを作る
個数が必要になるがいくつになるだろうか？
2^kループ
=#
function CreateVR(B,fr,N)
    VR = Vector{Int}[]
    Z = Set()
    K = setdiff(set(1:N),set(fr))
    for b in B
        U = Set(1:b-1)
        A = union(U,K)
        k = length(A)
        for i=1:2^k, j=1:N
            z = zeros(Int,N)
            if j == b
                z[j] = 1
            elseif j in A
                
            else

            end
        end
    end
end
