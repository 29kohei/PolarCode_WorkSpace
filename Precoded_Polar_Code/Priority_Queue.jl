#=
2024/12/11
Juliaで優先度付きキューを使う
パッケージをインストールする必要はあるか？
using DataStructures
でPriority Queueを使うことができる
アルゴリズムとかデータ構造ってほんとに面白い

ソートアルゴリズム好き

関数名は大文字のほうが好き

=#
using DataStructures

function MyPriority_Queue()
    #優先度付きキューの作成
    q = PriorityQueue{Matrix{Int},Int}()
    z = [1 1 1; 1 1 1]
    #キューに挿入する関数

    enqueue!(q,z,100)
    z = [1 0 1; 1 1 1]
    enqueue!(q,z,10)
    #最小値をキューから取り出す関数
    dequeue!(q)
end

#=
2024/12/11
ソートアルゴリズムで遊ぶ
アルゴリズムの評価　
    時間計算量  
    空間計算量　メモリをどれくらい使うか　少ないほうがいい
バブルソート

とにかく最小値を一番右（昇順の場合、最大値を一番右）にもっていきたい
=#
function Bubble_Sort(v=rand(5))
    #=
    バブルソートのアルゴリズム
    =#

    #vのベクトルの長さ
    N = length(v)

    #比較用にベクトルをコピー
    l = copy(v)

    #=
    計算量はO(N^2) 
    ループをみるとわかる

    jは今選択されている値、それが最小なら末尾までスワップされる
    バブルソートはスワップを繰り返すだけでソートが行える、楽
    選択ソートとかは最小値の選択とかいう部分が結構曖昧、最小値の選択放って色々ありそうだから
    バブルソートは操作が明確である

    =#

    for i=1:N-1, j=1:N-i
        if v[j] > v[j+1]
            v[j],v[j+1] = v[j+1],v[j]
        end
    end
    display([l v])

end

#=
2024/12/12
選択ソートは最小値を末尾に持ってくるだけでシンプルで直感的

降順にソートする
選択ソートは内容覚えてたわ、
=#
function SelectSort(v=rand(5))
    N = length(v)
    for l=N:-1:2
        #最小値をもとめる
        @views m = minimum(v[1:l])
        #最小値と等しい場所のインデックスを一つだけもとめる
        k = findfirst(isequal(m),v[1:l])
        v[k],v[l] = v[l],v[k]
    end
    v
end

#=
2024/12/13
プログラムがむずい
再帰関数でかけそうだが、なぜか？
上から実行されていく
=#
function MergeSort(v=rand(7))
    N = length(v)
    #=
    ベクトルの長さが1のときに返す
    =#
    if N==1
        return v
    end
    
    #=
    半分に分ける
    Nが奇数のとき、どうするか？
    =#
    d = ceil(Int,N/2)
    a = MergeSort(v[1:d])
    b = MergeSort(v[d+1:N])
    
    #=
    最後に結合する
    a,bはソートされている
    a,bを結合する
    aとbの先頭を比べる
    小さい方をxにpushする
    =#
    x = eltype(v)[]

    ai = bi = 1
    Na = length(a); Nb = length(b)
    while ai <= Na && bi <= Nb
        if a[ai] < b[bi]
            push!(x,a[ai])
            ai += 1
        else 
            push!(x,b[bi])
            bi += 1
        end
    end

    if ai <= Na
        return [x; a[ai:Na]]
    end
    
    if bi <= Nb
        return [x; b[bi:Nb]]
    end

    return x
end

#=
python とかだとdefだしなんだdefって、きもいな
=#
function QuickSort()
end


function InsertSort(v)
    
end










