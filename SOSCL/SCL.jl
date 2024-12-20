#=
2024/11/22
N=64だとうまくいく
N=1024だとうまくいかない
    現象 pathmetricの値が0.0になる

=#

#=
issue : L=4ではパスが消えていた, 正しいパスが消える
issue : Lを大きくするとエラーが増える！おかしい!
=#
#=
soft output SCL
LLRではなくでW(|0) W(|1)を計算するようにする
変更する点
　recursivelyCalcLLR
    recursivelyCalcPにする、これは確率を計算する関数
    
=#
include("SCLUtility.jl")
@inline function update_contForks!(q,pm::PathManager)
    L = getListSize(pm)
    #LLR based SCLの場合 probForksの値が小さい順にソートする
    #sort!(pm.index,by=x->pm.probForks[x])

    #確率の値が大きい順にソートする
    sort!(pm.index,by=x->pm.probForks[x],rev=true)
    
    pm.contForks .= false
    for i=1:min(q,L)
        pm.contForks[pm.index[i]] = true
    end
    nothing
end

@inline function update_probForks!(pm::PathManager)
    q = 0
    probForks = pm.probForks
    layer = getlayer(pm)
    L = getListSize(pm)
    for l=1:L
        #activeなpathだけ調べる
        if isactive(pm,l)
            P = getArrayPointer_P(pm,layer,l)
            #pathmetric
            #mm = getpathmetric(pm,l)
            #probForks[1,l] = pathmetricaprox(mm,P[1],0)
            #probForks[2,l] = pathmetricaprox(mm,P[1],1)
            probForks[1,l] = P[1][begin]
            probForks[2,l] = P[1][end]
            q += 2
        else
        #activeじゃないpath
        #=
        #LLR用
            probForks[1,l] = floatmax()
            probForks[2,l] = floatmax()
        =#

        #確率用
            probForks[1,l] = -1.0
            probForks[2,l] = -1.0
        end
    end
    return q
end

@inline function continuePaths_UnfrozenBit(pm::PathManager,phase)
    q = update_probForks!(pm)
    update_contForks!(q,pm)

    L = getListSize(pm)
    contForks = pm.contForks
    probForks = pm.probForks
    layer = getlayer(pm)

    for l=1:L
        if !isactive(pm,l); continue; end
        if !contForks[1,l] && !contForks[2,l]; killPath(pm,l); end
    end

    for l=1:L
        if !contForks[2l-1] && !contForks[2l]
            continue
        elseif contForks[2l-1] && contForks[2l]
            j = clonePath(pm,l)
            C = getArrayPointer_C(pm,layer,l)
            C[1][phasetoindex(phase)] = 0 
            C = getArrayPointer_C(pm,layer,j)
            C[1][phasetoindex(phase)] = 1
            getpathpath(pm,l)[phase] = 0
            getpathpath(pm,j)[phase] = 1
            pm.pathmetric[l] = probForks[2l-1]
            pm.pathmetric[j] = probForks[2l]
        elseif contForks[2l-1]
            C = getArrayPointer_C(pm,layer,l)
            C[1][phasetoindex(phase)] = 0
            getpathpath(pm,l)[phase] = 0
            pm.pathmetric[l] = probForks[2l-1]
        elseif contForks[2l]
            C = getArrayPointer_C(pm,layer,l)
            C[1][phasetoindex(phase)] = 1
            getpathpath(pm,l)[phase] = 1
            pm.pathmetric[l] = probForks[2l]
        end
    end
    nothing
end

@inline function SCL(pm::PathManager,isfrozen)
    layer = getlayer(pm)
    L = getListSize(pm)
    N = pm._N

    function recursivelyCalcP(lambda,phase)
        if isone(lambda)
            return nothing
        end

        if isodd(phase)
            recursivelyCalcP(lambda-1,phasetophai(phase))
        end
        
        for l = 1:L
            if !isactive(pm,l); continue; end
            P1 = getArrayPointer_P(pm,lambda,l)
            P2 = getArrayPointer_P(pm,lambda-1,l)
            C = getArrayPointer_C(pm,lambda,l)
            for b=1:2^(layer-lambda)
                #=
                変更が必要
                Pの仕様を確認
                    LLRを格納する場合、Pはベクトルのベクトル
                        P::Vector{Vector{Float64}}
                        P[layer]:Vector{Float64} 長さ,2^(m-layer)
                    0と1それぞれの尤度を格納する場合はどうなる？
                        P::Vector{Vecotr{Vector{FLoat64}}}
                        P[layer]:Vector{Vector{Float64}} 長さ,2^(m-layer)
                        P[layer][i]:Vector{Float64} 長さ2 [W(|0),W(|1)]
                            P[layer][i][begin] = W(|0)
                            P[layer][i][end] = W(|1)
                            W(|0),W(|1)をもっと正確に書く 
                        
                        Pの再帰的な計算方法
                        構造的帰納法？
                            

                計算方法
                    関数fと関数gを使ってLLRを計算している
                    fとgの引数 u_{i-1}と２つのLLR

                    LLRではなく確率を計算するには？
                        確率用のfとｇに対応する関数

                =#
                #P1[b] = isodd(phase) ? f(P2[2b-1],P2[2b]) : g(P2[2b-1],P2[2b],C[b][1])
                P1[b][begin] = isodd(phase) ? f(P2[2b-1],P2[2b],0) : g(P2[2b-1],P2[2b],C[b][1],0)
                P1[b][end] = isodd(phase) ? f(P2[2b-1],P2[2b],1) : g(P2[2b-1],P2[2b],C[b][1],1)
            end
        end
        nothing 
    end

    function recursivelyUpdateC(lambda,phase) 
        phai = phasetophai(phase)
        for l=1:L
            if !isactive(pm,l); continue; end

            C1 = getArrayPointer_C(pm,lambda,l)
            C2 = getArrayPointer_C(pm,lambda-1,l) 

            for b=1:2^(layer-lambda)
                C2[2b-1][phasetoindex(phai)] = mod(C1[b][1]+C1[b][2],0:1)
                C2[2b][phasetoindex(phai)] = C1[b][2]
            end
        end
        if iseven(phai)
            recursivelyUpdateC(lambda-1,phai) 
        end 
    end

    frozenbits_value = 0

    #pathmetricをどうするか？
    for i=1:N
        recursivelyCalcP(layer,i)
        
        if !isfrozen(i)
            continuePaths_UnfrozenBit(pm,i)
        else 
        for l=1:L
                if !isactive(pm,l); continue; end
                C = getArrayPointer_C(pm,layer,l)
                P = getArrayPointer_P(pm,layer,l)
                C[1][phasetoindex(i)] = frozenbits_value 
                #pm.pathmetric[l] = pathmetricaprox(pm.pathmetric[l],P[1],frozenbits_value)
                #確率用のpathmetric 
                #frozenbitsの値が0なので、begin
                pm.pathmetric[l] = P[1][begin]
                println(pm.pathmetric[l] )
                path_l = getpathpath(pm,l)
                path_l[i] = frozenbits_value
            end
        end
        if iseven(i); recursivelyUpdateC(layer,i); end
    end
    
    #=
    pathを選択する
        pathmetricを使う
        pathmetricの仕様に依存している
    =#

    #=
    p = floatmax(); resultpath = 0
    for l=1:L
        if !isactive(pm,l); continue; end
        if pm.pathmetric[l] < p
            p = pm.pathmetric[l]
            resultpath = l
        end
    end
    =#

    p = -1.0; resultpath = 0
    for l=1:L
        if !isactive(pm,l); continue; end
        if p < pm.pathmetric[l]
            resultpath = l
            p = pm.pathmetric[l]
        end
    end

    return resultpath
end

