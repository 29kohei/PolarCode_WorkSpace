#=
issue : L=4ではパスが消えていた, 正しいパスが消える
issue : Lを大きくするとエラーが増える！おかしい!
=#
include("SCLUtility.jl")
@inline function update_contForks!(q,pm::PathManager)
    L = getListSize(pm)
    sort!(pm.index,by=x->pm.probForks[x])
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
            mm = getpathmetric(pm,l)
            probForks[1,l] = pathmetricaprox(mm,P[1],0)
            probForks[2,l] = pathmetricaprox(mm,P[1],1)
            q += 2
        else 
        #activeじゃないpath
            probForks[1,l] = floatmax()
            probForks[2,l] = floatmax()
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

    function recursivelyCalcLLR(lambda,phase)
        if isone(lambda)
            return nothing
        end

        if isodd(phase)
            recursivelyCalcLLR(lambda-1,phasetophai(phase))
        end
        
        for l = 1:L
            if !isactive(pm,l); continue; end
            P1 = getArrayPointer_P(pm,lambda,l)
            P2 = getArrayPointer_P(pm,lambda-1,l)
            C = getArrayPointer_C(pm,lambda,l)
            for b=1:2^(layer-lambda)
                P1[b] = isodd(phase) ? f(P2[2b-1],P2[2b]) : g(P2[2b-1],P2[2b],C[b][1])
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
    for i=1:N
        recursivelyCalcLLR(layer,i) 
        
        if !isfrozen(i)
            continuePaths_UnfrozenBit(pm,i)
        else 
            for l=1:L
                if !isactive(pm,l); continue; end
                C = getArrayPointer_C(pm,layer,l)
                P = getArrayPointer_P(pm,layer,l)
                C[1][phasetoindex(i)] = frozenbits_value 
                pm.pathmetric[l] = pathmetricaprox(pm.pathmetric[l],P[1],frozenbits_value)
                path_l = getpathpath(pm,l)
                path_l[i] = frozenbits_value
            end
        end
        if iseven(i); recursivelyUpdateC(layer,i); end
    end
    
    p = floatmax(); resultpath = 0
    for l=1:L
        if !isactive(pm,l); continue; end
        if pm.pathmetric[l] < p
            p = pm.pathmetric[l]
            resultpath = l
        end
    end

    return resultpath
end