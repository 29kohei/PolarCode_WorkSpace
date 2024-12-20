include("PolarEncode.jl")

@inline function f(x, y)
    return 0.9375*sign(x)*sign(y)*min(abs(x), abs(y))
end

function Gcheck(v,u)
    polar_encoder_FN!(v)
    u == v
end

harddecision(::Type{T},s...) where {T<:Real} = sum(s) < 0.0 ? one(T) : zero(T)
harddecision(s...) = sum(s) < 0.0 ? one(Int) : zero(Int)

#=
依存関係
=#
"""
N 符号長 
L 
"""
function BP(
    N,mi,
    L,R,M,
    checkcond)

    n = Int(log2(N))

    for iter=1:mi
        @views for m=1:n
            @simd for i in M[:,m]
                R[i,m+1] = f(R[i,m],L[i+N>>m,m+1]+R[i+N>>m,m])
                R[i+N>>m,m+1] = R[i+N>>m,m]+f(R[i,m],L[i,m+1])
            end
        end
        
        @views for m=n:-1:1
            @simd for i in M[:,m]
                L[i,m] = f(L[i,m+1],L[i+N>>m,m+1]+R[i+N>>m,m])
                L[i+N>>m,m] = L[i+N>>m,m+1]+f(L[i,m+1],R[i,m])
            end
        end

        if checkcond()
            return true
        end
    end
    
    return false
end

#=
Lはなんでもいいのか？
Rはなんでもいいのか？
=#
"""
L BP復号でLのメッセージを格納する行列 \\
R BP復号でRのメッセージを格納する行列 \\
F 凍結ビットの集合 \\
llr LLRを格納しているベクトル 
"""
function initializeLandR!(L,R,F,llr)
    L .= 0.0
    R .= 0.0
    R[F,1] .= floatmax()
    L[:,end] .= llr
    nothing
end


