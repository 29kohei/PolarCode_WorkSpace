#=
issue fill
fillは使わない
=#

#=
確率用にPathを変更する
=#

struct Path
    pathindex::Int
    P::Vector{Vector{Vector{Float64}}}
    C::Vector{Vector{Vector{Int}}}
    path::Vector{Int}

    #=
    i,pathのインデックス
    N,符号長
    =#
    Path(i::Int,N::Int) =
    begin
        l = Int(log2(N)) + 1
        #P = map(i -> zeros(Float64,2^(i-1)), l:-1:1)
        P = map(i -> [[zero(Float64),zero(Float64)] for j=1:2^(i-1)], l:-1:1)
        C = map(i->[[zero(Int),zero(Int)] for j=1:2^(i-1)],l:-1:1)
        path = [-1 for _=1:N]
        new(i,P,C,path)
    end
end

function reset!(p::Path)
    P = p.P
    C = p.C

    #=
    PはVector{Vecotr{Vector{Float64}}}
    =#
    
    for x in eachindex(P)
        for y in eachindex(P[x])
            P[x][y] .= zero(eltype(P[x][y]))
            C[x][y] .= zero(eltype(C[x][y]))
        end
    end
    p.path .= -1
    nothing
end

#=
LLR
=#

struct PathManager
    patharr::Vector{Path}
    activePath::Vector{Bool}
    probForks::Matrix{Float64}
    contForks::Matrix{Bool}
    index::Vector{Int}
    pathmetric::Vector{Float64}
    _L::Int
    _N::Int
    _l::Int

    PathManager(N,L) = 
    begin
        patharr = map(i->Path(i,N),1:L)
        activePath = zeros(Bool,L)
        probForks = zeros(2,L)
        contForks = zeros(Bool,2,L)
        index = [1:2L;]
        pathmetric = zeros(L)
        layer = Int(log2(N)) + 1
        new(patharr,activePath,probForks,contForks,index,pathmetric,L,N,layer)
    end
end

function reset!(pm::PathManager) 
    L = getListSize(pm)
    for i in pm.patharr
        reset!(i)
    end
    pm.activePath .= false
    pm.probForks .= zero(eltype(pm.probForks))
    pm.contForks .= zero(eltype(pm.contForks))
    pm.pathmetric .= zero(eltype(pm.pathmetric))
    nothing
end

#=
LLR用の初期化
=#
function initialize!(pm::PathManager,llr)
    reset!(pm)
    activepath(pm,1)
    patharr = pm.patharr
    for i in eachindex(llr)
        x = exp(llr[i])
        p0 = x/(x+1.0)
        p1 = 1.0/(x+1.0)
        patharr[1].P[1][i][begin] = p0
        patharr[1].P[1][i][end] = p1
    end
    nothing
end

@inline inactivepathindex(pm::PathManager) = findfirst(!,pm.activePath)
@inline activepath(pm::PathManager,i) = (pm.activePath[i] = true)
@inline inactivepath(pm::PathManager,i) = (pm.activePath[i] = false)
@inline isactive(pm::PathManager,l) = pm.activePath[l]
@inline getPath(pm::PathManager,i)  = pm.patharr[i]
@inline getListSize(pm::PathManager) = pm._L
@inline getlayer(pm::PathManager) = pm._l
@inline killPath(pm::PathManager,i) = !inactivepath(pm,i)
@inline getpathpath(pm::PathManager,l) = getPath(pm,l).path
@inline getArrayPointer_C(pm::PathManager,layer,l) = pm.patharr[l].C[layer]
@inline getArrayPointer_P(pm::PathManager,layer,l) = pm.patharr[l].P[layer]
@inline getpathmetric(pm::PathManager,l) = pm.pathmetric[l]

#dにsをcopyする
@inline function copyPath(s::Path,d::Path)
    #path
    d.path .= s.path
    #P
    for k in eachindex(s.P)
        for l in eachindex(s.P[k])
            d.P[k][l] .= s.P[k][l]
        end
    end
    #C
    for k in eachindex(s.C)
        for l in eachindex(s.C[k])
            d.C[k][l] .= s.C[k][l]
        end
    end
    nothing
end

@inline function copyPath(pm::PathManager,i,j)
    s = getPath(pm,i); d = getPath(pm,j)
    copyPath(s,d)
    nothing
end

function clonePath(pm::PathManager,i)
    j = inactivepathindex(pm)
    activepath(pm,j)
    copyPath(pm,i,j)
    j
end


@inline phasetophai(phase,::Type{T}=Int) where {T<:Integer} = ceil(T,phase/2) 
@inline phasetoindex(phase) = mod(phase-1,2)+1
#=
#LLR用のf,g関数

@inline f(a,b) = min(abs(a),abs(b))sign(a)sign(b)
@inline g(a,b,u) = a*(-1)^u + b
@inline pathmetricaprox(m::Float64,k::Float64,u) = (u == (1-sign(k))/2) ? m : (m+abs(k))
=#

#=
a = P[layer][phase]
b = P[layer][phase]
=#

#確率用のf,g関数
"""
a [W(|0),W(|1)]
b [W(|0),W(|1)]
c 0,1

aに0.0が一つでもあるとどうなるか？

"""
@inline function f(a,b,c)
    w = zero(Float64)
    for i=0:1
        u = mod(i + c,2)
        ak = u + 1
        bk = i + 1
        w += a[ak]*b[bk]
    end
    1/2*w
end


"""
a [W(|0),W(|1)]
b [W(|0),W(|1)]
u1 0,1
u2 0,1
"""
@inline function g(a::Vector{Float64},b::Vector{Float64},u1,u2)
    ak = mod(u1+u2,0:1)+1
    bk = u2 + 1
    1/2*a[ak]*b[bk]
end

"""
SCL復号の誤りの個数

パスで、0でないビット（値が1であるビット）の総数
"""
@inline errors(pm::PathManager,i) = sum(getpathpath(pm,i))