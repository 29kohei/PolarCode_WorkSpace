#=
issue fill
fillは使わない
=#
struct Path
    pathindex::Int
    P::Vector{Vector{Float64}}
    C::Vector{Vector{Vector{Int}}}
    path::Vector{Int}

    Path(i::Int,N::Int) =
    begin
        l = Int(log2(N)) + 1
        P = map(i -> zeros(Float64,2^(i-1)), l:-1:1)
        C = map(i->[[zero(Int),zero(Int)] for j=1:2^(i-1)],l:-1:1)
        path = [-1 for _=1:N]
        new(i,P,C,path)
    end
end

function reset!(p::Path)
    P = p.P
    C = p.C
    for x in eachindex(p.P) 
        P[x] .= zero(eltype(P[x]))
        for y in eachindex(C[x])
            C[x][y] .= zero(eltype(C[x][y]))
        end
    end
    p.path .= -1
    nothing
end

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

function initialize!(pm::PathManager,llr)
    reset!(pm)
    activepath(pm,1)
    patharr = pm.patharr
    patharr[1].P[1] .= llr
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
        d.P[k] .= s.P[k]
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
@inline f(a,b) = min(abs(a),abs(b))sign(a)sign(b)
@inline g(a,b,u) = a*(-1)^u + b
@inline pathmetricaprox(m::Float64,k::Float64,u) = (u == (1-sign(k))/2) ? m : (m+abs(k))

@inline errors(pm::PathManager,i) = sum(getpathpath(pm,i))