using Random
struct AWGN
    N :: Int
    K :: Int
    eb::Float64
    snrdb :: Float64
    R :: Float64
    rng :: MersenneTwister
    llr ::Vector{Float64}

    AWGN(N,K,seed,eb) = begin 
        R = K/N
        snrdb = eb + 10.0*log10(R)
        rng = MersenneTwister(seed)
        llr = zeros(N)
        new(N,K,eb,snrdb,R,rng,llr)
    end
end

function updatellr!(ch::AWGN)
    llr = ch.llr; rng = ch.rng; snrdb = ch.snrdb
    for i in eachindex(llr)
        y = -10^(snrdb/20) + sqrt(1/2)*randn(rng)
        llr[i] = -4*y*(10^(snrdb/20))
    end
    nothing
end

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

struct Polar
    F::Vector{Int}
    A::Vector{Int}
end
getfrozenbits(p::Polar) = p.F
getinfobits(p::Polar) = p.A

struct SCLProblem
    N::Int
    K::Int
    pm::PathManager
    p::Polar
end
#=
SCLProblem(N,K,L) = begin
    p = Polar(BattGN,N,K)
    pm = PathManager(N,L)
    SCLProblem(N,K,pm,p)
end
=#

SCLProblem_eBCH_16_7(L) = begin
    F = [1,2,3,5,9,7,10,11,13]; sort!(F)
    A = filter(x->!in(x,F),1:16); sort!(A)
    p = Polar(F,A)
    pm = PathManager(N,L)
    SCLProblem(16,7,p,pm)
end

@inline function initialize!(t::SCLProblem,a::AWGN)
    updatellr!(a)
    initialize!(t.pm,a.llr)
    nothing
end

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

@inline function SCL(pm::PathManager,isfrozen,dynamic_frozen)
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
                path_l = getpathpath(pm,l)
                dynamic_frozen_value = dynamic_frozen(i,path_l)
                C = getArrayPointer_C(pm,layer,l)
                P = getArrayPointer_P(pm,layer,l)
                #C[1][phasetoindex(i)] = frozenbits_value 
                C[1][phasetoindex(i)] = dynamic_frozen_value
                #pm.pathmetric[l] = pathmetricaprox(pm.pathmetric[l],P[1],frozenbits_value)
                pm.pathmetric[l] = pathmetricaprox(pm.pathmetric[l],P[1],dynamic_frozen_value)
                #path_l = getpathpath(pm,l)
                #path_l[i] = frozenbits_value
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

function dynamic_eBCH_16_7_frozen(i,path)
    if i==7
        path[7] = path[4]
    elseif i==10
        path[10] = path[6]
    elseif i==11
        path[11] = mod(path[4]+path[6],0:1)
    elseif i==13
        path[13] = path[11]
    else
        #default
        path[i] = 0
    end
    return path[i]
end

struct AWGN
    N :: Int
    K :: Int
    eb::Float64
    snrdb :: Float64
    R :: Float64
    rng :: MersenneTwister
    llr ::Vector{Float64}

    AWGN(N,K,seed,eb) = begin 
        R = K/N
        snrdb = eb + 10.0*log10(R)
        rng = MersenneTwister(seed)
        llr = zeros(N)
        new(N,K,eb,snrdb,R,rng,llr)
    end
end

function updatellr!(ch::AWGN)
    llr = ch.llr; rng = ch.rng; snrdb = ch.snrdb
    for i in eachindex(llr)
        y = -10^(snrdb/20) + sqrt(1/2)*randn(rng)
        llr[i] = -4*y*(10^(snrdb/20))
    end
    nothing
end


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

struct Polar
    F::Vector{Int}
    A::Vector{Int}
end
getfrozenbits(p::Polar) = p.F
getinfobits(p::Polar) = p.A

struct SCLProblem
    N::Int
    K::Int
    pm::PathManager
    p::Polar
end
#=
SCLProblem(N,K,L) = begin
    p = Polar(BattGN,N,K)
    pm = PathManager(N,L)
    SCLProblem(N,K,pm,p)
end
=#

SCLProblem_eBCH_16_7(L) = begin
    F = [1,2,3,5,9,7,10,11,13]; sort!(F)
    A = filter(x->!in(x,F),1:16); sort!(A)
    p = Polar(F,A)
    pm = PathManager(16,L)
    SCLProblem(16,7,pm,p)
end

@inline function initialize!(t::SCLProblem,a::AWGN)
    updatellr!(a)
    initialize!(t.pm,a.llr)
    nothing
end

function solver(sim,t::SCLProblem,a::AWGN)
    v = [in(i,t.p.F) for i=1:t.N]
    @inline isfrozen(i)  = v[i]
    bler = 0; ber = 0
    for i=1:sim
        initialize!(t,a)
        p = SCL(t.pm,isfrozen,dynamic_eBCH_16_7_frozen)
        tmp = errors(t.pm,p)
        if !iszero(tmp)
            bler += 1
            ber += tmp
        end
    end
    return bler/sim,ber/sim/t.K
end

function main(sim,ebno,L)
    t = SCLProblem_eBCH_16_7(L)
    ch = AWGN(16,7,8,ebno)
    res = solver(sim,t,ch)
    open("result.txt","a") do file
        println(file,"ebno,$ebno")
        println(file,"L,$L")
        println(file,"N,16")
        println(file,"K,7")
        println(file,res)
    end
end





