dir = "C:\\Users\\meph_\\Kohei\\VScode\\polarcode\\PolarCode_WorkSpace\\Polar_Code_Utility"
include("$(dir)/IndexMatrix.jl")
include("$(dir)/BP.jl")
include("$(dir)/AWGN.jl")
include("$(dir)/Polar.jl")
#=

=#
struct BPProblem
    N::Int
    K::Int
    mi::Int
    M::Matrix{Int}
    L::Matrix{Float64}
    R::Matrix{Float64}
    polar::Polar
end

BPProblem(N,K,mi)  = begin 
    M = indexMatrix(N)
    n = Int(log2(N))
    L = zeros(N,n+1)
    R = zeros(N,n+1)
    polar = Polar(Batt,N,K)
    BPProblem(N,K,mi,M,L,R,polar)
end

#=
LLRを新しくする

=#
function initialize!(L,R,F,ch::AWGN)
    updatellr!(ch)
    initializeLandR!(L,R,F,ch.llr)
    nothing
end

function iserror(t::BPProblem)
    A = getinfobits(t.polar)
    @views !isnothing(findfirst(<(0.0),t.L[A,1]))
end

function errors(t::BPProblem)
    A = getinfobits(t.polar)
    sum(ifelse(t.L[i,1] <= 0.0,1,0) for i in A) 
end

function solver(sim,t::BPProblem,ch::AWGN)
    n = Int(log2(t.N)); F = getfrozenbits(t.polar); A = getinfobits(t.polar)
    u = zeros(Int,t.N); v = zeros(Int,t.N)
    
    #=
    BP復号の終了判定を行う関数
    =#
    function checkcond()
        for i=1:t.N
            #=
            harddecisionの
            =#
            v[i] = harddecision(Int,t.L[i,end],t.R[i,end])
            u[i] = harddecision(Int,t.L[i,1],t.R[i,1])
        end 
        Gcheck(u,v)
    end

    bler = 0; ber = 0; errornum = zeros(Int,sim)
    
    for i=1:sim
        #=
        初期化
        =#
        initialize!(t.L,t.R,F,ch)
        #=
        BP復号
        =#
        success = BP(t.N,t.mi,t.L,t.R,t.M,checkcond)
        
        #=
        =#
        if !iserror(t)
            continue
        end
        

        errornum[i] = errors(t)
        bler += 1
        ber += errornum[i]
    end
    return bler/sim,ber/(sim*t.K), errornum
end