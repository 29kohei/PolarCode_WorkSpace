#=
include("BPProblem.jl")
include("/AWGN.jl")
include("../BP.jl")
=#


#=
checkcond,終了判定をする返す true or falseを返す
=#



function solver(sim,t::BPProblem,ch::AWGN)
    #N=2^n, F:凍結ビットの集合 Vector{Int}, A:情報ビットの集合 Vector{Int}
    n = Int(log2(t.N)); F = getfrozenbits(t.polar); A = getinfobits(t.polar)
    
    #=
    Gcheckを使う
    =#
    u = zeros(Int,t.N); v = zeros(Int,t.N)
    function checkcond()
        for i=1:t.N
            v[i] = harddecision(Int,t.L[i,end],t.R[i,end])
            u[i] = harddecision(Int,t.L[i,1],t.R[i,1])
        end 
        Gcheck(u,v)
    end

    bler = 0; ber = 0; errornum = zeros(Int,sim)
    
    for i=1:sim
        initialize!(t.L,t.R,F,ch)
        success = BP(t.N,t.mi,t.L,t.R,t.M,checkcond)
        if !iserror(t)
            continue
        end
        errornum[i] = errors(t)
        bler += 1
        ber += errornum[i]
    end
    return bler/sim,ber/(sim*t.K), errornum
end