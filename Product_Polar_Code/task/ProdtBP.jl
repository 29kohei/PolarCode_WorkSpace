using SQLite
using Random
include("../../Encoder/PolarEncode.jl")
include("../../IndexMatrix.jl")
include("../FrozenBitsSelect.jl")

struct Channel
    N :: Int
    K :: Int
    eb::Float64
    snrdb :: Float64
    R :: Float64
    rng :: MersenneTwister
    LLR ::Vector{Float64}

    Channel(N,K,seed,eb) = begin 
        R = K/N
        snrdb = ebtosnrdb(eb,R)
        rng = MersenneTwister(seed)
        LLR  = zeros(N)
        new(N,K,eb,snrdb,R,rng,LLR)
    end
end

struct PolarProblem
    F::Vector{Int}
    A::Vector{Int}
end

struct ProductBPProblem
    N :: Int
    K :: Int
    L :: Matrix{Float64}
    R :: Matrix{Float64}
    M :: Matrix{Int}
    sim::Int
    Imax::Int
    alpha :: Float64
    eb :: Float64
    Le :: Matrix{Float64}
    La :: Matrix{Float64}
    Lapp :: Matrix{Float64}
    Lch_plus_a :: Matrix{Float64}
    ProductBP(N,K,sim,Imax,alpha,eb) = begin
        n = Int(log2(N))
        L = zeros(N,n+1)
        R = zeros(N,n+1)
        M = indexMatrix(N)
        Le = zeros(N,N)
        La = zeros(N,N)
        Lapp = zeros(N,N)
        Lch_plus_a = zeros(N,N)
        new(N,K,L,R,M,sim,Imax,alpha,eb,Le,La,Lapp,Lch_plus_a)
    end
end


function solver(product_bp_problem :: ProductBPProblem, ch::Channel,polar::PolarProblem)
    @unpack N1,K1,L,R,alphaN,K,L,R,M,sim,Imax,alpha,eb,Le,La,Lapp,Lch_plus_a = product_bp_problem
    @unpack N,K,
    @unpack F,A = polar
    if N == N1^2; println("not product code N is not equal N1^2"); return -2.0; end
    n = Int(log2(N))

    cache_checkcondforBP = zeros(Int,N)
    function checkcondforBP()
        @views harddecision(cache_checkcondforBP,L[:,n+1])
        polar_encoder_FN!(cache_checkcondforBP)
        @views if hardcheck(R[:,1],cache_checkcondforBP)
            return true
        else 
            return false
        end
    end

    rec = fill(false,N); bler = 0
    for t = 1:sim
        initializeLLR!(rng,snrdb,llr)
        for _=1:Imax
            #row
            for i=1:N
                @views initializeLandR!(L,R,F,Lch_plus_a[i,:])
                res = BP(N,mi,L,R,M,checkcond); rec[i] = res
                @views updateLapp!(i,0,Lapp,L,R)
                @views updateLeandLa!(i,0,Lapp,Lch,La,Le,alpha)
            end

            if isnothing(findfirst(isequal(false),rec))
                break
            end

            updateLch_plus_La!(Lch_plus_a,Lch,La)
            #column
            for j=1:N
                @views initializeLandR!(L,R,F,Lch_plus_a[:,j])
                res = BP(N,mi,L,R,M,checkcond); rec[i] = res
                @views updateLapp!(0,j,Lapp,L,R)
                @views updateLeandLa!(0,j,Lapp,Lch,La,Le,alpha)
            end

            if isnothing(findfirst(isequal(false),rec))
                break
            end
        end

        if isbperror(Lapp); bler+=1 ; end
        
    end
    return bler
end

function main(::Type{ProductBPProblem};N,K,sim,alpha)
    product_code_N = Int(sqrt(N));product_code_K = Int(sqrt(K))
    pbp = PruductBP(product_code_N,product_code_K,sim,alpha,eb)
    ch = Channel(N,K,seed,eb)
    F = frozenbits_select_FN(product_code_N,product_code_K); A = infobits_select(F,product_code_N)
    polar = Polar(F,A)
    bler = solver(pbp,ch,polar)

    #DBに書き込み
    dn = "ProductBP"
    db = SQLite.DB("Simulation_Result.db")
    td = DataFrame(SQLite.tables(db))
    if !in(dn,td[:,:name])
        DBInterface.execute(db,"create table $dn(N,K,sim,alpha,eb,bler)")
    end
    df = DataFrame(DBInterface.execute(db,"select * from $dn where N=$N and K=$K and sim=$sim and alpha=$alpha and eb = $eb"))
    if size(df)[1] == 0
        DBInterface.execute(db,"insert into $dn values($N,$K,$sim,$alpha,$eb,$bler)")
    else 
        DBInterface.execute(db,"update $dn set bler=$bler where N=$N and K=$K and sim=$sim and alpha=$alpha and eb = $eb")
    end
    SQLite.close(db)
end

