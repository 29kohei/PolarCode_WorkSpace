
function frozenbits_select_FN(N,K)    
    #=
    designSNR:
    =#
    designSNR = 0; n = Int(log2(N)); S = 10^(designSNR/10); 
    z = zeros(N)
    #z[1] = exp(-S)
    z[1] = 0.32 #design snr = 0.3
    #z[1]=0.36787944117144233
    for j=1:n
        u=2^j
        for t=1:u>>1
            T=z[t]
            z[t]=2*T-T^2
            z[u>>1+t]=T^2
        end
    end
    
    frozenbits = [1:N;]
    sort!(frozenbits,by=x->z[x],rev=true)
    for _=1:K; pop!(frozenbits); end
    return frozenbits
end

function frozenbits_select_GN(N,K)
    FN=frozenbits_select_FN(N,K)
    @. FN = (bit_reversal_n(FN-1,Int64(log2(N))))+1
end

function bit_reversal_n(b,n) 
    b = (b & 0xffffffff00000000) >> 32 | (b & 0x00000000ffffffff) << 32
    b = (b & 0xffff0000ffff0000) >> 16 | (b & 0x0000ffff0000ffff) << 16
    b = (b & 0xff00ff00ff00ff00) >> 8  | (b & 0x00ff00ff00ff00ff) << 8
    b = (b & 0xf0f0f0f0f0f0f0f0) >> 4  | (b & 0x0f0f0f0f0f0f0f0f) << 4
    b = (b & 0xcccccccccccccccc) >> 2  | (b & 0x3333333333333333) << 2
    b = (b & 0xaaaaaaaaaaaaaaaa) >> 1  | (b & 0x5555555555555555) << 1
    Int(b>>(64-n))
end

"""
frozenbits
"""
@inline infobits_select(frozenbits,N) = filter(p->!in(p,frozenbits),1:N)

#=
Z値を計算する関数
Arugument
N.符号長

Return
Z値が格納されたベクトル

Z値とはなにか？
チャンネルの誤り率である

SCの誤り率を計算するには？
=#
function Zvalue(N)
    designSNR = 0; n = Int(log2(N)); S = 10^(designSNR/10); 
    z = zeros(N)
    #z[1] = exp(-S)
    z[1] = 0.32 #design snr = 0.3
    #z[1]=0.36787944117144233
    for j=1:n
        u=2^j
        for t=1:u>>1
            T=z[t]
            z[t]=2*T-T^2
            z[u>>1+t]=T^2
        end
    end
    z
end

struct ReedMuller; end
"""
selectmethodは凍結ビットの選択方法を示す
#Arguments
selectmethod:凍結ビットの選択方法
N:符号長
K:符号長
r:最小距離
"""
function frozenbits_select_FN(selectmethod::Type{ReedMuller},N,K,r)
    z = Zvalue(N); n = Int(log2(N))
    G = kron_generate_G(n)
    index = filter(i -> r <= sum(G[i,:]),1:N)
    sort!(index,by=x->z[x],rev=true)
    [i for i=1:N if in(i,index)]
end

function frozenbits_select_GN(::Type{ReedMuller},N,K,r)
    FN = frozenbits_select_FN(RM,N,K,r)
    @. FN = (bit_reversal_n(FN-1,Int64(log2(N))))+1
end


#=
function rKRM(N,r)
    m = Int(log2(N))
    sum(binomial(m,i) for i=0:r)
end

function rRateRM(N,K)
    m = Int(log2(N))
    for i=0:m
        if rKRM(N,i) > K
            return i
        end
    end
end
=#