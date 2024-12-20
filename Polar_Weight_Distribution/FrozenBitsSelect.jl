"""
バタチャリア距離を用いた凍結ビットの選択法

FN:BP復号用
"""
function frozenbits_select_FN(N,K)    
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

"""
バタチャリア距離を用いた凍結ビットの選択法

GN:SC復号用

FNと何が違うか？
FN[i]はi番目の凍結ビットを表す
iをビット反転させた数をi_revとする
ビット反転とはなにか？
    例えば4を3ビットで表現する 普通にやると100となるが、
    1インデックス、つまり1,2,3,4,..と先頭
    0インデックスで計算するので1引いた数をビットで表現する
    4-1 = 3 で3を3ビットで表現すると 011になる
    そしてこれを逆側から読むのがビット反転である
    なので、011は110になる
"""
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

"""
Zvalue

"""
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

function Zvalue_GN(N)
    z =Zvalue(N)
    @. z = (bit_reversal_n(FN-1,Int64(log2(N))))+1
    z
end

struct ReedMuller; end
"""
selectmethodは凍結ビットの選択方法を示す

**Arguments**

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
