polar_encoder_FN!(seq) = polar_encoder!(seq,true)
polar_encoder_GN!(seq) = polar_encoder!(seq,false)

ope_xor(x::Integer,y::Integer)= xor(x,y)
ope_xor(x::T,y::T) where {T<:AbstractVector}= ifelse(isequal(x,y),zero(T),one(T)) 

function polar_encoder!(seq,V)
    N = length(seq); n = Int(log2(N))
    for i = (V ? (1:n) : (n:-1:1))
        kernel = N >> i
        @simd for j=1:N
            if iszero(mod(div(j-1,kernel),2))
                seq[j] = ope_xor(seq[j],seq[j+kernel])
            end
        end
    end
    nothing
end

"""

"""
function polar_parity_check(seq,frozenbits,cache_seq=similar(seq))
    @. cache_seq = seq
    polar_encoder_GN!(cache_seq)
    error = 0
    @simd for i in frozenbits
        error += (cache_seq[i]!=0)
    end
    (error == 0)
end

view_polar_infoseq(seq,ib) = @view seq[ib]

function kron_generate_G(n)
    i = 1
    kernel = [1 0; 1 1]
    G = [1]
    while i <= n
        G = kron(kernel,G)
       i += 1
    end
    G
end

function polar_systematic_encoder(seq,N,A)
    n = Int(log2(N))
    G = kron_generate_G(n)
    uA = mod.(G[A,A]*seq,2)
    u = zeros(eltype(uA),N)
    u[A] .= uA
    polar_encoder_FN!(u)
    u
end