#=
依存関係

=#
function indexMatrix(N)
    x = [1:N;]
    n = Int(log2(N))
    M = zeros(eltype(x),N-1, n)
    for k = 0:n-1
        for i = 0:1<<(k+1):N
            if (1<<k+i) < N
                for j=1:2^k; M[i+j,n-k] = x[i+j]; end;
            end
        end
    end
    A = reshape(M[M.>0],N>>1,n)
end


