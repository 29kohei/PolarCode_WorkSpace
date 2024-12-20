#=
Reed Muller符号の
=#
function calcK(m,r)
    comb(n,k) = prod(1:n)/prod(1:k)/prod(1:n-k)
    sum(comb(m,i) for i=0:r)
end

"""
Reed Muller符号の最小距離を調べる

m:2^mの符号長

r:monomialの多項式の最大次数
"""
calc_mindistance(m,r) = 1<<(m-r)

dat = "dat"
function RMcode_Output_File(N)
    m = Int(log2(N))
    l = [(Int(calcK(m,r)), calc_mindistance(m,r)) for r=0:m]

    string_K_mindistance = String[]
    #for文の中がlに依存している
    for i=0:m
        j = i+1
        K = l[j][1]
        minimumdistance = l[j][2]
        entry = "($K,$(minimumdistance))"
        push!(string_K_mindistance,entry)
    end
 
    table_K_minimumdistance = ""
    for i=0:m
        j = i+1
        entry = string_K_mindistance[j]
        table_K_minimumdistance = string(table_K_minimumdistance,entry,"\n")
    end

    output = """
    Reed Muller Codeの最小距離と情報ビット数(K)を符号長($N)から計算
    (K,minimum distance) 
    $(table_K_minimumdistance)
    """
    
    open("$(dat)/reedmuller$N.dat","w") do file
        println(file,output)
    end
    println(output)
    nothing
end