#=
=#

function MySet()
    set1 = Set([1,2,3,3,4])
    #Setに追加する
    push!(set1,1)

    #{6,7}
    set2 = Set([6,7])

    #和集合
    union(set1,set2) |> println

    #差集合
    setdiff(set1,set2) |> println

    issubset(set1,set2) |> println

    intersect(set1,set2) |> println
end