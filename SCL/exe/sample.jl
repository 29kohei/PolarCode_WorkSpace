#=
2024/12/06
依存関係

=#
function sample()
    if length(ARGS) > 0
        for i in eachindex(ARGS)
            sim = parse(Int,ARGS[i])
            println("sim $sim")
        end
    else
        println("No arguments provided")
    end
end

sample()