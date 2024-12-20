include("FrozenBitsSelect.jl")

"""
凍結ビットと情報ビットを格納する構造体
"""
struct Polar
    F::Vector{Int}
    A::Vector{Int}
end
getfrozenbits(p::Polar) = p.F
getinfobits(p::Polar) = p.A

abstract type Batt end

Polar(::Type{Batt},N,K) = 
begin
    F = frozenbits_select_FN(N,K); sort!(F)
    A = infobits_select(F,N); sort!(A)
    Polar(F,A)
end

#=
Polar(::Type{RM},N,K,r) = 
begin
    F = frozenbits_select_FN(RM,N,K,r); sort!(F)
    A = infobits_select(F,N); sort!(A)
    Polar(F,A)
end
=#

function ()
end