using Parameters
struct ProductTurboProblem
    Imax::Int
    Lapp::Matrix{Float64}
    Lch::Matrix{Float64}
    La::Matrix{Float64}
    Le::Matrix{Float64}
end

struct Row; end
struct Column; end

function solver(sim,productturboproblem::ProductTurboProblem)
    @unpack Imax,Lapp,Lch,La,Le = productturboproblem
    for i=1:sim
        
    end
end