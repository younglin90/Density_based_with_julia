#= include("./module.jl")
include("./controls.jl")
using .mesh =#

function timestep!(
    👉::controls, 
    cell::Vector{mesh.Cell},
    Δx::Float64,
    Δy::Float64,
    Δz::Float64
)

    👉.Δt = 1.e10
    for i in cell

        U = sqrt(i.var[👉.u]^2+i.var[👉.v]^2+i.var[👉.w]^2)
        i.var[👉.Vᵣ] = U
        i.var[👉.Vᵣ] = max(i.var[👉.Vᵣ], 👉.Lco/π/👉.Δt)
        i.var[👉.Vᵣ] = min(i.var[👉.Vᵣ], i.var[👉.c])



        i.var[👉.Vᵣ] = i.var[👉.c]

        
		ϵ = (i.var[👉.Vᵣ]/i.var[👉.c])^2.0
        λu = 0.5*U*(1.0+ϵ)
        λc = 0.5*√(U^2*(1.0-ϵ)^2.0+4.0*i.var[👉.Vᵣ]^2.0)

        👉.Δt = min(👉.Δt,👉.CFL * min(Δx,Δy,Δz) / (λu+λc))

    end
    
end


