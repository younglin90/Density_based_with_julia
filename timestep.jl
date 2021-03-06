#= include("./module.jl")
include("./controls.jl")
using .mesh =#

function timestep!(
    ๐::controls, 
    cell::Vector{mesh.Cell},
    ฮx::Float64,
    ฮy::Float64,
    ฮz::Float64
)

    ๐.ฮt = 1.e10
    for i in cell

        U = sqrt(i.var[๐.u]^2+i.var[๐.v]^2+i.var[๐.w]^2)
        i.var[๐.Vแตฃ] = U
        i.var[๐.Vแตฃ] = max(i.var[๐.Vแตฃ], ๐.Lco/ฯ/๐.ฮt)
        i.var[๐.Vแตฃ] = min(i.var[๐.Vแตฃ], i.var[๐.c])



        i.var[๐.Vแตฃ] = i.var[๐.c]

        
		ฯต = (i.var[๐.Vแตฃ]/i.var[๐.c])^2.0
        ฮปu = 0.5*U*(1.0+ฯต)
        ฮปc = 0.5*โ(U^2*(1.0-ฯต)^2.0+4.0*i.var[๐.Vแตฃ]^2.0)

        ๐.ฮt = min(๐.ฮt,๐.CFL * min(ฮx,ฮy,ฮz) / (ฮปu+ฮปc))

    end
    
end


