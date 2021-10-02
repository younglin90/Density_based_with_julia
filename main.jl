include("./structured_grid_uniform.jl")
include("./constant.jl")
include("./controls.jl")
include("./timestep.jl")
include("./construct_A_matrix.jl")
include("./face_left_right.jl")
include("./flux.jl")
include("./linear_solver.jl")
include("./update_conservative.jl")
include("./update_primitive.jl")
include("./real_time_terms.jl")
include("./residual_norm.jl")
include("./EOS.jl")

using Plots

function main()


        #🎲
        #⬜
        #◽

    Nx = 100
    Ny = 1
    Nz = 1
    Lx = 1.0
    Ly = 1.0
    Lz = 0.5
    realMaxIter = 1000000
    pseudoMaxIter = 30
    pseudoMaxResidual = -4.0


    iterCFL0 = 0
    CFL0 = 0.0001
    CFL = 0.1
    Δt = 1.e-7
    Lco = 1.0
    Uco = 500.0

    👉 = controls(
        Nx,Ny,Nz, Lx,Ly,Lz, 
        realMaxIter,pseudoMaxIter,pseudoMaxResidual, 
        CFL, Δt, Lco, Uco,
        0.0, 0, 0, 0.0, 
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
    )

    cells = Vector{mesh.Cell}(undef, 0)
    faces = Vector{mesh.Face}(undef, 0)
    faces_internal = Vector{mesh.Face}(undef, 0)
    faces_boundary = Vector{mesh.Face}(undef, 0)
    faces_boundary_top = Vector{mesh.Face}(undef, 0)
    faces_boundary_bottom = Vector{mesh.Face}(undef, 0)
    faces_boundary_left = Vector{mesh.Face}(undef, 0)
    faces_boundary_right = Vector{mesh.Face}(undef, 0)

    structured_grid_uniform!(
        👉,
        cells,
        faces,
        faces_internal,
        faces_boundary,
        faces_boundary_top,
        faces_boundary_bottom,
        faces_boundary_left,
        faces_boundary_right
    )

    println(" Cell size = ",length(cells))
    println(" Face size = ",length(faces))
    println(" Face internal size = ",length(faces_internal))
    println(" Face boundary size = ",length(faces_boundary))
    println(" Face boundary top size = ",length(faces_boundary_top))
    println(" Face boundary bottom size = ",length(faces_boundary_bottom))
    println(" Face boundary left size = ",length(faces_boundary_left))
    println(" Face boundary right size = ",length(faces_boundary_right))

    # initialization
    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 0.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
    end
#=
    # 1D cavitation
    half_cell_num::Int32 = round(length(cells)/2)
    for i in 1:half_cell_num
        cells[i].var[👉.p] = 1.e6
        cells[i].var[👉.u] = -2500.0
        cells[i].var[👉.v] = 0.0
        cells[i].var[👉.w] = 0.0
        cells[i].var[👉.T] = 1000.0
        cells[i].var[👉.Y₁] = 0.0
    end
    
    for i in half_cell_num+1:length(cells)
        cells[i].var[👉.p] = 1.e6
        cells[i].var[👉.u] = 2500.0
        cells[i].var[👉.v] = 0.0
        cells[i].var[👉.w] = 0.0
        cells[i].var[👉.T] = 950.0
        cells[i].var[👉.Y₁] = 0.0
    end
=#

#=
    # 1D high pressure air & low pressure water
    for cell in cells
        if cell.x < 0.5
            cell.var[👉.p] = 1.e9
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 1.e-5
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 1.e4
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 1.0 - 1.e-5
            cell.var[👉.α₁] = 1.0
        end
    end
=#

#=

    # 1D high pressure water & low pressure air
    for cell in cells
        if cell.x < 0.5
            cell.var[👉.p] = 1.0
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.003484
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 0.1
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.002787
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end
    end

=#

    # One-dimensional helium-bubble in air
    for cell in cells
        
        if cell.x < 0.3
            cell.var[👉.p] = 1.245e5
            cell.var[👉.u] = 55.33
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 319.48
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 1.e5
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end

        if 0.5 < cell.x < 0.7
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end



    

    # EOS
    EOS!(👉, cells)


    # solver
    👉.time = 0.0
    👉.realIter = 1
    total_iter = 1.0
    
    plt2 = plot([],[])

    while(
        👉.realIter <= 👉.realMaxIter
    )
        println("real-time Step: $(👉.realIter) \t Time: $(👉.time)")

        if 👉.realIter < iterCFL0
            👉.CFL = CFL0
        else
            👉.CFL = CFL
        end
        
        for cell in cells
            cell.Qᵐ[1] = cell.var[👉.p]
            cell.Qᵐ[2] = cell.var[👉.u]#*cell.var[👉.u]
            cell.Qᵐ[3] = cell.var[👉.v]#*cell.var[👉.v]
            cell.Qᵐ[4] = cell.var[👉.w]#*cell.var[👉.w]
            cell.Qᵐ[5] = cell.var[👉.T]#*cell.var[👉.Hₜ]-cell.var[👉.p]
            cell.Qᵐ[6] = cell.var[👉.Y₁]#*cell.var[👉.Y₁]
        end


        for RK3 in 1:1

            # time-step
            timestep!(👉, cells, Lx/Nx, Ly/Ny, Lz/Nz)

            # face left, Right
            face_left_right!(👉, cells, faces, faces_internal, faces_boundary, 
                faces_boundary_top, faces_boundary_bottom, faces_boundary_left, faces_boundary_right)

            # right hand side
            B = zeros(Float64, length(cells), 6)

            # flux
            flux!(👉, B, cells, faces_internal, faces_boundary)

            #println(B)

            # sparse A matrix
            A = zeros(Float64, length(cells), 6, 6)
            #construct_A_matrix_implicit!(👉, A, cells, faces)
            construct_A_matrix_explicit!(👉, A, cells)

            # Ax=B
            x = zeros(Float64, length(cells), 6)
            #linear_solver_implicit!(A, x, B)
            linear_solver_explicit!(A, x, B)

            # update primitive
            update_primitive!(👉, x, cells, RK3)
            
            # EOS
            EOS!(👉, cells)
            
            gr()
            X = zeros(Float64, length(cells), 1)
            Y = zeros(Float64, length(cells), 6)
            for i in 1:length(cells)
                X[i] = cells[i].x
                Y[i,1] = cells[i].var[👉.p]
                Y[i,2] = cells[i].var[👉.u]
                Y[i,3] = cells[i].var[👉.T]
                Y[i,4] = cells[i].var[👉.Y₁]
                Y[i,5] = cells[i].var[👉.ρ]
                Y[i,6] = cells[i].var[👉.c]
                
            end
            #push!(plt2,total_iter,👉.residual-residual0)
            plt = plot(X,Y,layout = 
            grid(3, 2),
            label = ["p" "u" "T" "Y₁" "ρ" "c"] )
            plot(plt,plt2,layout = 
            grid(2, 1, heights=[0.8 ,0.2]))

            gui()
        end


        👉.realIter += 1
        👉.time += 👉.Δt

    end





    #Δt = CFL * Δx
end



# calculation main
main()


