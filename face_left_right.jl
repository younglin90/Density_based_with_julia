function face_left_right!(
    π::controls, 
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face}
)

    # internal faces
    for face in faces_internal

        face.varβ[π.p] = cells[face.owner].var[π.p]
        face.varα΅£[π.p] = cells[face.neighbour].var[π.p]

        face.varβ[π.u] = cells[face.owner].var[π.u]
        face.varα΅£[π.u] = cells[face.neighbour].var[π.u]
        
        face.varβ[π.v] = cells[face.owner].var[π.v]
        face.varα΅£[π.v] = cells[face.neighbour].var[π.v]
        
        face.varβ[π.w] = cells[face.owner].var[π.w]
        face.varα΅£[π.w] = cells[face.neighbour].var[π.w]
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = cells[face.neighbour].var[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = cells[face.neighbour].var[π.Yβ]

    end


    # Top
    #bc_sub_outlet!(π, faces_boundary_top, cells, 101325.0)
    #bc_moving_wall!(π, faces_boundary_top, cells, 1.0, 0.0, 0.0)
    bc_slip_wall!(π, faces_boundary_top, cells)

    # Bottom
    #bc_wall!(π, faces_boundary_bottom, cells)
    bc_slip_wall!(π, faces_boundary_bottom, cells)

    # Left
    #bc_slip_wall!(π, faces_boundary_left, cells)
    #bc_wall!(π, faces_boundary_left, cells)
    bc_sup_outlet!(π, faces_boundary_left, cells)

    # Right
    #bc_slip_wall!(π, faces_boundary_right, cells)
    #bc_wall!(π, faces_boundary_right, cells)
    bc_sup_outlet!(π, faces_boundary_right, cells)


    # EOS
    for face in faces
        face.varβ[π.Ο], face.varβ[π.Hβ], face.varβ[π.c] = 
        faceEOS!(face.varβ[π.p], 
        face.varβ[π.u], face.varβ[π.v], face.varβ[π.w], face.varβ[π.T], face.varβ[π.Yβ])

        
        face.varα΅£[π.Ο], face.varα΅£[π.Hβ], face.varα΅£[π.c] = 
        faceEOS!(face.varα΅£[π.p], 
        face.varα΅£[π.u], face.varα΅£[π.v], face.varα΅£[π.w], face.varα΅£[π.T], face.varα΅£[π.Yβ])

    end


end


function bc_wall!(
    π::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)
    for face in faces_boundary

        face.varβ[π.p] = cells[face.owner].var[π.p]
        face.varα΅£[π.p] = face.varβ[π.p]

        face.varβ[π.u] = 0.0
        face.varα΅£[π.u] = 0.0
        
        face.varβ[π.v] = 0.0
        face.varα΅£[π.v] = 0.0
        
        face.varβ[π.w] = 0.0
        face.varα΅£[π.w] = 0.0
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = face.varβ[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = face.varβ[π.Yβ]

    end

end

function bc_slip_wall!(
    π::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)

    for face in faces_boundary

        face.varβ[π.p] = cells[face.owner].var[π.p]
        face.varα΅£[π.p] = face.varβ[π.p]

        Un = 0.0
        Un += cells[face.owner].var[π.u]*face.nΜ[1]
        Un += cells[face.owner].var[π.v]*face.nΜ[2]
        Un += cells[face.owner].var[π.w]*face.nΜ[3]

        invU = cells[face.owner].var[π.u] - Un * face.nΜ[1]
        invV = cells[face.owner].var[π.v] - Un * face.nΜ[2]
        invW = cells[face.owner].var[π.w] - Un * face.nΜ[3]

        face.varβ[π.u] = invU
        face.varα΅£[π.u] = invU
        
        face.varβ[π.v] = invV
        face.varα΅£[π.v] = invV
        
        face.varβ[π.w] = invW
        face.varα΅£[π.w] = invW
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = face.varβ[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = face.varβ[π.Yβ]

    end

end


    
function bc_moving_wall!(
    π::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell},
    u::Float64, v::Float64, w::Float64
)

    for face in faces_boundary

        face.varβ[π.p] = cells[face.owner].var[π.p]
        face.varα΅£[π.p] = face.varβ[π.p]

        face.varβ[π.u] = u
        face.varα΅£[π.u] = u
        
        face.varβ[π.v] = v
        face.varα΅£[π.v] = v
        
        face.varβ[π.w] = w
        face.varα΅£[π.w] = w
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = face.varβ[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = face.varβ[π.Yβ]

    end

end


function bc_sub_outlet!(
    π::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell},
    set_p::Float64
)

    for face in faces_boundary

        face.varβ[π.p] = set_p
        face.varα΅£[π.p] = set_p

        face.varβ[π.u] = cells[face.owner].var[π.u]
        face.varα΅£[π.u] = face.varβ[π.u]
        
        face.varβ[π.v] = cells[face.owner].var[π.v]
        face.varα΅£[π.v] = face.varβ[π.v]
        
        face.varβ[π.w] = cells[face.owner].var[π.w]
        face.varα΅£[π.w] = face.varβ[π.w]
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = face.varβ[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = face.varβ[π.Yβ]

    end


end
    


function bc_sup_outlet!(
    π::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)

    for face in faces_boundary

        face.varβ[π.p] = cells[face.owner].var[π.p]
        face.varα΅£[π.p] = face.varβ[π.p]

        face.varβ[π.u] = cells[face.owner].var[π.u]
        face.varα΅£[π.u] = face.varβ[π.u]
        
        face.varβ[π.v] = cells[face.owner].var[π.v]
        face.varα΅£[π.v] = face.varβ[π.v]
        
        face.varβ[π.w] = cells[face.owner].var[π.w]
        face.varα΅£[π.w] = face.varβ[π.w]
        
        face.varβ[π.T] = cells[face.owner].var[π.T]
        face.varα΅£[π.T] = face.varβ[π.T]
        
        face.varβ[π.Yβ] = cells[face.owner].var[π.Yβ]
        face.varα΅£[π.Yβ] = face.varβ[π.Yβ]

    end


end
    

