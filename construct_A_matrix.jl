function construct_A_matrix_explicit!(
    ๐::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell}
)

    for i in 1:length(cells)

        ฯ = cells[i].var[๐.ฯ]
        u = cells[i].var[๐.u]
        v = cells[i].var[๐.v]
        w = cells[i].var[๐.w]
        T = cells[i].var[๐.T]
        Yโ = cells[i].var[๐.Yโ]
        Hโ = cells[i].var[๐.Hโ]
        โฯโp = cells[i].var[๐.โฯโp]
        โฯโT = cells[i].var[๐.โฯโT]
        โHโโp = cells[i].var[๐.โHโโp]
        โHโโT = cells[i].var[๐.โHโโT]
        โฯโYโ = cells[i].var[๐.โฯโYโ]
        โHโโYโ = cells[i].var[๐.โHโโYโ]

        T = zeros(Float64, 6, 6)
        T[1, 1] = โฯโp
        T[2, 1] = โฯโp*u
        T[3, 1] = โฯโp*v
        T[4, 1] = โฯโp*w
        T[5, 1] = โฯโp*Hโ + ฯ*โHโโp - 1.0
        
        T[2, 2] = ฯ
        T[5, 2] = ฯ*u
        
        T[3, 3] = ฯ
        T[5, 3] = ฯ*v
        
        T[4, 4] = ฯ
        T[5, 4] = ฯ*w
        
        T[1, 5] = โฯโT
        T[2, 5] = โฯโT*u
        T[3, 5] = โฯโT*v
        T[4, 5] = โฯโT*w
        T[5, 5] = โฯโT*Hโ + ฯ*โHโโT

        T[6, 1] = โฯโp*Yโ
        T[6, 5] = โฯโT*Yโ

        T[6, 6] = โฯโYโ*Yโ
        T[6, 6] += ฯ

        T[1, 6] = โฯโYโ
        T[2, 6] = โฯโYโ*u
        T[3, 6] = โฯโYโ*v
        T[4, 6] = โฯโYโ*w
        T[5, 6] = โฯโYโ*Hโ + ฯ*โHโโYโ

        P = deepcopy(T)
        ฮฒ = 1.0/cells[i].var[๐.Vแตฃ]^2.0
            + โฯโp - 1.0/cells[i].var[๐.c]^2.0
        P[1, 1] = ฮฒ
        P[2, 1] = ฮฒ*u
        P[3, 1] = ฮฒ*v
        P[4, 1] = ฮฒ*w
        P[5, 1] = ฮฒ*Hโ + ฯ*โHโโp - 1.0
        P[6, 1] = ฮฒ*Yโ

        A[i, :, :] = T[:, :]./๐.ฮt


    end


end



function construct_A_matrix_implicit!(
    ๐::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face}
)




end