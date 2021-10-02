function update_primitive!(
    👉::controls, 
    x::Array{Float64},
    cells::Vector{mesh.Cell},
    RK3::Int64
)

    if RK3 == 1
        for i in 1:length(cells)
            cells[i].var[👉.p] += x[i, 1]
            cells[i].var[👉.u] += x[i, 2]
            cells[i].var[👉.v] += x[i, 3]
            cells[i].var[👉.w] += x[i, 4]
            cells[i].var[👉.T] += x[i, 5]
            cells[i].var[👉.Y₁] += x[i, 6]
    
        end

    elseif RK3 == 2
        for i in 1:length(cells)
            cells[i].var[👉.p] = 3.0/4.0*cells[i].Qᵐ[👉.p] + 1.0/4.0*cells[i].var[👉.p] + 1.0/4.0*x[i, 1]
            cells[i].var[👉.u] = 3.0/4.0*cells[i].Qᵐ[👉.u] + 1.0/4.0*cells[i].var[👉.u] + 1.0/4.0*x[i, 2]
            cells[i].var[👉.v] = 3.0/4.0*cells[i].Qᵐ[👉.v] + 1.0/4.0*cells[i].var[👉.v] + 1.0/4.0*x[i, 3]
            cells[i].var[👉.w] = 3.0/4.0*cells[i].Qᵐ[👉.w] + 1.0/4.0*cells[i].var[👉.w] + 1.0/4.0*x[i, 4]
            cells[i].var[👉.T] = 3.0/4.0*cells[i].Qᵐ[👉.T] + 1.0/4.0*cells[i].var[👉.T] + 1.0/4.0*x[i, 5]
            cells[i].var[👉.Y₁] = 3.0/4.0*cells[i].Qᵐ[👉.Y₁] + 1.0/4.0*cells[i].var[👉.Y₁] + 1.0/4.0*x[i, 6]
    
        end

    else
        for i in 1:length(cells)
            cells[i].var[👉.p] = 1.0/3.0*cells[i].Qᵐ[👉.p] + 2.0/3.0*cells[i].var[👉.p] + 2.0/3.0*x[i, 1]
            cells[i].var[👉.u] = 1.0/3.0*cells[i].Qᵐ[👉.u] + 2.0/3.0*cells[i].var[👉.u] + 2.0/3.0*x[i, 2]
            cells[i].var[👉.v] = 1.0/3.0*cells[i].Qᵐ[👉.v] + 2.0/3.0*cells[i].var[👉.v] + 2.0/3.0*x[i, 3]
            cells[i].var[👉.w] = 1.0/3.0*cells[i].Qᵐ[👉.w] + 2.0/3.0*cells[i].var[👉.w] + 2.0/3.0*x[i, 4]
            cells[i].var[👉.T] = 1.0/3.0*cells[i].Qᵐ[👉.T] + 2.0/3.0*cells[i].var[👉.T] + 2.0/3.0*x[i, 5]
            cells[i].var[👉.Y₁] = 1.0/3.0*cells[i].Qᵐ[👉.Y₁] + 2.0/3.0*cells[i].var[👉.Y₁] + 2.0/3.0*x[i, 6]
    
        end

    end

end
