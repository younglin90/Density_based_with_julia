function update_primitive!(
    ๐::controls, 
    x::Array{Float64},
    cells::Vector{mesh.Cell},
    RK3::Int64
)

    if RK3 == 1
        for i in 1:length(cells)
            cells[i].var[๐.p] += x[i, 1]
            cells[i].var[๐.u] += x[i, 2]
            cells[i].var[๐.v] += x[i, 3]
            cells[i].var[๐.w] += x[i, 4]
            cells[i].var[๐.T] += x[i, 5]
            cells[i].var[๐.Yโ] += x[i, 6]
    
        end

    elseif RK3 == 2
        for i in 1:length(cells)
            cells[i].var[๐.p] = 3.0/4.0*cells[i].Qแต[๐.p] + 1.0/4.0*cells[i].var[๐.p] + 1.0/4.0*x[i, 1]
            cells[i].var[๐.u] = 3.0/4.0*cells[i].Qแต[๐.u] + 1.0/4.0*cells[i].var[๐.u] + 1.0/4.0*x[i, 2]
            cells[i].var[๐.v] = 3.0/4.0*cells[i].Qแต[๐.v] + 1.0/4.0*cells[i].var[๐.v] + 1.0/4.0*x[i, 3]
            cells[i].var[๐.w] = 3.0/4.0*cells[i].Qแต[๐.w] + 1.0/4.0*cells[i].var[๐.w] + 1.0/4.0*x[i, 4]
            cells[i].var[๐.T] = 3.0/4.0*cells[i].Qแต[๐.T] + 1.0/4.0*cells[i].var[๐.T] + 1.0/4.0*x[i, 5]
            cells[i].var[๐.Yโ] = 3.0/4.0*cells[i].Qแต[๐.Yโ] + 1.0/4.0*cells[i].var[๐.Yโ] + 1.0/4.0*x[i, 6]
    
        end

    else
        for i in 1:length(cells)
            cells[i].var[๐.p] = 1.0/3.0*cells[i].Qแต[๐.p] + 2.0/3.0*cells[i].var[๐.p] + 2.0/3.0*x[i, 1]
            cells[i].var[๐.u] = 1.0/3.0*cells[i].Qแต[๐.u] + 2.0/3.0*cells[i].var[๐.u] + 2.0/3.0*x[i, 2]
            cells[i].var[๐.v] = 1.0/3.0*cells[i].Qแต[๐.v] + 2.0/3.0*cells[i].var[๐.v] + 2.0/3.0*x[i, 3]
            cells[i].var[๐.w] = 1.0/3.0*cells[i].Qแต[๐.w] + 2.0/3.0*cells[i].var[๐.w] + 2.0/3.0*x[i, 4]
            cells[i].var[๐.T] = 1.0/3.0*cells[i].Qแต[๐.T] + 2.0/3.0*cells[i].var[๐.T] + 2.0/3.0*x[i, 5]
            cells[i].var[๐.Yโ] = 1.0/3.0*cells[i].Qแต[๐.Yโ] + 2.0/3.0*cells[i].var[๐.Yโ] + 2.0/3.0*x[i, 6]
    
        end

    end

end
