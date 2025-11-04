using DataFrames
import UncertaintyQuantification: evaluate!

struct HeadsModel <: UQModel
    f::Function
end

function evaluate!(m::HeadsModel, df::DataFrame)
    m.f(df)
end

function heads_values(df::DataFrame)
    v = mapreduce(x -> collect(range(10.323, 0, length=61))' .* x, vcat, df.head_factor)
    for (i, col) in enumerate(eachcol(v))
        df[!, Symbol("head_$i")] = col
    end
    return nothing
end


function head_formats(elements::Int=61)
    format = Dict([Symbol("head_" * string(i)) => ".4f" for i in range(1, elements)])
    return format
end

heads_model = HeadsModel(heads_values)