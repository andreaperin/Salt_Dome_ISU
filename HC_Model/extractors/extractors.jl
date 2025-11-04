using DelimitedFiles
using DataFrames

function C_Rn_max_single_cell(res)
    maximum(maximum.(res.C_Rn_sup))
end

function C_Rn_sum_surface(res)
    maximum(sum.(res.C_Rn_sup))
end

function surface_concentrations(res, number_of_layers_from_surface::Int64; x::Int64=61, z::Int64=31)
    df_conc = DataFrame(years=parse.(Float64, collect(keys(res))), Concentrations=collect(values(res)))
    df_conc = sort(df_conc, :years)

    df_superficial_conc = map(row -> (row.years, reshape(row.Concentrations.C_Rn, (z, x))[z-number_of_layers_from_surface:end, :]), eachrow(df_conc))

    return DataFrame(df_superficial_conc, [:years, :C_Rn_sup])
end


function concentrations(output_file::String)
    regexs = Dict(
        "variable_regex" => r"(?<=\")[^,]*?(?=\")",
        "day_regex" => r"\d*\.\d{2,5}",
        "x_regex" => r"(?<=i=).*?(?=[,j=])",
        "z_regex" => r"(?<=j=).*?(?=,)",
    )
    extractors = Extractor(
        base -> begin
            file = joinpath(base, output_file)
            result, var, x, z = _concentrationplt2dict(
                file,
                regexs["variable_regex"],
                regexs["day_regex"],
                regexs["x_regex"],
                regexs["z_regex"]
            )
            return result
        end,
        Symbol("Concentrations"),
    )
    return extractors
end

function _concentrationplt2dict(file_path::String, var_regex::Regex, year_regex::Regex, x_regex::Regex, z_regex::Regex)
    data = readdlm(file_path, '\t', skipstart=0)
    result = Dict{Any,DataFrame}()
    parsed_dict = Dict()
    year = []
    x = []
    z = []
    vars = []
    t = 0
    for i in 1:size(data)[1]
        if occursin("variables", string(data[i]))
            push!(vars, [match.match for match in eachmatch(var_regex, data[1])])
            continue
        elseif occursin("text", string(data[i]))
            year = match(year_regex, data[i]).match
            parsed_dict[year] = []
            continue
        elseif occursin("zone", string(data[i]))
            push!(x, parse(Float64, match(x_regex, data[i]).match))
            push!(z, parse(Float64, match(z_regex, data[i]).match))
            continue
        else
            if typeof(data[i]) == Float64
                push!(parsed_dict[year], data[i])
            elseif typeof(data[i]) == SubString{String}
                for val in split(data[i])
                    push!(parsed_dict[year], parse(Float64, val))
                end
            else
                print("ERROR")
            end
        end
    end
    x = unique(x)[1]
    z = unique(z)[1]
    ### Checking matrix dimension
    for (key, val) in parsed_dict
        if length(parsed_dict[key]) != 3 * x * z && length(parsed_dict[key]) != 5 * x * z
            println("error parsing values")
        end
    end
    t = 0
    ## Model Matrix x-z
    x_values = []
    z_values = []
    for (key, val) in parsed_dict
        if length(parsed_dict[key]) == 5 * x * z && isempty(x_values) && isempty(z_values)
            append!(x_values, parsed_dict[key][1:5:end])
            append!(z_values, parsed_dict[key][2:5:end])
        end
    end
    ## Result dictionary
    for (key, val) in parsed_dict
        if length(parsed_dict[key]) == 5 * x * z
            result["$key"] = DataFrame()
            result["$key"][!, "C_salt"] = parsed_dict[key][3:5:end]
            result["$key"][!, "C_Rn"] = parsed_dict[key][4:5:end]
            result["$key"][!, "Spn"] = parsed_dict[key][5:5:end]
            result["$key"][!, :x] = x_values
            result["$key"][!, :z] = z_values
        elseif length(parsed_dict[key]) == 3 * x * z
            result["$key"] = DataFrame()
            result["$key"][!, "C_salt"] = parsed_dict[key][1:3:end]
            result["$key"][!, "C_Rn"] = parsed_dict[key][2:3:end]
            result["$key"][!, "Spn"] = parsed_dict[key][3:3:end]

            result["$key"][!, :x] = x_values
            result["$key"][!, :z] = z_values
        end
    end
    return result, vars, Int(x), Int(z)
end

function _sum_surface_radionuclide(input::Dict{Any,DataFrame}, number_of_layers::Int64)
    maximum_concentratation = []
    for v in values(input)
        matrix = reshape(v.C_Rn, (length(unique(v.z)), length(unique(v.x))))
        z_cell_value = length(unique(v.z)) - number_of_layers + 1
        surface_matrix = matrix[z_cell_value:length(unique(v.z)), :]
        push!(maximum_concentratation, sum(surface_matrix))
    end
    return maximum(maximum_concentratation)
end

function concentrations_surface_radionuclide(output_file::String, number_of_layers::Int64)
    regexs = Dict(
        "variable_regex" => r"(?<=\")[^,]*?(?=\")",
        "day_regex" => r"\d*\.\d{2,5}",
        "x_regex" => r"(?<=i=).*?(?=[,j=])",
        "z_regex" => r"(?<=j=).*?(?=,)",
    )
    extractors = Extractor(
        base -> begin
            file = joinpath(base, output_file)
            result, var, x, z = _concentrationplt2dict(
                file,
                regexs["variable_regex"],
                regexs["day_regex"],
                regexs["x_regex"],
                regexs["z_regex"]
            )
            return _sum_surface_radionuclide(result, number_of_layers)
        end,
        Symbol("C_Rn_surface"),
    )
    return extractors
end
