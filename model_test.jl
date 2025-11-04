using UncertaintyQuantification

include("HC_Model/extractors/heads_model.jl")
include("HC_Model/extractors/extractors.jl")

const number_of_layers = 10
const cleanup = true

Kx = RandomVariable(truncated(Normal(9.81 * 10^(-6), 10^(-4)); lower=0, upper=0.001), :Kx)
Kz = RandomVariable(truncated(Normal(9.81 * 10^(-6), 10^(-4)); lower=0, upper=0.001), :Kz)
disp_long_h = RandomVariable(Uniform(10, 60), :disp_lng)
disp_long_v = RandomVariable(Uniform(1, 6), :disp_trnsv)
diff_coefficient = RandomVariable(truncated(Normal(2(10^(-8)), 10^(-7)); lower=0, upper=0.00001), :diff_coeff)
porosity = RandomVariable(truncated(Normal(0.15, 0.1); lower=0, upper=0.6), :porosity)
sim_duration = Parameter(1_000, :sim_duration)
head_factor = RandomVariable(Uniform(0.8, 1.2), :head_factor)

inputs = [head_factor, Kz, Kx, disp_long_h, disp_long_v, porosity, diff_coefficient, sim_duration]

output_file = "smoker_cxz.plt"
cwd = pwd()
path = joinpath(cwd, "HC_Model/smoker_compiled/smokerV3TC")
solver = Solver(path, "", "")
sourcedir = joinpath(cwd, "HC_Model/smoker_compiled")
sources = ["smoker.data"]
workdir = joinpath(cwd, "runs/temp")

extractors = [concentrations(output_file)]
model_HC = ExternalModel(sourcedir, sources, extractors, solver; workdir=workdir, cleanup=cleanup)
model_superficial_concentration = Model(df -> surface_concentrations.(df.Concentrations, number_of_layers), :C_Rn_surface)
model_max_single_cell_over_time = Model(df -> C_Rn_max_single_cell.(df.C_Rn_surface), :C_Max_single_cell)
model_sum_surface_over_time = Model(df -> C_Rn_sum_surface.(df.C_Rn_surface), :C_Rn_cumulative)
models = [heads_model, model_HC, model_superficial_concentration, model_max_single_cell_over_time, model_sum_surface_over_time]

samples = sample(inputs, 2)
@time evaluate!(models, samples)
