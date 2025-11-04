using UncertaintyQuantification
using JLD2
using DataFrames
# using Plots
include("HC_Model/extractors/heads_model.jl")
include("HC_Model/extractors/extractors.jl")

const number_of_layers = 20
const sobols_simulation = 2
const cleanup = true

Kx = RandomVariable(Uniform(10e-10, 0.001), :Kx)
Kz = RandomVariable(Uniform(10e-10, 0.001), :Kz)
disp_long_h = RandomVariable(Uniform(5, 40), :disp_lng)
disp_long_v = RandomVariable(Uniform(0.5, 4), :disp_trnsv)
diff_coefficient = RandomVariable(Uniform(10e-11, 0.000001), :diff_coeff)
porosity = RandomVariable(Uniform(10e-3, 0.6), :porosity)
sim_duration = Parameter(10_000, :sim_duration)
head_factor = RandomVariable(Uniform(0.8, 1.2), :head_factor)

inputs = [head_factor, Kz, Kx, disp_long_h, disp_long_v, porosity, diff_coefficient, sim_duration]

output_file = "smoker_cxz.plt"
cwd = pwd()
path = joinpath(cwd, "HC_Model/smoker_compiled/smokerV3TC")
solver = Solver(path, "", "")
sourcedir = joinpath(cwd, "HC_Model/smoker_compiled/")
sources = ["smoker.data"]
workdir = joinpath(cwd, "runs/temp")
extractors = [concentrations(output_file), concentrations_surface_radionuclide(output_file, number_of_layers)]

model_HC = ExternalModel(sourcedir, sources, extractors, solver; workdir=workdir, cleanup=cleanup)
models = [heads_model, model_HC]

# df = sample(inputs, 2)
# evaluate!(models, df)

sis = sobolindices(models, inputs, [:C_Max_single_cell, :C_Rn_cumulative], MonteCarlo(sobols_simulation))

