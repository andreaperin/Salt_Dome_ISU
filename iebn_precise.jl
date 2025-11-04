using EnhancedBayesianNetworks
# import EnhancedBayesianNetworks: evaluate!
using Cairo
using PGFPlotsX
using JLD2
using Dates

include("HC_Model/extractors/heads_model.jl")
include("HC_Model/extractors/extractors.jl")

const sim = MonteCarlo(10)
const cleanup = true
const number_of_layers = 20

earthquake_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:EarthQuake)
earthquake_cpt[:EarthQuake=>:yes_eq] = 10e-4
earthquake_cpt[:EarthQuake=>:no_eq] = 1 - 10e-4
earthquake_node = DiscreteNode(:EarthQuake, earthquake_cpt)

climate_change_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:ClimCh)
climate_change_cpt[:ClimCh=>:warmer] = 0.7
climate_change_cpt[:ClimCh=>:colder] = 0.2
climate_change_cpt[:ClimCh=>:stable] = 0.1
climate_change_node = DiscreteNode(:ClimCh, climate_change_cpt)

extreme_rain_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:Rain)
extreme_rain_cpt[:Rain=>:extreme_rain] = 0.4
extreme_rain_cpt[:Rain=>:normal_rain] = 0.6
extreme_rain_parameters = Dict(
    :extreme_rain => [Parameter(1.2, :head_factor)],
    :normal_rain => [Parameter(0.8, :head_factor)]
)
extreme_rain_node = DiscreteNode(:Rain, extreme_rain_cpt, extreme_rain_parameters)

time_scenario_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:Time)
time_scenario_cpt[:Time=>:short] = 0.5
time_scenario_cpt[:Time=>:long] = 0.5
time_scenario_parameters = Dict(
    :short => [Parameter(700.0, :sim_duration)],
    :long => [Parameter(7_000.0, :sim_duration)]
)
time_scenario_node = DiscreteNode(:Time, time_scenario_cpt, time_scenario_parameters)

long_disp_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
long_disp_cpt[] = Uniform(5, 40)
long_disp_node = ContinuousNode(:disp_lng, long_disp_cpt)

transv_disp_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
transv_disp_cpt[] = Uniform(0.5, 4)
transv_disp_node = ContinuousNode(:disp_trnsv, transv_disp_cpt)

hydro_cond_x_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
hydro_cond_x_cpt[] = truncated(Normal(9.81 * 10^(-6), 10^(-4)); lower=0, upper=0.001)
hydro_cond_x_node = ContinuousNode(:Kx, hydro_cond_x_cpt)

hydro_cond_z_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
hydro_cond_z_cpt[] = truncated(Normal(9.81 * 10^(-6), 10^(-4)); lower=0, upper=0.001)
hydro_cond_z_node = ContinuousNode(:Kz, hydro_cond_x_cpt)

porosity_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EarthQuake)
porosity_cpt[:EarthQuake=>:no_eq] = truncated(Normal(0.1, 0.05); lower=0, upper=4)
porosity_cpt[:EarthQuake=>:yes_eq] = truncated(Normal(0.3, 0.05); lower=0, upper=4)
porosity_node = ContinuousNode(:porosity, porosity_cpt)

diffusion_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ClimCh)
diffusion_cpt[:ClimCh=>:warmer] = truncated(Normal(2 * 10^(-7), 10^(-6)); lower=0, upper=0.00001)
diffusion_cpt[:ClimCh=>:colder] = truncated(Normal(2 * 10^(-8), 10^(-7)); lower=0, upper=0.00001)
diffusion_cpt[:ClimCh=>:stable] = truncated(Normal(2 * 10^(-9), 10^(-6)); lower=0, upper=0.00001)
diffusion_node = ContinuousNode(:diff_coeff, diffusion_cpt)


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
performance = df -> -df.C_Rn_surface

HC_node = DiscreteFunctionalNode(:HC, models, performance, sim)

nodes = [earthquake_node, climate_change_node, extreme_rain_node, time_scenario_node, long_disp_node, transv_disp_node, hydro_cond_x_node, hydro_cond_z_node, porosity_node, diffusion_node, HC_node]

ebn = EnhancedBayesianNetwork(nodes)
add_child!(ebn, :EarthQuake, :porosity)
add_child!(ebn, :ClimCh, :diff_coeff)
add_child!(ebn, :diff_coeff, :HC)
add_child!(ebn, :porosity, :HC)
add_child!(ebn, :Rain, :HC)
add_child!(ebn, :Time, :HC)
add_child!(ebn, :disp_lng, :HC)
add_child!(ebn, :disp_trnsv, :HC)
add_child!(ebn, :Kx, :HC)
add_child!(ebn, :Kz, :HC)
order!(ebn)

fig_path = joinpath(cwd, "imgs")
ebn_fig = gplot(ebn, NODESIZEFACTOR=0.17, NODELABELSIZE=3.6)
# draw(PDF(joinpath(fig_path, "1_ebn_precise.pdf"), 16cm, 16cm), ebn_fig)

EnhancedBayesianNetworks.evaluate!(ebn, true, false)

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
path = joinpath(pwd(), "results/iebn")
name = timestamp * "_" * string(sim) * ".jld2"
@save joinpath(path, name) ebn