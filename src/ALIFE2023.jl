module ALIFE2023

using GeneticLogicGraph
using ModelingToolkit
using JumpProcesses
using Graphs
using StatsBase
using Catalyst
using Plots
using GraphRecipes
using ClusterManagers


include("problems.jl")
export model_from_graph
export problem_from_model
export change_input_levels
export sample_from_problem
export snr_incremental


include("algorithm.jl")
export random_genome
export random_population
export mutation_operator
export decoder
export Algorithm
export evolve!


include("Feedback.jl")
include("Inverter.jl")
include("Xor.jl")
include("And.jl")
include("Bandpass.jl")
include("Amplify.jl")

end # module ALIFE2023
