using ALIFE2023
using ALIFE2023.And
using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
using Graphs
using StatsBase
using GraphRecipes
using Random
using CSV
using Tables


function build_and_experiment(M)
    N = length(And.components())
    n = 2 * ((N - 2)^2 + N - 1)
    MR = 4 / n
    return Algorithm(
        random_population(M, N, n, 2, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 2, 1),
        decoder(N, 2, 1),
        And.selector(),
    )
end

function evolve_and(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        scores[g, :] .= alg.scores
        @show g, scores[g, :]
        _, i = findmax(alg.scores)
        model = Circuit(alg.decode(alg.population[i, :]), And.components(); name=:model)
        output = (@nonamespace model.YFP).monomer
        problem = problem_from_model(model)
        a, b, c, d = And.distributions(problem, model, 256)
        plt = StatsPlots.density(a.u, label="00")
        StatsPlots.density!(b.u, label="01")
        StatsPlots.density!(c.u, label="10")
        StatsPlots.density!(d.u, label="11")
        display(plt)
    end
    return scores
end

@show And.benchmark_score(256)
alg = build_and_experiment(1)
