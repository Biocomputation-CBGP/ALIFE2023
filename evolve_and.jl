using ALIFE2023
using ALIFE2023.And
using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
using Graphs
using StatsBase
using StatsPlots
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

function plot_distributions(G)
    model = Circuit(G, And.components(); name=:model)
    problem = problem_from_model(model)
    a, b, c, d = And.distributions(problem, model, 256)
    plt = StatsPlots.density(a.u, label="00")
    StatsPlots.density!(b.u, label="01")
    StatsPlots.density!(c.u, label="10")
    StatsPlots.density!(d.u, label="11")
    plot!(plt, legend=:best, xlabel="Abundance of output", ylabel="Frequency density")
    plot!(plt, tickfontsize=8, legendfontsize=8, guidefontsize=8)
    plot!(plt, size=(350, 275))
    return plt
end

function evolve_and(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        @time evolve!(alg, 4)
        scores[g, :] .= alg.scores
        @show g, scores[g, :]
        if g % 8 == 1
            _, i = findmax(alg.scores)
            plt = plot_distributions(alg.decode(alg.population[i, :]))
            savefig(plt, "evolve-and-$g-dists.svg")
            display(plt)
        end
    end
    return scores
end

@show And.benchmark_score(256)
# a, b, c, d = And.distributions(And.benchmark_problem(), And.benchmark_model(), 256)

# Plots.theme(:dao)
# plt = plot_distributions(And.benchmark_graph())
# display(plt)

function runner(G)
    alg = build_and_experiment(1)
    scores = evolve_and(alg, G)
    CSV.write("and-population.csv", Tables.table(alg.population), writeheader=false)
    CSV.write("and-scores.csv", Tables.table(scores), writeheader=false)
end




