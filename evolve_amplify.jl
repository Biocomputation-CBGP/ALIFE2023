using ALIFE2023
using ALIFE2023.Amplify
using ModelingToolkit
using Plots
using Graphs
using StatsBase
using StatsPlots
using GraphRecipes
using CSV
using Tables

function build_amplify_experiment(M)
    N = length(Amplify.components())
    n = 2 * ((N - 2)^2 + N - 1)
    MR = 3 / n
    return Algorithm(
        random_population(M, N, n, 1, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 1, 1),
        decoder(N, 1, 1),
        Amplify.selector(),
    )
end

function evolve_amplify(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        @show g, alg.scores
        scores[g, :] .= alg.scores
    end
    return scores
end

@show Amplify.benchmark_score(512)
alg = build_amplify_experiment(1)
