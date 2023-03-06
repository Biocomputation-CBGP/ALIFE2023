using ALIFE2023
using ALIFE2023.Bandpass
using DifferentialEquations
using Plots
using Graphs
using StatsBase
using GraphRecipes


function build_bandpass_experiment(M)
    N = length(Bandpass.components())
    n = 2 * ((N - 2)^2 + N - 1)
    MR = 3 / n
    return Algorithm(
        random_population(M, N, n, 2, 1),
        mutation_operator(N, MR, 2, 1),
        decoder(N),
        Bandpass.selector(),
    )
end

function evolve_bandpass(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        x = ALIFE2023.scores(alg, Bandpass.score)
        @show g, x
        scores[g, :] .= x
    end
    return scores
end

alg = build_bandpass_experiment(1)
