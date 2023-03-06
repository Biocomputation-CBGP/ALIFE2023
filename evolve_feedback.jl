using ALIFE2023
using ALIFE2023.Feedback
using ModelingToolkit
using Plots
using Graphs
using StatsBase
using StatsPlots
using GraphRecipes

function build_feedback_experiment(M)
    N = length(Feedback.components())
    n = 2 * ((N - 1)^2) + N - 1
    MR = 3 / n
    
    return Algorithm(
        random_population(M, N, n, 0, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 0, 1),
        decoder(N),
        Feedback.selector(),
    )
end

function evolve_feedback(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        scores[g, :] .= alg.scores
        @show g, scores[g, :]
    end
    return scores
end

function random_iterations(N)
    i = 1
    D = decoder(N)
    n = 2 * (N - 1)^2 + (N - 1)
     while Feedback.score(D(random_genome(N, n, 0, 1)), 512) < 2
        i = i + 1
    end
    return i
end

function cgp_convergence_iterations(N)
    i = 0
    alg = build_feedback_experiment(components()[end-N+1:end], 4)
    while score(alg.decode(alg.population[1,:]), 512) < 2
        i = i + 1
        evolve!(alg, 4)
    end
    return i
end

@time f1 = Feedback.ordinary_score(2048)
@time f2 = Feedback.super_ordinary_score(2048)
@time f3 = Feedback.benchmark_score(2048)
@time f4 = Feedback.super_benchmark_score(2048)

alg = build_feedback_experiment(1)
