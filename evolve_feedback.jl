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
    MR = 4 / n
    
    return Algorithm(
        random_population(M, N, n, 0, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 0, 1),
        decoder(N, 0, 1),
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

@time f1 = Feedback.ordinary_score(2048)
@time f2 = Feedback.super_ordinary_score(2048)
@time f3 = Feedback.benchmark_score(2048)
@time f4 = Feedback.super_benchmark_score(2048)
@show f1, f2, f3, f4

function runner(G)
    alg = build_feedback_experiment(1)
    scores = evolve_feedback(alg, G)
    CSV.write("feedback-population.csv", Tables.table(alg.population), writeheader=false)
    CSV.write("feedback-scores.csv", Tables.table(scores), writeheader=false)
end
