using ALIFE2023
using ALIFE2023.Inverter
using ModelingToolkit
using Plots
using Graphs
using StatsBase
using StatsPlots
using GraphRecipes
using CSV
using Tables

function build_inverter_experiment(M, N)
    N = length(Inverter.components(N))
    n = 2 * ((N - 2)^2 + N - 1)
    MR = 4 / n
    return Algorithm(
        random_population(M, N, n, 1, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 1, 1),
        decoder(N, 1, 1),
        Inverter.selector(N),
    )
end

function evolve_inverter(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        @show g, alg.scores
        scores[g, :] .= alg.scores
    end
    return scores
end

@time bscore = Inverter.benchmark_score(2048)
@show bscore

alg = build_inverter_experiment(1, 7)
D = decoder(length(Inverter.components()), 1, 1)

function cgp_convergence_iterations(N)
    i = 0
    alg = build_inverter_experiment(1, N)
    bscore = Inverter.benchmark_score(128)
    @show bscore
    while true
        i = i + 1
        evolve!(alg, 4)
        @show i, alg.scores
        g = alg.decode(alg.population[1, :]) 
        if g.ne > 0
            score = Inverter.score(g, 128)
            @show score
            score >= bscore && return i
        end
    end
    return i
end

function random_graph_convergence_iterations(N)
    i = 0
    N = length(Inverter.components(N))
    alg = build_inverter_experiment(1, N)
    n = 2 * ((N - 2)^2 + N - 1)
    while true
        i = i + 1
        g = random_genome(N, n, 1, 1)
        if alg.decode(g).ne > 0
            score, swap = alg.select(alg.decode(g), Inverter.benchmark_graph(N))
            swap && return i
        end
    end
end

function average_random_graph_convergence(n, N)
    is = zeros(Int, n)
    for i in eachindex(is)
        @show i
        is[i] = random_graph_convergence_iterations(N)
    end
    return is
end
