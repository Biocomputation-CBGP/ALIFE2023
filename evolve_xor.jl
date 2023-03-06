using ALIFE2023
using ALIFE2023.Xor
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


function build_xor_experiment(M)
    N = length(Xor.components())
    n = 2 * ((N - 2)^2 + N - 1)
    MR = 3 / n
    return Algorithm(
        random_population(M, N, n, 2, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 2, 1),
        decoder(N, 2, 1),
        Xor.selector(),
    )
end

function evolve_xor(alg, G)
    scores = Matrix{Float64}(undef, G, size(alg, 1))
    for g in 1:G
        evolve!(alg, 4)
        @show g, alg.scores
        scores[g, :] .= alg.scores
        _, i = findmax(alg.scores)
        model = JumpCircuit(alg.decode(alg.population[i, :]), Xor.components(); name=:model)
        output = (@nonamespace model.YFP).monomer
        problem = problem_from_model(model)
        a, b, c, d = Xor.distributions(problem, model, 256)
        plt = StatsPlots.density(a.u, label="00")
        StatsPlots.density!(b.u, label="01")
        StatsPlots.density!(c.u, label="10")
        StatsPlots.density!(d.u, label="11")
        display(plt)
    end
    return scores
end

@show Xor.benchmark_score(2048)

alg = build_xor_experiment(1)

# evolve_xor(alg, 32)

# save_alg(alg, "xor_cgp_1_32_3_256")



