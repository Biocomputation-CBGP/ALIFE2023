using ALIFE2023
using ECCOCGP
using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
using Graphs
using StatsBase
using StatsPlots
using GraphRecipes
using CSV
using Tables
using Distributions

function SomeOperon(;name)
    R = Dimer(1, 1, 1/60, name=name)
    return RegulatedPromoter(
        0.01, 0.5, R, 0.1, 1.0;
        name=Symbol("p", name)
    )
end

function components()
    @named YFP = Monomer(1)
    @named pYFP = RegulatedPromoter(0, 0, YFP, 0, 0)
    @named SrpR = SomeOperon()
    @named PhiF = SomeOperon()
    @named AmeR = SomeOperon()
    @named QacR = SomeOperon()
    @named BetI = SomeOperon()
    @named Lmra = SomeOperon()
    return pSrpR, pPhiF, pAmeR, pQacR, pBetI, pLmra, pYFP
end

function bin(x::Vector{T}, upper::Int) where {T<:Integer}
    binned = zeros(T, upper + 1)
    for i in 0:upper
        binned[i+1] = count(isequal(i), x)
    end
    return binned
end

function bin(x::Vector{T}, upper::Int) where {T<:AbstractFloat}
    return bin(Int.(floor.(x)), upper)
end

function earthmovers(x::Vector{T1}, y::Vector{T2}) where {T1<:Real, T2<:Real}
    upper = max(maximum(floor.(x)), maximum(floor.(y)))
    P = bin(x, Int(upper))
    Q = bin(y, Int(upper))
    return sum(abs.(cumsum(P .- Q)))
end

function score(graph, systems)
    model = model_from_graph(graph, systems)
    problem = problem_from_model(model)
    N = 1024
    samples = sample_from_problem(problem, model, YFP.monomer, N)
    targets = rand(Binomial(100, 0.5), N)
    return -earthmovers(targets, samples.u)
end

function scorer(systems)
    let systems=systems
        function (G::SimpleDiGraph{Int})
            G.ne == 0 && return -Inf
            return score(G, systems)
        end
    end
end

function build_dist_experiment(systems, M, n, MR)
    N = length(systems)
    return CGP(
        random_population(M, N, n, 0, 1),
        fill(-Inf, M),
        mutation_operator(N, MR, 0, 1),
        decoder(N, 0, 1),
        scorer(systems),
    )
end

function evolve_dist(alg, G)
    for g in 1:G
        ECCOCGP.step!(alg, 4)
        @show g, alg.scores
    end
end

