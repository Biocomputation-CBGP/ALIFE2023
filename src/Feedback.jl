module Feedback

using GeneticLogicGraph
using Graphs
using ModelingToolkit
using Catalyst
using StatsBase
using SciMLBase
using JumpProcesses
using ..ALIFE2023

function SomeOperon(;name)
    R = Dimer(1, 1, 1/60, name=name)
    return RegulatedPromoter(1/50, 1, R, 1/10, 1; name=Symbol("p", name))
end

function components()
    @named YFP = Monomer(1/10)
    @named pYFP = RegulatedPromoter(0, 0, YFP, 0, 0)
    @named SrpR = SomeOperon()
    @named PhiF = SomeOperon()
    @named AmeR = SomeOperon()
    @named QacR = SomeOperon()
    @named BetI = SomeOperon()
    @named Lmra = SomeOperon()
    return SrpR, PhiF, AmeR, QacR, BetI, Lmra, pYFP
end

function benchmark_graph(i::Int=1)
    N = length(components())
    G = SimpleDiGraph(N)
    add_edge!(G, i, N)
    add_edge!(G, i, i)
    return G
end

function super_benchmark_graph()
    N = length(components())
    G = SimpleDiGraph(N)
    for i in 1:N-1
        add_edge!(G, i, N)
        add_edge!(G, i, i)
    end
    return G
end

benchmark_model() = Circuit(benchmark_graph(), components(); name=:x)
super_benchmark_model() = Circuit(super_benchmark_graph(), components(); name=:x)

benchmark_problem() = problem_from_model(benchmark_model())
super_benchmark_problem() = problem_from_model(super_benchmark_model())

function ordinary_model()
    @variables t
    @species N(t)=0
    @parameters α μ
    @named pJS2321 = PromoterRegion(0.5)
    @named YFP = Monomer(1)
    rxs = [
        transcription(pJS2321, YFP);
        translation(YFP);
        mrna_degradation(YFP, α);
        [Reaction(μ, nothing, [N])];
    ]
    return ReactionSystem(
        rxs, t, [N], [μ, α];
        name=:x, systems=[pJS2321, YFP], connection_type=(Circuit,)
    )
end

function super_ordinary_model()
    @variables t
    @species N(t)=0
    @parameters α μ
    promoters = [PromoterRegion(0.5; name=Symbol(:p, i)) for i in 1:6]
    @named YFP = Monomer(1)
    rxs = [
        reduce(vcat, transcription(p, YFP) for p in promoters);
        translation(YFP);
        mrna_degradation(YFP, α);
        [Reaction(μ, nothing, [N])];
    ]
    return ReactionSystem(
        rxs, t, [N], [μ, α];
        name=:x, systems=[promoters; YFP], connection_type=(Circuit,)
    )
end

ordinary_problem() = problem_from_model(ordinary_model())
super_ordinary_problem() = problem_from_model(super_ordinary_model())

function ordinary_distribution(N)
    return sample_from_problem(
        ordinary_problem(),
        ordinary_model(),
        (@nonamespace ordinary_model().YFP).monomer,
        N,
    )
end

function super_ordinary_distribution(N)
    return sample_from_problem(
        super_ordinary_problem(),
        super_ordinary_model(),
        (@nonamespace super_ordinary_model().YFP).monomer,
        N,
    )
end

function benchmark_distribution(N)
    return sample_from_problem(
        benchmark_problem(),
        benchmark_model(),
        (@nonamespace benchmark_model().YFP).monomer,
        N,
    )
end

function super_benchmark_distribution(N)
    return sample_from_problem(
        super_benchmark_problem(),
        super_benchmark_model(),
        (@nonamespace super_benchmark_model().YFP).monomer,
        N,
    )
end

function ordinary_score(N)
    distribution = ordinary_distribution(N)
    return mean(distribution) / std(distribution)
end

function super_ordinary_score(N)
    distribution = super_ordinary_distribution(N)
    return mean(distribution) / std(distribution)
end

function benchmark_score(N)
    distribution = benchmark_distribution(N)
    return mean(distribution) / std(distribution)
end

function super_benchmark_score(N)
    distribution = super_benchmark_distribution(N)
    return mean(distribution) / std(distribution)
end

function score(G::DiGraph, N::Int=256)
    @named model = Circuit(G, components())
    output = (@nonamespace model.YFP).monomer
    problem = problem_from_model(model)
    samples = sample_from_problem(problem, model, output, N)
    return mean(samples) / std(samples)
end

function selector()
    let systems=components()
        function (G1::T, G2::T) where {T<:DiGraph{Int}}
            G1.ne == 0 && return G2
            G2.ne == 0 && return G1

            @named M1 = Circuit(G1, systems)
            @named M2 = Circuit(G2, systems)
            P1 = remake(problem_from_model(M1), tspan=(0, 1440 + rand() * 180))
            P2 = remake(problem_from_model(M2), tspan=(0, 1440 + rand() * 180))
            i = findfirst(isequal((@nonamespace M1.YFP).monomer), states(M1))
            j = findfirst(isequal((@nonamespace M2.YFP).monomer), states(M2))

            f = let prob=P1, i=i
                () -> solve(prob, SSAStepper())[end][i]
            end
            g = let prob=P2, i=j
                () -> solve(prob, SSAStepper())[end][i]
            end
            
            return snr_incremental(f, g)
        end
    end
end

end
