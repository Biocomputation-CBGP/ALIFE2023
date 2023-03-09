module Inverter

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

function components(N::Int=7)
    @named YFP = Monomer(1)
    @named pYFP = RegulatedPromoter(0, 0, YFP, 0, 0)
    @named LacI = InputSpecies(1 * log(2) / 90)
    @named pTac = RegulatedPromoter(1/50, 1, LacI, 1/10, 1)
    @named SrpR = SomeOperon()
    @named PhiF = SomeOperon()
    @named AmeR = SomeOperon()
    @named QacR = SomeOperon()
    @named BetI = SomeOperon()
    @named AmtR = SomeOperon()
    @named Lmra = SomeOperon()
    return Tuple(vcat([pTac], [SrpR, PhiF, AmeR, QacR, BetI, AmtR, Lmra][1:N-2], [pYFP]))
end

function benchmark_graph(N::Int=7)
    N = length(components(N))
    G = SimpleDiGraph(N)
    add_edge!(G, 1, N)
    return G
end

benchmark_model() = Circuit(benchmark_graph(), components(); name=:x)
benchmark_problem() = problem_from_model(benchmark_model())
benchmark_distribution(N) = complement(benchmark_problem(), benchmark_model(), N)

function benchmark_score(N)
    distribution = benchmark_distribution(N)
    return mean(distribution) / std(distribution)
end

function complement(problem, model, N)
    input = (@nonamespace model.LacI).λ
    output = (@nonamespace model.YFP).monomer

    problem = change_input_levels(problem, [input], [10 * log(2) / 90])
    los = sample_from_problem(problem, model, output, N).u
    
    problem = change_input_levels(problem, [input], [100 * log(2) / 90])
    his = sample_from_problem(problem, model, output, N).u
    
    return los .- his
end

function score(G::DiGraph, n::Int=256)
    N = first(size(G))
    model = Circuit(G, components(N); name=:x)
    problem = problem_from_model(model)
    distribution = complement(problem, model, n)
    return mean(distribution) / std(distribution)
end

function difference_sampler(problem, inputs, levels, output)
    idx = findfirst(isequal(output), states(problem.prob.f.sys))
    let problem=problem, inputs=inputs, levels=levels, idx=idx
        function ()
            idx == nothing && return 0.0
            problem = change_input_levels(problem, inputs, levels[1] .* log(2) ./ 90)
            T = 1440 + rand() * 180
            problem = remake(problem, tspan=(0, T))
            hi = solve(problem, SSAStepper())[end][idx]
            problem = change_input_levels(problem, inputs, levels[2] .* log(2) ./ 90)
            lo = solve(problem, SSAStepper())[end][idx]
            return hi - lo
        end
    end
end

function selector(N)
    let systems=components(N)
        function (G::T, H::T) where {T<:DiGraph}
            x = Circuit(G, systems; name=:x)
            y = Circuit(H, systems; name=:y)
            input  = (@nonamespace x.LacI).λ
            output = (@nonamespace x.YFP).monomer
            
            x = problem_from_model(x)
            y = problem_from_model(y)

            f = difference_sampler(x, [input], [[10], [100]], output)
            g = difference_sampler(y, [input], [[10], [100]], output)

            return snr_incremental(f, g)
        end
    end
end

function scorer(N)
    let systems=components(N)
        function (G::T) where {T<:DiGraph}
            x = Circuit(G, systems; name=:x)
            input  = (@nonamespace x.LacI).λ
            output = (@nonamespace x.YFP).monomer
            
            x = problem_from_model(x)

            f = difference_sampler(x, [input], [[10], [100]], output)
            return f()
        end
    end
end

end
