module And

using GeneticLogicGraph
using Graphs
using ModelingToolkit
using Catalyst
using Random
using SciMLBase
using JumpProcesses
using StatsBase
using StatsPlots
using ..ALIFE2023

function SomeOperon(;name)
    R = Dimer(1, 1, 1/60, name=name)
    return RegulatedPromoter(1/50, 1, R, 1/10, 1; name=Symbol("p", name))
end

function components()
    @named YFP = Monomer(1)
    @named pYFP = RegulatedPromoter(0, 0, YFP, 0, 0)
    @named LacI = InputSpecies(1 * log(2) / 90)
    @named pTac = RegulatedPromoter(1/50, 1, LacI, 1/10, 1)
    @named aTc  = InputSpecies(1 * log(2) / 90)
    @named pTet = RegulatedPromoter(1/50, 1, aTc, 1/10, 1)
    @named pSrpR = SomeOperon()
    @named pPhiF = SomeOperon()
    @named pAmeR = SomeOperon()
    @named pQacR = SomeOperon()
    @named pBetI = SomeOperon()
    @named pLmra = SomeOperon()
    @named pAmtR = SomeOperon()
    return pTac, pTet, pSrpR, pPhiF, pAmeR, pQacR, pBetI, pLmra, pYFP
end

function benchmark_graph()
    N = length(components())
    G = SimpleDiGraph(N)
    add_edge!(G, 1, 3)
    add_edge!(G, 2, 4)
    add_edge!(G, 3, 5)
    add_edge!(G, 4, 5)
    add_edge!(G, 5, 6)
    add_edge!(G, 6, N)
    return G
end

benchmark_model() = Circuit(benchmark_graph(), components(); name=:x)
benchmark_problem() = problem_from_model(benchmark_model())
benchmark_distribution(N) = complement(benchmark_problem(), benchmark_model(), N)

function benchmark_score(N)
    distribution = benchmark_distribution(N)
    return mean(distribution) / std(distribution)
end

function distributions(problem, model, N)
    inputs = [(@nonamespace model.LacI).λ, (@nonamespace model.aTc).λ]
    a = change_input_levels(problem, inputs, [1, 1] .* log(2) ./ 90)
    b = change_input_levels(problem, inputs, [1, 100] .* log(2) ./ 90)
    c = change_input_levels(problem, inputs, [100, 1] .* log(2) ./ 90)
    d = change_input_levels(problem, inputs, [100, 100] .* log(2) ./ 90)

    output = (@nonamespace model.YFP).monomer
    a = sample_from_problem(a, model, output, N)
    b = sample_from_problem(b, model, output, N)
    c = sample_from_problem(c, model, output, N)
    d = sample_from_problem(d, model, output, N)

    return a, b, c, d
end

function complement(problem, model, N)
    a, b, c, d = distributions(problem, model, N)
    supposed_to_be_lows = vcat(a.u, b.u, c.u)
    supposed_to_be_lows = sample(supposed_to_be_lows, N; replace=false)
    
    supposed_to_be_highs = d.u
    
    return supposed_to_be_highs .- supposed_to_be_lows
end

function score(G::DiGraph, N::Int=256)
    model = Circuit(G, components(); name=:x)
    output = (@nonamespace model.YFP).monomer
    problem = problem_from_model(model)
    distribution = complement(problem, model, N)
    return mean(distribution) / std(distribution)
end

function difference_sampler(problem, inputs, levels, output)
    idx = findfirst(isequal(output), states(problem.prob.f.sys))
    let problem=problem, inputs=inputs, levels=levels .* log(2) ./ 90, idx=idx
        function ()
            T = 1440 + rand() * 180
            problem = remake(problem, tspan=(0, T))

            lo_input = rand([1, 2, 3])
            hi_input = 4
            problem = change_input_levels(problem, inputs, levels[hi_input])
            hi = solve(problem, SSAStepper())[end][idx]
            
            problem = change_input_levels(problem, inputs, levels[lo_input])
            lo = solve(problem, SSAStepper())[end][idx]

            return hi - lo
        end
    end
end

function selector()
    let systems=components()
        function (G::T, H::T) where {T<:DiGraph}
            G.ne == 0 && return H
            H.ne == 0 && return G
            
            x = Circuit(G, systems; name=:x)
            y = Circuit(H, systems; name=:y)
            inputs = [(@nonamespace x.LacI).λ, (@nonamespace x.aTc).λ]
            output = (@nonamespace x.YFP).monomer
            
            x = problem_from_model(x)
            y = problem_from_model(y)
            levels = [[1, 1], [1, 100], [100, 1], [100, 100]]
            f = difference_sampler(x, inputs, levels, output)
            g = difference_sampler(y, inputs, levels, output)

            return snr_incremental(f, g)
        end
    end
end

end
