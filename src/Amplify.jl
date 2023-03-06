module Amplify

using GeneticLogicGraph
using Graphs
using ModelingToolkit
using Catalyst
using Random
using StatsBase
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
    @named pSrpR = SomeOperon()
    @named pPhiF = SomeOperon()
    @named pAmeR = SomeOperon()
    @named pQacR = SomeOperon()
    @named pBetI = SomeOperon()
    return pTac, pSrpR, pPhiF, pAmeR, pQacR, pBetI, pYFP
end

function complement(problem, model, N)
    input = (@nonamespace model.LacI).λ
    output = (@nonamespace model.YFP).monomer

    problem = change_input_levels(problem, [input], [100 * log(2) / 90])
    his = sample_from_problem(problem, model, output, N)
    
    problem = change_input_levels(problem, [input], [10 * log(2) / 90])
    los = sample_from_problem(problem, model, output, N)
    
    return his .- los
end

function score(G::DiGraph, N::Int=256)
    model = model_from_graph(G, components())
    output = (@nonamespace model.YFP).monomer
    problem = problem_from_model(model)
    distribution = complement(problem, model, N)
    return mean(distribution) / std(distribution)
end

function difference_sampler(problem, inputs, levels, output)
    idx = findfirst(isequal(output), states(problem.prob.f.sys))
    let problem=problem, inputs=inputs, levels=levels, idx=idx
        function ()
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

function selector()
    let systems=components()
        function (G::T, H::T) where {T<:DiGraph}
            G.ne == 0 && return H
            H.ne == 0 && return G
            
            x = model_from_graph(G, systems)
            y = model_from_graph(H, systems)
            input  = (@nonamespace x.LacI).λ
            output = (@nonamespace x.YFP).monomer
            
            x = problem_from_model(x)
            y = problem_from_model(y)

            f = difference_sampler(x, [input], [[100], [10]], output)
            g = difference_sampler(y, [input], [[100], [10]], output)

            snr, xbigger = snr_incremental(f, g)
            return xbigger
        end
    end
end

end
