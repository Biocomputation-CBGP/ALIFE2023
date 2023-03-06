using ECCOCGP
using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
using Graphs
using StatsBase
using GraphRecipes

function SomeOperon(;name)
    R = Dimer(1, 0.25, 0.01, 1, 0.1, name=name)
    return RegulatedPromoter(
        0.1, 1, R, 0.1, 0.1;
        name=Symbol("p", name)
    )
end


const YFP = Monomer(1, 0.25, 0.01; name=:YFP)
const pYFP = RegulatedPromoter(0, 0, YFP, 0, 0; name=:pYFP)
const LacI = InputSpecies(1, 0.1; name=:LacI)
const pTac = RegulatedPromoter(0.01, 1, LacI, 0.1, 0.1; name=:pTac)
const pSrpR = SomeOperon(name=:SrpR)
const pPhiF = SomeOperon(name=:PhiF)
const pAmeR = SomeOperon(name=:AmeR)
const pQacR = SomeOperon(name=:QacR)
const pBetI = SomeOperon(name=:BetI)
const pLmra = SomeOperon(name=:Lmra)
const pAmtR = SomeOperon(name=:AmtR)

const components = (pTac, pSrpR, pPhiF, pAmeR, pQacR, pBetI, pLmra, pAmtR, pYFP);

function model_from_graph(graph, systems)
    return Circuit(graph, systems; name=gensym("circuit"))
end

function problem_from_model(model)
    u0 = reduce(merge, [randu0(sys) for sys in model.systems], init=Dict())
    problem = DiscreteProblem(model, u0, (0.0, 3000.0), Dict(:μ => 1/90))
    cb = CallbackSet(make_doubling_callback(model), make_termination_callback(1e-12))
    return JumpProblem(model, problem, Direct(), callback=cb, save_positions=(false, false))
end

function benchmark_graph(systems)
    G = SimpleDiGraph(length(systems))
    add_edge!(G, 1, 2)
    add_edge!(G, 2, length(systems))
    return G
end

function benchmark_problem(systems)
    G = benchmark_graph(systems)
    return problem_from_model(model_from_graph(G, systems))
end

function change_input_level(problem, input, level)
    parameter_idx = findfirst(isequal(input), parameters(problem.prob.f.sys))
    p = copy(problem.prob.p)
    p[parameter_idx] = level
    return remake(problem, p=p)
end

function sample_from_problem(problem, output, Ts::Vector{<:Real})::Vector{Int}
    output_idx = findfirst(isequal(output), states(problem.prob.f.sys))
    problem = remake(problem, tspan=(0.0, maximum(Ts)))
    solution = solve(problem, SSAStepper(), saveat=sort(Ts))
    return [solution(T)[output_idx] for T in Ts]
end

function complement(graph, systems, N)
    problem = problem_from_model(model_from_graph(graph, systems))
    problemlo = change_input_level(problem, LacI.λ, 1)
    problemhi = change_input_level(problem, LacI.λ, 10)
    Ts = 1500 .+ rand(N) .* 45 .* N

    los = sample_from_problem(problemlo, YFP.monomer, Ts)
    any(isequal(typemax(Int)), los) && return typemin(Int)

    his = sample_from_problem(problemhi, YFP.monomer, Ts)
    any(isequal(typemax(Int)), his) && return typemin(Int)

    return los .- his
end

function score(graph, systems)
    distribution = complement(graph, systems, 1024)
    return mean(distribution) / std(distribution)
end

function random_genome(N::Int, n::Int)
    genome = Vector{Int}(undef, n)
    genome[1] = 1
    genome[n] = N
    genome[2:2:n-1] .= rand(2:N, (n - 2) ÷ 2)
    genome[3:2:n-1] .= rand(1:N-1, (n - 2) ÷ 2)
    return genome
end

function random_population(M::Int, N::Int, n::Int)
    population = Matrix{Int}(undef, M, n)
    for i in 1:M
        population[i, :] .= random_genome(N, n)
    end
    return population
end

function mutation_operator(N::Int, MR::T) where {T<:AbstractFloat}
    let N=N, MR=MR
        function Op(genome::Vector{Int})
            n = length(genome)
            mutations = random_genome(N, n)
            R = rand(n) .< MR
            proposal = Vector{Int}(undef, n)
            for i in 1:n
                proposal[i] = R[i] ? mutations[i] : genome[i]
            end
            return proposal
        end
    end
end

function shave(graph, input, output)
    A = adjacency_matrix(graph)
    A[output, :] .= 0
    graph = SimpleDiGraph(A)
    graph = prune(graph, output)
    if outdegree(graph, input) == 0
        graph = SimpleDiGraph(size(A, 1))
    end
    return graph
end


function decoder(N::Int)
    let N=N
        function D(genome::Vector{Int})
            edges = zip(genome[1:2:end], genome[2:2:end])
            graph = SimpleDiGraph(N)
            for edge in edges
                add_edge!(graph, edge)
            end
            return shave(graph, 1, N)
        end
    end
end

function scorer(systems)
    let systems=systems
        (G::SimpleDiGraph{Int}) -> G.ne == 0 ? -Inf : score(G, systems)
    end
end

function build_identity_experiment(systems, M, n, MR)
    N = length(systems)
    return CGP(
        random_population(M, N, n),
        fill(-Inf, M),
        mutation_operator(N, MR),
        decoder(N),
        scorer(systems),
    )
end

function evolve_identity(alg, G)
    for g in 1:G
        ECCOCGP.step!(alg, 4)
        @show g, alg.scores
    end
end
