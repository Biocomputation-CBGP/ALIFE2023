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
using Distributions
using LaTeXStrings
using GeneticLogicGraph

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
        @named model = Circuit(alg.decode(alg.population[1,:]), Inverter.components())
        problem = problem_from_model(model)
        d = Inverter.complement(problem, model, 1024)
        CSV.write("not-$g-complement.csv", Tables.table(d), writeheader=false)
    end
    return scores
end

@time bscore = Inverter.benchmark_score(2048)
@show bscore
alg = build_inverter_experiment(1, 7)

function cgp_convergence_iterations(N)
    i = 0
    alg = build_inverter_experiment(1, N)
    while true
        i = i + evolve!(alg, 4)
        g = alg.decode(alg.population[1, :])
        if g.ne > 0
            score, swap = alg.select(g, Inverter.benchmark_graph(N))
            swap && return i
        end
    end
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
    Threads.@threads for i in eachindex(is)
        is[i] = random_graph_convergence_iterations(N)
    end
    return is
end

function record_averages(n, Ns, fn)
    iterations = Matrix{Int}(undef, n, length(Ns))
    for i in eachindex(Ns)
        @show Ns[i]
        Threads.@threads for j in 1:n
            iterations[j, i] = random_graph_convergence_iterations(Ns[i])
        end
    end
    CSV.write("iterations-$(fn).csv", Tables.table(iterations), writeheader=false)
end

function load_averages(fn)
    CSV.File(fn, header=false) |> Tables.matrix
end

function record_cgp_averages(n, Ns, fn)
    iterations = Matrix{Int}(undef, n, length(Ns))
    for i in eachindex(Ns)
        @show Ns[i]
        Threads.@threads for j in 1:n
            iterations[j, i] = cgp_convergence_iterations(Ns[i])
        end
    end
    CSV.write("cgp-iterations-$(fn).csv", Tables.table(iterations), writeheader=false)
end

function timecourse_of_benchmark()
    inverter = Inverter.benchmark_model()
    problem = Inverter.benchmark_problem()
    lprob = change_input_levels(problem, [(@nonamespace inverter.LacI).??], [10 * log(2) / 90])
    hprob = change_input_levels(problem, [(@nonamespace inverter.LacI).??], [100 * log(2) / 90])

    Plots.theme(:dao)
    plt = plot()
    plot!(plt, solve(lprob, SSAStepper(), saveat=1), idxs=[8], label="Low Input")
    plot!(plt, solve(hprob, SSAStepper(), saveat=1), idxs=[8], label="High Input")
    plot!(
        plt,
        xlabel="Time (hours)",
        ylabel="Abundance of output",
        legend=:best,
        size=(350, 275),
        tickfontsize=8,
        legendfontsize=8,
        guidefontsize=8,
        xticks=([720, 1440, 2160, 2880], ["12", "24", "36", "48"])
    )
    return plt
end

function complement_of_benchmark()
    inverter = Inverter.benchmark_model()
    problem = Inverter.benchmark_problem()

    comp = Inverter.benchmark_distribution(1024)
    Plots.theme(:dao)
    plt = plot()
    histogram!(plt, comp)
    plot!(
        plt,
        xlabel="High - Low output",
        ylabel="Frequency",
        legend=false,
        size=(350, 275),
        tickfontsize=8,
        legendfontsize=8,
        guidefontsize=8,
    )
    return plt
end

function convergence_to_benchmark()
    ifn = "/home/lewis/sauces/julia/GeneticLogicGraph/src/iterations.csv"
    cfn = "/home/lewis/sauces/julia/GeneticLogicGraph/src/cgp-iterations.csv"
    A = load_averages(ifn)
    B = load_averages(cfn)
    Plots.theme(:dao)
    x = repeat([3, 4, 5, 6, 7, 8], inner=64)
    plt = plot()
    A??s = mean(A', dims=2)
    B??s = mean(B', dims=2)
    A??s = std(A', dims=2)
    B??s = std(B', dims=2)
    A??s = A??s / ???64
    B??s = B??s / ???64

    @show A??s, A??s, A??s
    @show B??s, B??s, B??s

    groupedboxplot!(
        plt,
        vcat(x, x),
        vcat(A[:], B[:]),
        group=vcat(zeros(Int, length(x)), ones(Int, length(x))),
        markersize=3,
        markershape=:x,
        markerstrokewidth=0,
        outliers=false,
        labels=["Random search" "Algorithm search"]
    )
    scatter!(
        plt,
        x .- rand(Normal(0.25, 0.05), length(x)),
        A[:],
        color=:black,
        markersize=1.5,
        markeralpha=0.5,
        label=false
    )
    scatter!(
        plt,
        x .+ rand(Normal(0.25, 0.05), length(x)),
        B[:],
        color=:black,
        markersize=1.5,
        markeralpha=0.5,
        label=false
    )
    plot!(
        plt,
        minorgrid=false,
        minorticks=false,
        yscale=:log10,
        legend=:best,
        size=(350, 275),
        tickfontsize=8,
        legendfontsize=8,
        guidefontsize=8,
        xlabel=L"Graph size $|V|$",
        ylabel="Function evaluations",
    )
    
    return plt
end

function runner(N, fn)
    record_averages(8, N, fn)
    record_cgp_averages(8, N, fn)
end
