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

# @time bscore = Inverter.benchmark_score(2048)
# @show bscore
# alg = build_inverter_experiment(1, 7)

function cgp_convergence_iterations(N)
    i = 0
    alg = build_inverter_experiment(1, N)
    while true
        i = i + 1
        evolve!(alg, 4)
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

function record_averages(N)
    Ns = collect(3:N)
    iterations = Matrix{Int}(undef, 64, length(Ns))
    for i in eachindex(Ns)
        @show i
        iterations[:, i] .= average_random_graph_convergence(64, Ns[i])
    end
    print(iterations)
    CSV.write("iterations.csv", Tables.table(iterations), writeheader=false)
end

function load_averages(fn)
    CSV.File(fn) |> Tables.matrix
end

function record_cgp_averages(N)
    Ns = collect(3:N)
    iterations = Matrix{Int}(undef, 64, length(Ns))
    for i in eachindex(Ns)
        @show Ns[i]
        Threads.@threads for j in 1:64
            iterations[j, i] = cgp_convergence_iterations(Ns[i])
        end
        @show mean(iterations[:, i])
    end
    print(iterations)
    CSV.write("cgp-iterations.csv", Tables.table(iterations .* 4), writeheader=false)
end

function timecourse_of_benchmark()
    inverter = Inverter.benchmark_model()
    problem = Inverter.benchmark_problem()
    lprob = change_input_levels(problem, [(@nonamespace inverter.LacI).λ], [10 * log(2) / 90])
    hprob = change_input_levels(problem, [(@nonamespace inverter.LacI).λ], [100 * log(2) / 90])

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
    x = repeat([3, 4, 5, 6], inner=63)
    plt = plot()
    μs = mean(A', dims=2)
    groupedboxplot!(plt, vcat(x, x), vcat(A[:], B[:]), group=vcat(zeros(Int, length(A)), ones(Int, length(A))), markersize=3, markershape=:x, markerstrokewidth=0, outliers=false, labels=["Random" "Algorithm"])
    scatter!(plt, x .- rand(Normal(0.25, 0.05), 63 * 4), A[:], color=:black, markersize=2, label=false)
    scatter!(plt, x .+ rand(Normal(0.25, 0.05), 63 * 4), B[:], color=:black, markersize=2, label=false)
    plot!(plt, minorgrid=false, minorticks=false, yscale=:log10)
    return plt
end

function runner(N)
    record_averages(N)
    record_cgp_averages(N)
end
