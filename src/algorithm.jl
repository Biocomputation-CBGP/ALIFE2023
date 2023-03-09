"""
    `random_genome(N, n, nᵢ, nₒ)`

Generate a random genome consisting of `N` vertices, and `n ÷ 2`
edges. Vertices `1:nᵢ` are considered source vertices, and vertices
`n-nₒ:n` are considered sink vertices.

# Example
"""
function random_genome(N::T, n::T, nᵢ::T, nₒ::T) where {T<:Integer}
    genome = Vector{Int}(undef, n)

    sources  = 1:nᵢ
    sinks = N-nₒ+1:N
    inputs = 1:N-nₒ
    outputs = nᵢ+1:N
    
    genome[1:2:end] .= rand(inputs, n÷2) # the first components of edges
    genome[2:2:end] .= rand(outputs, n÷2) # the second components of edges

    for i in 1:nᵢ
        genome[2i - 1] = sources[i]
    end

    for o in 1:nₒ
        genome[end - 2o + 2] = sinks[o]
    end
    
    return genome
end

"""
    `random_population(M, N, n, nᵢ, nₒ)`

Generate a random population of M individuals, each with `N` vertices,
and `n ÷ 2` edges. Vertices `1:nᵢ` are considered source vertices, and
vertices `n-nₒ:n` are considered sink vertices.

# Example
"""
function random_population(M::Int, N::Int, n::Int, nᵢ::Int, nₒ::Int)
    population = Matrix{Int}(undef, M, n)
    for i in 1:M
        population[i, :] .= random_genome(N, n, nᵢ, nₒ)
    end
    return population
end

function mutation_operator(N::Int, MR::T, nᵢ::Int, nₒ::Int) where {T<:AbstractFloat}
    let N=N, MR=MR, nᵢ=nᵢ, nₒ=nₒ
        function Op(genome::Vector{Int})
            n = length(genome)
            mutations = random_genome(N, n, nᵢ, nₒ)
            R = rand(n) .< MR
            proposal = Vector{Int}(undef, n)
            for i in 1:n
                proposal[i] = R[i] ? mutations[i] : genome[i]
            end
            return proposal
        end
    end
end

"""
    `shave(graph, inputs, outputs)`

First remove vertices and edges from `graph` that do connect to the
sink vertices in `outputs`. Then if the resultant graph does not
contain any source vertices from `inputs`, return an empty graph,
otherwise, return the graph of vertices and edges connected to the
outputs.

# Example
"""
function shave(graph::T, inputs, outputs) where {T<:DiGraph}
    graph = prune(graph, outputs)
    if length(inputs) > 0 && any(outdegree(graph, input) == 0 for input in inputs)
        graph = T(nv(graph))
    end
    return graph
end

function decoder(N::Int, nᵢ::Int, nₒ::Int)
    let N=N, nᵢ=nᵢ, nₒ=nₒ
        function D(genome::Vector{Int})
            edges = zip(genome[1:2:end], genome[2:2:end])
            graph = SimpleDiGraph(N)
            for edge in edges
                add_edge!(graph, edge)
            end
            return shave(graph, 1:nᵢ, N-nₒ+1:N)
        end
    end
end

function decoder(N::Int)
    let N=N
        function D(genome::Vector{Int})
            edges = zip(genome[1:2:end], genome[2:2:end])
            graph = SimpleDiGraph(N)
            for edge in edges
                add_edge!(graph, edge)
            end
            return graph
        end
    end
end

struct Algorithm{T}
    population::Matrix{Int}
    scores::Vector{T}
    mutate::Function
    decode::Function
    select::Function
end

Base.size(x::Algorithm, args...) = Base.size(x.population, args...)
Base.getindex(x::Algorithm, args...) = Base.getindex(x.population, args...)

function evolve!(x::Algorithm, λ::Int, i::Int)
    parent_thing = x.decode(x[i, :])
    best_thing = x.decode(x[i, :])
    mutant_genomes = [x.mutate(x[i, :]) for _ in 1:λ]

    k = 0
    evals = 0
    
    for j in 1:λ
        mutant_thing = x.decode(mutant_genomes[j])
        if k == 0 && mutant_thing == parent_thing
            k = j
        elseif mutant_thing.ne > 0
            score, swap = x.select(mutant_thing, best_thing)
            evals = evals + 1
            if swap
                x.scores[i] = score
                best_thing = mutant_thing
                k = j
            end
        end
    end

    if k > 0
        x.population[i, :] .= mutant_genomes[k]
    end
    
    return evals
end

function evolve!(x::Algorithm, λ::Int)
    evals = 0
    for i in 1:size(x, 1)
        evals = evals + evolve!(x, λ, i)
    end
    
    _, i = findmax(x.scores)
    _, j = findmin(x.scores)
    if rand() < 0.5
        x.population[j, :] .= x.population[i, :]
        x.scores[j] = x.scores[i]
    end
    
    return evals
end
