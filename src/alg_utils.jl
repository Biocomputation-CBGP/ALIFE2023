function decoder(N::Int, outputs::Vector{Int})
    let N=N, outputs=outputs
        function D(genome::Vector{Int})
            edges = zip(genome[1:2:end], genome[2:2:end])
            graph = SimpleDiGraph(N)
            for edge in edges
                add_edge!(graph, edge)
            end
            return shave(graph, outputs)
        end
    end
end

function decoder(N::Int, inputs::Vector{Int}, outputs::Vector{Int})
    let N=N
        function D(genome::Vector{Int})
            edges = zip(genome[1:2:end], genome[2:2:end])
            graph = SimpleDiGraph(N)
            for edge in edges
                add_edge!(graph, edge)
            end
            return shave(graph, inputs, outputs)
        end
    end
end
