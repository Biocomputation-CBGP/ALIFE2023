function shave(graph::T, inputs, outputs) where {T<:DiGraph}
    graph = prune(graph, outputs)
    if all(outdegree(graph, input) == 0 for input in inputs)
        graph = SimpleDiGraph(graph.nv)
    end
    return graph
end
