function problem_from_model(::Type{Circuit}, model)
    u0 = reduce(merge, [randu0(sys) for sys in model.systems], init=Dict())
    p = Dict(:μ => 1/90, :α => 1)
    problem = DiscreteProblem(model, u0, (0.0, 3000.0), p)
    cb = CallbackSet(
        make_doubling_callback(model),
        make_species_limit_callback(9999),
    )
    return JumpProblem(
        model, problem, Direct();
        callback=cb, save_positions=(false, false),
    )
end

function problem_from_model(model::ReactionSystem)
    return problem_from_model(GeneticLogicGraph.component_type(model), model)
end

function change_input_levels(problem, params, levels)
    function find_input(p)
        return findfirst(isequal(p), parameters(problem.prob.f.sys))
    end                         
    parameter_idxs = find_input.(params)
    p = copy(problem.prob.p)
    for (i, idx) in enumerate(parameter_idxs)
        if idx != nothing
            p[idx] = levels[i]
        else
            @warn "Input not found!"
        end
    end
    return remake(problem, p=p)
end

function update_statistics(Aₖ₋₁, Qₖ₋₁, xₖ, k)
    Aₖ = Aₖ₋₁ + (xₖ - Aₖ₋₁) / k
    Qₖ = Qₖ₋₁ + (xₖ - Aₖ₋₁) * (xₖ - Aₖ)
    return Aₖ, Qₖ
end

function keep_going(Ax, Ay, Qx, Qy, k)
    k < 16 && return true
    σx = sqrt(Qx / (k - 1))
    σy = sqrt(Qy / (k - 1))
    x = Ax / σx
    y = Ay / σy
    clear_winner = abs(x - y) > 2 / √k
    return !clear_winner && k < 256
end

function test_incremental(samples)
    A = 0.0; Q = 0.0
    for i in eachindex(samples)
        A, Q = update_statistics(A, Q, samples[i], i) 
    end
    return A, sqrt(Q / (length(samples) - 1))
end

function snr_incremental(f, g)
    k = 1
    Ax, Ay, Qx, Qy = 0.0, 0.0, 0.0, 0.0
    while keep_going(Ax, Ay, Qx, Qy, k)
        sx, sy = f(), g()
        Ax, Qx = update_statistics(Ax, Qx, sx, k)
        Ay, Qy = update_statistics(Ay, Qy, sy, k)
        k = k + 1
    end
    x̂ = Ax / sqrt(Qx / (k - 1))
    ŷ = Ay / sqrt(Qy / (k - 1))
    return max(x̂, ŷ), x̂ >= ŷ
end

function sample_from_problem(problem, model, output, N)
    output_idx = findfirst(isequal(output), states(model))
    prob_func = let model=model
        function (prob, i, repeat)
            T = 1440 + rand() * 180
            udict = reduce(merge, [randu0(sys) for sys in model.systems], init=Dict())
            u0 = copy(prob.prob.u0)
            for (k, v) in collect(udict)
                u0[findfirst(isequal(k), states(model))] = v
            end
            return remake(prob, u0=u0, tspan=(0, T))
        end
    end

    sol_func = let output_idx=output_idx
        function (sol, i)
            return sol[end][output_idx], false
        end
    end

    problem = prob_func(problem, nothing, nothing)
    ensemble_problem = EnsembleProblem(
        problem;
        output_func=sol_func,
        prob_func=prob_func
    )
    return solve(ensemble_problem, SSAStepper(), EnsembleThreads(); trajectories=N)
end

