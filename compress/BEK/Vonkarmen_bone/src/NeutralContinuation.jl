module NeutralContinuation

using LinearAlgebra

export continue_modes, find_neutral_R, refined_beta_step, recovered_beta_step

"""
    refined_beta_step(current_step, minimum_step; enabled=true, factor=0.5)

Return a smaller continuation step after a failed neutral correction. `nothing`
means that refinement is disabled or that the next reduction would cross the
configured minimum. This gives callers a finite, explicit stopping condition.
"""
function refined_beta_step(
    current_step::Real, minimum_step::Real;
    enabled::Bool=true, factor::Real=0.5,
)
    current_step > 0 || throw(ArgumentError("current_step must be positive"))
    minimum_step > 0 || throw(ArgumentError("minimum_step must be positive"))
    minimum_step <= current_step || throw(ArgumentError(
        "minimum_step must not exceed current_step",
    ))
    0 < factor < 1 || throw(ArgumentError("factor must lie in (0,1)"))
    enabled || return nothing

    proposed = current_step * factor
    tolerance = 10eps(Float64) * max(abs(current_step), abs(minimum_step), 1.0)
    proposed + tolerance >= minimum_step || return nothing
    return max(Float64(proposed), Float64(minimum_step))
end

"""Increase a reduced beta step after enough consecutive successful points."""
function recovered_beta_step(
    current_step::Real, nominal_step::Real,
    consecutive_successes::Integer, required_successes::Integer;
    factor::Real=2.0,
)
    0 < current_step <= nominal_step || throw(ArgumentError(
        "current_step must lie in (0,nominal_step]",
    ))
    required_successes >= 1 || throw(ArgumentError(
        "required_successes must be positive",
    ))
    consecutive_successes >= 0 || throw(ArgumentError(
        "consecutive_successes must be nonnegative",
    ))
    factor > 1 || throw(ArgumentError("factor must exceed one"))
    if current_step == nominal_step || consecutive_successes < required_successes
        return Float64(current_step)
    end
    return min(Float64(nominal_step), Float64(current_step * factor))
end

function mode_overlap(left, right)
    denominator = norm(left) * norm(right)
    denominator > eps(Float64) || return 0.0
    return abs(dot(left, right)) / denominator
end

function phase_align!(vector, seed)
    phase = dot(seed, vector)
    abs(phase) > eps(Float64) && (vector .*= conj(phase) / abs(phase))
    return vector
end

function candidate_assignment(previous_vectors, values, vectors, mode_count)
    candidate_count = min(length(values), size(vectors, 2))
    candidate_count >= mode_count || return nothing
    overlaps = [
        mode_overlap(view(previous_vectors, :, old), view(vectors, :, new))
        for old in 1:mode_count, new in 1:candidate_count
    ]

    if mode_count == 1
        candidate = argmax(view(overlaps, 1, :))
        return [candidate], [overlaps[1, candidate]]
    end

    best_score = -Inf
    best_assignment = nothing
    for first in 1:candidate_count, second in 1:candidate_count
        first == second && continue
        score = overlaps[1, first] + overlaps[2, second]
        if score > best_score
            best_score = score
            best_assignment = [first, second]
        end
    end
    best_assignment === nothing && return nothing
    selected_overlaps = [overlaps[index, best_assignment[index]] for index in 1:2]
    return best_assignment, selected_overlaps
end

function finalize_modes(previous_vectors, values, vectors, assignment)
    mode_count = length(assignment)
    next_values = Vector{ComplexF64}(undef, mode_count)
    next_vectors = Matrix{ComplexF64}(undef, size(previous_vectors, 1), mode_count)
    overlaps = Vector{Float64}(undef, mode_count)
    for mode_index in 1:mode_count
        candidate = assignment[mode_index]
        seed = view(previous_vectors, :, mode_index)
        vector = ComplexF64.(vectors[:, candidate])
        phase_align!(vector, seed)
        next_values[mode_index] = values[candidate]
        next_vectors[:, mode_index] .= vector
        overlaps[mode_index] = mode_overlap(seed, vector)
    end
    return next_values, next_vectors, overlaps
end

"""
    continue_modes(solve_at, R, beta, previous_values, previous_vectors; kwargs...)

Continue at most two eigenpairs. The fast path runs one IAR solve for all modes
using a combined seed and assigns candidates by maximum total eigenvector
overlap. Separate warm-started solves are used only when that assignment fails.
"""
function continue_modes(
    solve_at, R, beta, previous_values, previous_vectors;
    min_overlap::Real=0.60,
    duplicate_tol::Real=1.0e-7,
)
    mode_count = min(length(previous_values), size(previous_vectors, 2), 2)
    mode_count > 0 || error("No eigenpair is available for continuation")
    size(previous_vectors, 1) > 0 || error("Tracked eigenvectors are empty")

    combined_seed = vec(sum(previous_vectors[:, 1:mode_count]; dims=2))
    combined_norm = norm(combined_seed)
    if combined_norm <= eps(Float64)
        combined_seed .= view(previous_vectors, :, 1)
        combined_norm = norm(combined_seed)
    end
    combined_seed ./= combined_norm
    combined_target = sum(previous_values[1:mode_count]) / mode_count

    try
        values, vectors = solve_at(
            R, beta, combined_target, mode_count; v0=combined_seed,
        )
        assignment = candidate_assignment(
            previous_vectors, values, vectors, mode_count,
        )
        if assignment !== nothing
            indices, selected_overlaps = assignment
            distinct = mode_count == 1 ||
                abs(values[indices[1]] - values[indices[2]]) >= duplicate_tol
            if distinct && minimum(selected_overlaps) >= min_overlap
                return finalize_modes(previous_vectors, values, vectors, indices)
            end
        end
    catch exception
        exception isa InterruptException && rethrow()
    end

    next_values = Vector{ComplexF64}(undef, mode_count)
    next_vectors = Matrix{ComplexF64}(undef, size(previous_vectors, 1), mode_count)
    overlaps = zeros(Float64, mode_count)
    for mode_index in 1:mode_count
        seed = view(previous_vectors, :, mode_index)
        values, vectors = solve_at(
            R, beta, previous_values[mode_index], 1; v0=seed,
        )
        isempty(values) && error(
            "IAR returned no eigenvalue for tracked mode $mode_index",
        )
        candidate = 1
        first_overlap = mode_overlap(seed, view(vectors, :, candidate))
        duplicate = mode_index > 1 && any(
            abs(values[candidate] - next_values[index]) < duplicate_tol
            for index in 1:mode_index-1
        )

        if duplicate || first_overlap < min_overlap
            values, vectors = solve_at(
                R, beta, previous_values[mode_index], 2; v0=seed,
            )
            candidate_overlaps = [
                mode_overlap(seed, view(vectors, :, index))
                for index in eachindex(values)
            ]
            available = [
                index for index in eachindex(values) if all(
                    abs(values[index] - next_values[old]) >= duplicate_tol
                    for old in 1:mode_index-1
                )
            ]
            candidate = isempty(available) ? argmax(candidate_overlaps) :
                available[argmax(candidate_overlaps[available])]
        end

        vector = ComplexF64.(vectors[:, candidate])
        phase_align!(vector, seed)
        next_values[mode_index] = values[candidate]
        next_vectors[:, mode_index] .= vector
        overlaps[mode_index] = mode_overlap(seed, vector)
        if overlaps[mode_index] < min_overlap
            @warn "Low eigenvector overlap during continuation" R beta mode_index overlap=overlaps[mode_index]
        end
    end
    return next_values, next_vectors, overlaps
end

function evaluate_state(
    solve_at, R, beta, seed_values, seed_vectors, active_index;
    min_overlap, duplicate_tol, on_evaluation,
)
    values, vectors, overlaps = continue_modes(
        solve_at, R, beta, seed_values, seed_vectors;
        min_overlap=min_overlap, duplicate_tol=duplicate_tol,
    )
    active_index <= length(values) || error("The active tracked mode was lost")
    state = (
        R=Float64(R), beta=Float64(beta), values=values, vectors=vectors,
        overlaps=overlaps, active_index=active_index,
        residual=imag(values[active_index]),
    )
    on_evaluation(state)
    return state
end

function try_evaluate(args...; kwargs...)
    try
        return evaluate_state(args...; kwargs...)
    catch exception
        exception isa InterruptException && rethrow()
        return nothing
    end
end

function refine_bracket(
    solve_at, beta, left, right;
    neutral_tol, R_tol, max_refine, min_overlap, duplicate_tol,
    on_evaluation,
)
    left.R <= right.R || ((left, right) = (right, left))
    left.residual * right.residual <= 0 || error(
        "Neutral refinement requires a sign-changing bracket",
    )

    for iteration in 1:max_refine
        abs(left.residual) <= neutral_tol && return merge(left, (refine_iterations=iteration - 1,))
        abs(right.residual) <= neutral_tol && return merge(right, (refine_iterations=iteration - 1,))
        width = right.R - left.R
        width <= R_tol && break

        denominator = right.residual - left.residual
        trial_R = abs(denominator) > eps(Float64) ?
            (left.R * right.residual - right.R * left.residual) / denominator :
            0.5 * (left.R + right.R)
        margin = 0.1 * width
        trial_R = clamp(trial_R, left.R + margin, right.R - margin)
        seed = abs(trial_R - left.R) <= abs(right.R - trial_R) ? left : right
        trial = try_evaluate(
            solve_at, trial_R, beta, seed.values, seed.vectors,
            seed.active_index;
            min_overlap=min_overlap, duplicate_tol=duplicate_tol,
            on_evaluation=on_evaluation,
        )
        if trial === nothing
            trial_R = 0.5 * (left.R + right.R)
            seed = abs(trial_R - left.R) <= abs(right.R - trial_R) ? left : right
            trial = evaluate_state(
                solve_at, trial_R, beta, seed.values, seed.vectors,
                seed.active_index;
                min_overlap=min_overlap, duplicate_tol=duplicate_tol,
                on_evaluation=on_evaluation,
            )
        end
        abs(trial.residual) <= neutral_tol && return merge(trial, (refine_iterations=iteration,))
        if left.residual * trial.residual <= 0
            right = trial
        else
            left = trial
        end
    end

    best = abs(left.residual) <= abs(right.residual) ? left : right
    abs(best.residual) <= neutral_tol || error(
        "Neutral refinement stopped at residual=$(best.residual), above $neutral_tol",
    )
    return merge(best, (refine_iterations=max_refine,))
end

"""
    find_neutral_R(solve_at, beta, R_guess, seed_values, seed_vectors,
                   active_index; kwargs...)

Locate `imag(alpha_active)=0` near a predicted Reynolds number. The search uses
warm-started mode continuation, a local sign-changing bracket, and safeguarded
secant refinement. It never accepts a point above `neutral_tol`.
"""
function find_neutral_R(
    solve_at, beta, R_guess, seed_values, seed_vectors, active_index;
    R_step::Real=0.5,
    preferred_direction::Integer=0,
    max_scan_steps::Integer=80,
    max_refine::Integer=30,
    neutral_tol::Real=1.0e-7,
    R_tol::Real=1.0e-4,
    R_bounds=(1.0e-6, 700.0),
    max_R_deviation::Real=Inf,
    min_overlap::Real=0.60,
    duplicate_tol::Real=1.0e-7,
    on_evaluation=(state -> nothing),
)
    R_step > 0 || throw(ArgumentError("R_step must be positive"))
    preferred_direction in (-1, 0, 1) || throw(ArgumentError(
        "preferred_direction must be -1, 0, or 1",
    ))
    lower_R, upper_R = R_bounds
    lower_R < R_guess < upper_R || throw(ArgumentError("R_guess lies outside R_bounds"))

    center = evaluate_state(
        solve_at, R_guess, beta, seed_values, seed_vectors, active_index;
        min_overlap=min_overlap, duplicate_tol=duplicate_tol,
        on_evaluation=on_evaluation,
    )
    abs(center.residual) <= neutral_tol && return merge(center, (scan_iterations=0, refine_iterations=0,))

    roots = NamedTuple[]
    neighbors = NamedTuple[]
    for direction in (-1, 1)
        next_R = R_guess + direction * R_step
        lower_R < next_R < upper_R || continue
        neighbor = try_evaluate(
            solve_at, next_R, beta, center.values, center.vectors,
            active_index;
            min_overlap=min_overlap, duplicate_tol=duplicate_tol,
            on_evaluation=on_evaluation,
        )
        neighbor === nothing && continue
        neighbor = merge(neighbor, (direction=direction, scan_iterations=1,))
        if center.residual * neighbor.residual <= 0
            root = refine_bracket(
                solve_at, beta, center, neighbor;
                neutral_tol=neutral_tol, R_tol=R_tol,
                max_refine=max_refine, min_overlap=min_overlap,
                duplicate_tol=duplicate_tol, on_evaluation=on_evaluation,
            )
            abs(root.R - R_guess) <= R_step && return merge(root, (scan_iterations=1,))
            push!(roots, merge(root, (scan_iterations=1,)))
        else
            push!(neighbors, neighbor)
        end
    end
    isempty(neighbors) && isempty(roots) && error(
        "Both initial R-continuation directions failed at beta=$beta",
    )

    sort!(neighbors; by=neighbor -> (
        preferred_direction == 0 || neighbor.direction == preferred_direction ? 0 : 1,
        abs(neighbor.residual),
    ))
    for neighbor in neighbors
        previous = center
        current = neighbor
        direction = neighbor.direction
        for iteration in 2:max_scan_steps
            next_R = current.R + direction * R_step
            lower_R < next_R < upper_R || break
            next = try_evaluate(
                solve_at, next_R, beta, current.values, current.vectors,
                active_index;
                min_overlap=min_overlap, duplicate_tol=duplicate_tol,
                on_evaluation=on_evaluation,
            )
            next === nothing && break
            if current.residual * next.residual <= 0
                root = refine_bracket(
                    solve_at, beta, current, next;
                    neutral_tol=neutral_tol, R_tol=R_tol,
                    max_refine=max_refine, min_overlap=min_overlap,
                    duplicate_tol=duplicate_tol, on_evaluation=on_evaluation,
                )
                push!(roots, merge(root, (scan_iterations=iteration,)))
                break
            end
            previous, current = current, next
            if iteration >= 5 &&
               abs(current.residual) > abs(previous.residual) > abs(center.residual)
                break
            end
        end
    end
    isempty(roots) && error(
        "Failed to bracket a neutral point at beta=$beta from R_guess=$R_guess",
    )
    root = roots[argmin([abs(candidate.R - R_guess) for candidate in roots])]
    deviation = abs(root.R - R_guess)
    deviation <= max_R_deviation || error(
        "Nearest neutral root is $deviation away from R_guess=$R_guess",
    )
    abs(root.residual) <= neutral_tol || error(
        "Corrected neutral residual $(root.residual) exceeds $neutral_tol",
    )
    return root
end

end
