function communication_barrier(rejection_rate, β_current)
    x = reverse(β_current)
    rejection_rate_cum_sum = cumsum(rejection_rate)
    y = [0.0; rejection_rate_cum_sum[1:(end-1)]]
    spline = Interpolations.interpolate(x, y, Interpolations.FritschCarlsonMonotonicInterpolation())
    Λ_fun(β) = spline(β)
    return Λ_fun
end

function update_βs(β_current, Λ_, num_replicas)
    # rejection_rate here is the average rejection rate over n_scan iters
    β = zeros(eltype(β_current), num_replicas)
    β[1] = 1.0
    β[N] = 0.0

    Λ = Λ_(1)

    for n in 2:(N-1)
        f(x) = Λ_(x) - Λ * n / (N - 1)
        β[n] = Roots.find_zero(f, (0.0, 1.0), Roots.Bisection())
    end
    return β
end
