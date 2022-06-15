function communication_barrier(rejection_rate, β_current)
    x = reverse(β_current)
    rejection_rate_cum_sum = cumsum(rejection_rate)
    y = [0.0; rejection_rate_cum_sum[1:(end-1)]]
    spline = Interpolations.interpolate(x, y, Interpolations.FritschCarlsonMonotonicInterpolation())
    Λ_fun(β) = spline(β)
    return Λ_fun
end

function update_βs(β_current, rejection_rate)
    # rejection_rate here is the average rejection rate over n_scan iters
    N = length(β_current)
    β = zeros(length(β_current))
    β[1] = 0.0
    β[N] = 1.0
    Λ_ = communication_barrier(rejection_rate, β_current)
    Λ = Λ_(1)

    for n in 2:(N-1)
        f(x) = Λ_(x) - Λ / (N - 1)
        β[n] = Roots.find_zero(f, (0.0, 1.0), Roots.Bisection())
    end
    return β
end
