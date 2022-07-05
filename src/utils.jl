function communication_barrier(rejection_rate, β_current)
    x = β_current
    rejection_rate_cum_sum = cumsum(rejection_rate)
    y = [0.0; rejection_rate_cum_sum[1:(end-1)]]
    spline = Interpolations.interpolate(x, y, Interpolations.FritschCarlsonMonotonicInterpolation())
    Λ_fun(β) = spline(β)
    return Λ_fun
end

function update_βs(β_current, Λ_)
    # rejection_rate here is the average rejection rate over n_scan iters
    β_update = zeros(eltype(β_current), length(β_current))
    N = length(β_current)
    β_update[1] = 1.0
    β_update[N] = 0.0

    Λ = Λ_(1)

    for n in 2:(N-1)
        f(x) = Λ_(x) - Λ* (n / (N - 1))
        β_update[n] = Roots.find_zero(f, (max(0.0, β_update[n - 1] - 0.1), 1.0), Roots.Bisection())
    end
    return β_update
end
