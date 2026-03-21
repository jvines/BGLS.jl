"""
    sbgls(t, y, err, freqs; n_min=10) -> SBGLSResult

Compute the Stacked BGLS periodogram.

Implementation following Mortier & Collier Cameron (2017):
"The stacked Bayesian general Lomb-Scargle periodogram"
A&A 601, A110

Computes BGLS on successively longer subsets of the time series,
producing a 2D map of normalized log-probability vs frequency and
observation count. Real planetary signals grow monotonically as
data accumulates; aliases and activity signals appear/disappear.

# Arguments
- `t::AbstractVector{<:Real}`: observation times
- `y::AbstractVector{<:Real}`: observed values
- `err::AbstractVector{<:Real}`: measurement uncertainties (1σ)
- `freqs::AbstractVector{<:Real}`: test frequencies (1/day)

# Keyword Arguments
- `n_min::Int=10`: minimum number of observations for the first row

# Returns
- `SBGLSResult` with fields `freqs`, `obs_times`, and `matrix`
  where `matrix` is `(N - n_min + 1) × length(freqs)`, row-normalized to [0, 1]
"""
function sbgls(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
               err::AbstractVector{<:Real}, freqs::AbstractVector{<:Real};
               n_min::Int=10)
    N = length(t)
    nf = length(freqs)

    if N < n_min
        throw(ArgumentError("Need at least n_min=$n_min observations, got N=$N"))
    end

    # Sort by time
    perm = sortperm(t)
    t_sorted = collect(Float64, t[perm])
    y_sorted = collect(Float64, y[perm])
    err_sorted = collect(Float64, err[perm])
    freqs_f64 = collect(Float64, freqs)

    n_rows = N - n_min + 1
    matrix = Matrix{Float64}(undef, n_rows, nf)

    @inbounds for row in 1:n_rows
        k = n_min + row - 1  # number of observations in this subset

        result = bgls(t_sorted[1:k], y_sorted[1:k], err_sorted[1:k], freqs_f64)
        logp = result.power

        # Row-normalize to [0, 1]
        lo = minimum(logp)
        hi = maximum(logp)
        span = hi - lo
        if span > 0.0
            @simd for j in 1:nf
                matrix[row, j] = (logp[j] - lo) / span
            end
        else
            matrix[row, :] .= 0.0
        end
    end

    # obs_times: the time of the last observation in each subset
    obs_times_vec = t_sorted[n_min:N]

    return SBGLSResult(freqs_f64, obs_times_vec, matrix)
end


"""
    sbgls_auto(t, y, err; min_freq=0.001, max_freq=2.0, samples_per_peak=5,
               n_min=10, max_frequencies=2000) -> SBGLSResult

Compute Stacked BGLS with an automatically generated frequency grid.

The frequency grid is capped at `max_frequencies` to keep the 2D result
matrix at a manageable size.
"""
function sbgls_auto(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                    err::AbstractVector{<:Real};
                    min_freq::Real=0.001,
                    max_freq::Real=2.0,
                    samples_per_peak::Integer=5,
                    n_min::Int=10,
                    max_frequencies::Int=2000)
    T = Float64(maximum(t) - minimum(t))
    T > 0.0 || throw(ArgumentError("Time baseline must be positive, got T=$T"))

    df = 1.0 / (T * samples_per_peak)
    nf = max(1, ceil(Int, (Float64(max_freq) - Float64(min_freq)) / df))
    nf = min(nf, max_frequencies)

    freqs = collect(range(Float64(min_freq), Float64(max_freq), length=nf))

    return sbgls(t, y, err, freqs; n_min=n_min)
end
