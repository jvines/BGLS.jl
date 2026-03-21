"""
    autofrequencies(t; min_freq=0.001, max_freq=2.0, samples_per_peak=5)

Build an evenly-spaced frequency grid appropriate for a time series with
timestamps `t`. The spacing is `df = 1 / (T × samples_per_peak)` where
`T = maximum(t) - minimum(t)`.

Returns a `Vector{Float64}`.
"""
function autofrequencies(t::AbstractVector{<:Real};
                         min_freq::Real=0.001,
                         max_freq::Real=2.0,
                         samples_per_peak::Integer=5)
    T = Float64(maximum(t) - minimum(t))
    T > 0.0 || throw(ArgumentError("Time baseline must be positive, got T=$T"))

    df = 1.0 / (T * samples_per_peak)
    nf = max(1, ceil(Int, (Float64(max_freq) - Float64(min_freq)) / df))
    return collect(range(Float64(min_freq), Float64(max_freq), length=nf))
end
