"""
    BGLSResult

Result from a BGLS periodogram computation.

# Fields
- `freqs::Vector{Float64}`: test frequencies (1/day)
- `power::Vector{Float64}`: unnormalized log-probability at each frequency
"""
struct BGLSResult
    freqs::Vector{Float64}
    power::Vector{Float64}
end

Base.iterate(r::BGLSResult, state=1) = state == 1 ? (r.freqs, 2) :
                                        state == 2 ? (r.power, 3) : nothing
Base.length(::BGLSResult) = 2

"""
    frequencies(r::BGLSResult)

Return the frequency grid.
"""
frequencies(r::BGLSResult) = r.freqs

"""
    periods(r::BGLSResult)

Return periods (1 ./ frequencies).
"""
periods(r::BGLSResult) = 1.0 ./ r.freqs

"""
    power(r::BGLSResult)

Return the log-probability spectrum.
"""
power(r::BGLSResult) = r.power

"""
    best_frequency(r::BGLSResult)

Return the frequency with the highest log-probability.
"""
best_frequency(r::BGLSResult) = r.freqs[argmax(r.power)]

"""
    best_period(r::BGLSResult)

Return the period corresponding to the highest log-probability.
"""
best_period(r::BGLSResult) = 1.0 / best_frequency(r)


"""
    SBGLSResult

Result from a Stacked BGLS periodogram computation.

# Fields
- `freqs::Vector{Float64}`: test frequencies (1/day)
- `obs_times::Vector{Float64}`: time of the last observation in each subset
- `matrix::Matrix{Float64}`: `(N - n_min + 1) × length(freqs)`, row-normalized to [0, 1]
"""
struct SBGLSResult
    freqs::Vector{Float64}
    obs_times::Vector{Float64}
    matrix::Matrix{Float64}
end

Base.iterate(r::SBGLSResult, state=1) = state == 1 ? (r.freqs, 2) :
                                          state == 2 ? (r.obs_times, 3) :
                                          state == 3 ? (r.matrix, 4) : nothing
Base.length(::SBGLSResult) = 3

"""
    frequencies(r::SBGLSResult)

Return the frequency grid.
"""
frequencies(r::SBGLSResult) = r.freqs

"""
    periods(r::SBGLSResult)

Return periods (1 ./ frequencies).
"""
periods(r::SBGLSResult) = 1.0 ./ r.freqs

"""
    obs_times(r::SBGLSResult)

Return the observation times for each row of the stacked matrix.
"""
obs_times(r::SBGLSResult) = r.obs_times

"""
    power(r::SBGLSResult)

Return the stacked periodogram matrix.
"""
power(r::SBGLSResult) = r.matrix
