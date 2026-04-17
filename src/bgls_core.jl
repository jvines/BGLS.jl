"""
    bgls(t, y, err, freqs) -> BGLSResult

Compute the Bayesian Generalized Lomb-Scargle log-probability spectrum.

Implementation of Mortier et al. (2015):
"BGLS: A Bayesian formalism for the generalised Lomb-Scargle periodogram"
A&A 573, A101

# Arguments
- `t::AbstractVector{<:Real}`: observation times (BJD or similar)
- `y::AbstractVector{<:Real}`: observed values (RV, etc.)
- `err::AbstractVector{<:Real}`: measurement uncertainties (1σ)
- `freqs::AbstractVector{<:Real}`: test frequencies (1/day)

# Returns
- `BGLSResult` with fields `freqs` and `power` (unnormalized log-probability)

# Algorithm (Mortier et al. 2015, Eq. 7-12)
For each frequency ω = 2πf, build the weighted design matrix
D = [1, sin(ωt), cos(ωt)] and compute the marginal likelihood by
analytically integrating over the linear parameters (offset, amplitudes).
"""
function bgls(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
              err::AbstractVector{<:Real}, freqs::AbstractVector{<:Real})
    N = length(t)
    nf = length(freqs)
    logp = Vector{Float64}(undef, nf)

    # Precompute weights
    w = 1.0 ./ (err .^ 2)
    W = sum(w)

    # Weighted mean and centered data
    Y_w = dot(w, y) / W
    y_centered = y .- Y_w

    # Null-model log-likelihood (flat model = weighted mean)
    chi2_null = dot(w, y_centered .^ 2)
    log_L0 = -0.5 * chi2_null

    @inbounds for k in 1:nf
        omega = 2.0 * π * freqs[k]

        # Trig vectors
        s = sin.(omega .* t)
        c = cos.(omega .* t)

        # Weighted sums (Mortier Eq. 8-11)
        Ss  = dot(w, s)
        Sc  = dot(w, c)
        YS  = dot(w, y_centered .* s)
        YC  = dot(w, y_centered .* c)
        SS  = dot(w, s .* s)
        CC  = dot(w, c .* c)
        CS  = dot(w, c .* s)

        # Determinant of the normal equations sub-block
        D = SS * CC - CS * CS

        if D < 1e-30
            # Degenerate frequency — skip
            logp[k] = log_L0
            continue
        end

        # Amplitude estimators (Eq. 12)
        delta_chi2 = (YS * YS * CC - 2.0 * YS * YC * CS + YC * YC * SS) / D

        # Log marginal likelihood (up to constant shared across frequencies):
        #   logp = 0.5*Δχ² - 0.5*log(D)
        logp[k] = 0.5 * delta_chi2 - 0.5 * log(D)
    end

    return BGLSResult(collect(Float64, freqs), logp)
end


"""
    bgls_auto(t, y, err; min_freq=0.001, max_freq=2.0, samples_per_peak=5) -> BGLSResult

Compute BGLS with an automatically generated frequency grid.
"""
function bgls_auto(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                   err::AbstractVector{<:Real};
                   min_freq::Real=0.001,
                   max_freq::Real=2.0,
                   samples_per_peak::Integer=5)
    freqs = autofrequencies(t; min_freq, max_freq, samples_per_peak)
    return bgls(t, y, err, freqs)
end
