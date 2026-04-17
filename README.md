# BayesianGLS.jl

Bayesian Generalized Lomb-Scargle periodogram for Julia, implementing the
methods of [Mortier et al. (2015)](https://doi.org/10.1051/0004-6361/201424908)
and the stacked variant from [Mortier & Collier Cameron (2017)](https://doi.org/10.1051/0004-6361/201630092).

## Features

- **BGLS**: Bayesian periodogram that analytically marginalizes over amplitude, phase, and offset
- **Stacked BGLS (sBGLS)**: 2D map of signal evolution as observations accumulate — real signals grow monotonically, aliases and activity signals appear/disappear
- Zero non-stdlib dependencies (only `LinearAlgebra`)
- Accepts any `AbstractVector{<:Real}` inputs

## Installation

```julia
using Pkg
Pkg.add("BayesianGLS")
```

## Quick Start

```julia
using BayesianGLS

# Synthetic RV data
t = sort(rand(100) .* 200.0)
y = 20.0 .* sin.(2π .* 0.1 .* t) .+ randn(100) .* 2.0
err = fill(2.0, 100)

# BGLS periodogram
result = bgls_auto(t, y, err)
println("Best period: ", best_period(result), " days")

# Stacked BGLS
sresult = sbgls_auto(t, y, err; n_min=15)
# sresult.matrix is a 2D array: rows = accumulating observations, cols = frequencies
```

## API

### BGLS

| Function | Description |
|----------|-------------|
| `bgls(t, y, err, freqs)` | Compute BGLS for a given frequency grid |
| `bgls_auto(t, y, err; min_freq, max_freq, samples_per_peak)` | BGLS with automatic frequency grid |

### Stacked BGLS

| Function | Description |
|----------|-------------|
| `sbgls(t, y, err, freqs; n_min)` | Compute stacked BGLS for a given frequency grid |
| `sbgls_auto(t, y, err; min_freq, max_freq, samples_per_peak, n_min, max_frequencies)` | Stacked BGLS with automatic frequency grid |

### Utilities

| Function | Description |
|----------|-------------|
| `autofrequencies(t; min_freq, max_freq, samples_per_peak)` | Build a frequency grid from time baseline |
| `best_frequency(result)` | Frequency of highest log-probability |
| `best_period(result)` | Period of highest log-probability |
| `frequencies(result)` | Frequency grid |
| `periods(result)` | Period grid (1 ./ frequencies) |
| `power(result)` | Log-probability spectrum (or matrix for sBGLS) |

All result types support tuple destructuring:
```julia
freqs, logp = bgls(t, y, err, freqs)
freqs, obs_times, matrix = sbgls(t, y, err, freqs)
```

## Citation

If you use this package, please cite:

> Mortier, A., Faria, J. P., Correia, C. M., Santerne, A., & Santos, N. C.
> (2015). BGLS: A Bayesian formalism for the generalised Lomb-Scargle
> periodogram. *Astronomy & Astrophysics*, 573, A101.

> Mortier, A. & Collier Cameron, A. (2017). The stacked Bayesian general
> Lomb-Scargle periodogram. *Astronomy & Astrophysics*, 601, A110.

```bibtex
@article{Mortier2015,
  author  = {Mortier, A. and Faria, J. P. and Correia, C. M. and Santerne, A. and Santos, N. C.},
  title   = {{BGLS}: A {B}ayesian formalism for the generalised {L}omb-{S}cargle periodogram},
  journal = {Astronomy \& Astrophysics},
  volume  = {573},
  pages   = {A101},
  year    = {2015},
  doi     = {10.1051/0004-6361/201424908}
}

@article{Mortier2017,
  author  = {Mortier, A. and Collier Cameron, A.},
  title   = {The stacked {B}ayesian general {L}omb-{S}cargle periodogram},
  journal = {Astronomy \& Astrophysics},
  volume  = {601},
  pages   = {A110},
  year    = {2017},
  doi     = {10.1051/0004-6361/201630092}
}
```

## See Also

- [LombScargle.jl](https://github.com/JuliaAstro/LombScargle.jl) -- Classical Lomb-Scargle periodogram
- [L1Periodograms.jl](https://github.com/jvines/L1Periodograms.jl) -- LASSO-based sparse periodogram (Hara et al. 2017)

## License

MIT
