# BayesianGLS.jl

Bayesian Generalized Lomb-Scargle periodogram for Julia.

This package implements the BGLS periodogram from
[Mortier et al. (2015)](https://doi.org/10.1051/0004-6361/201424908)
and the stacked BGLS variant from
[Mortier & Collier Cameron (2017)](https://doi.org/10.1051/0004-6361/201630092).

## Quick Start

```julia
using BayesianGLS

# Observation data
t = sort(rand(100) .* 200.0)
y = 20.0 .* sin.(2π .* 0.1 .* t) .+ randn(100) .* 2.0
err = fill(2.0, 100)

# BGLS periodogram with automatic frequency grid
result = bgls_auto(t, y, err)
println("Best period: ", best_period(result), " days")

# Stacked BGLS to see signal evolution
sresult = sbgls_auto(t, y, err; n_min=15)
```
