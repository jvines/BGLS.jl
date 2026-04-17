using Test
using BayesianGLS
const BGLS = BayesianGLS
using Random
using Statistics: mean

function make_single_sinusoid(; seed=42)
    rng = MersenneTwister(seed)
    n_obs = 100; true_freq = 0.1; amplitude = 20.0; noise_level = 2.0
    t = sort(rand(rng, n_obs) .* 200.0)
    y = amplitude .* sin.(2π .* true_freq .* t) .+ randn(rng, n_obs) .* noise_level
    err = fill(noise_level, n_obs)
    return t, y, err, true_freq
end

function make_noise_only(; seed=99)
    rng = MersenneTwister(seed)
    n_obs = 80; noise_level = 5.0
    t = sort(rand(rng, n_obs) .* 100.0)
    y = randn(rng, n_obs) .* noise_level
    err = fill(noise_level, n_obs)
    return t, y, err
end

@testset "BGLS.jl" begin

    @testset "Single sinusoid recovery" begin
        t, y, err, true_freq = make_single_sinusoid()
        freqs = collect(range(0.01, 0.5, length=500))
        result = bgls(t, y, err, freqs)
        peak_freq = best_frequency(result)
        df = freqs[2] - freqs[1]
        @test abs(peak_freq - true_freq) <= df
    end

    @testset "Noise-only: low dynamic range" begin
        t, y, err = make_noise_only()
        result = bgls_auto(t, y, err)
        pwr = power(result)
        @test maximum(pwr) - minimum(pwr) < 50.0
    end

    @testset "Output shapes and types" begin
        t, y, err, _ = make_single_sinusoid()
        freqs = collect(range(0.01, 0.5, length=200))
        result = bgls(t, y, err, freqs)
        @test length(frequencies(result)) == 200
        @test length(power(result)) == 200
        @test all(isfinite, power(result))
    end

    @testset "Destructuring" begin
        t, y, err, _ = make_single_sinusoid()
        freqs = collect(range(0.01, 0.5, length=100))
        f, p = bgls(t, y, err, freqs)
        @test length(f) == 100
        @test length(p) == 100
    end

    @testset "autofrequencies" begin
        t = [0.0, 50.0, 100.0]
        freqs = autofrequencies(t; samples_per_peak=5)
        @test length(freqs) > 0
        @test freqs[1] >= 0.001
        @test freqs[end] <= 2.0
        @test issorted(freqs)
    end

    @testset "periods accessor" begin
        t, y, err, _ = make_single_sinusoid()
        result = bgls_auto(t, y, err)
        p = periods(result)
        f = frequencies(result)
        @test all(p .≈ 1.0 ./ f)
    end

    @testset "best_period consistency" begin
        t, y, err, true_freq = make_single_sinusoid()
        freqs = collect(range(0.01, 0.5, length=500))
        result = bgls(t, y, err, freqs)
        @test best_period(result) ≈ 1.0 / best_frequency(result)
    end

    @testset "Edge cases" begin
        @testset "Minimum data (3 points)" begin
            t = [0.0, 1.0, 5.0]
            y = [1.0, -1.0, 0.5]
            err = [0.1, 0.1, 0.1]
            freqs = collect(range(0.01, 0.5, length=20))
            result = bgls(t, y, err, freqs)
            @test length(power(result)) == 20
            @test all(isfinite, power(result))
        end

        @testset "Single frequency" begin
            t, y, err, _ = make_single_sinusoid()
            result = bgls(t, y, err, [0.1])
            @test length(power(result)) == 1
        end

        @testset "Zero baseline throws" begin
            @test_throws ArgumentError bgls_auto([1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
        end

        @testset "View inputs" begin
            t, y, err, _ = make_single_sinusoid()
            freqs = collect(range(0.01, 0.5, length=200))
            result = bgls(view(t, 1:50), view(y, 1:50), view(err, 1:50), freqs)
            @test length(power(result)) == 200
        end
    end

    @testset "Stacked BGLS" begin
        @testset "Output dimensions" begin
            t, y, err, _ = make_single_sinusoid()
            freqs = collect(range(0.01, 0.5, length=50))
            result = sbgls(t, y, err, freqs; n_min=10)
            n_rows = length(t) - 10 + 1
            @test size(power(result)) == (n_rows, 50)
            @test length(obs_times(result)) == n_rows
            @test length(frequencies(result)) == 50
        end

        @testset "Row normalization to [0, 1]" begin
            t, y, err, _ = make_single_sinusoid()
            freqs = collect(range(0.01, 0.5, length=50))
            result = sbgls(t, y, err, freqs; n_min=10)
            mat = power(result)
            for row in eachrow(mat)
                @test minimum(row) >= -1e-12
                @test maximum(row) <= 1.0 + 1e-12
            end
        end

        @testset "Signal grows with data" begin
            rng = MersenneTwister(42)
            n_obs = 80; true_freq = 0.1; amplitude = 30.0; noise_level = 1.0
            t = sort(rand(rng, n_obs) .* 200.0)
            y = amplitude .* sin.(2π .* true_freq .* t) .+ randn(rng, n_obs) .* noise_level
            err = fill(noise_level, n_obs)

            freqs = collect(range(0.05, 0.15, length=50))
            result = sbgls(t, y, err, freqs; n_min=15)
            mat = power(result)

            # Find the column closest to true_freq
            _, col_idx = findmin(abs.(frequencies(result) .- true_freq))

            # The peak column should have high values in later rows
            last_quarter = mat[end-div(size(mat, 1), 4):end, col_idx]
            @test mean(last_quarter) > 0.5
        end

        @testset "n_min validation" begin
            t = collect(1.0:5.0)
            y = zeros(5)
            err = ones(5)
            freqs = [0.1]
            @test_throws ArgumentError sbgls(t, y, err, freqs; n_min=10)
        end

        @testset "Destructuring" begin
            t, y, err, _ = make_single_sinusoid()
            freqs = collect(range(0.01, 0.5, length=30))
            f, ot, mat = sbgls(t, y, err, freqs; n_min=10)
            @test length(f) == 30
            @test length(ot) == length(t) - 10 + 1
            @test size(mat, 2) == 30
        end

        @testset "sbgls_auto max_frequencies cap" begin
            t, y, err, _ = make_single_sinusoid()
            result = sbgls_auto(t, y, err; max_frequencies=50, n_min=10)
            @test length(frequencies(result)) <= 50
        end
    end
end
