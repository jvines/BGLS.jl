using Documenter
using BayesianGLS

makedocs(
    sitename = "BayesianGLS.jl",
    modules = [BayesianGLS],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://jvines.github.io/BayesianGLS.jl",
    ),
)

deploydocs(
    repo = "github.com/jvines/BayesianGLS.jl.git",
    devbranch = "main",
)
