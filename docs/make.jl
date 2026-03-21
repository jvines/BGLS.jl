using Documenter
using BGLS

makedocs(
    sitename = "BGLS.jl",
    modules = [BGLS],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://jvines.github.io/BGLS.jl",
    ),
)

deploydocs(
    repo = "github.com/jvines/BGLS.jl.git",
    devbranch = "main",
)
