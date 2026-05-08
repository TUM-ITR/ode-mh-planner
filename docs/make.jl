using Documenter, OdeMHPlanner

makedocs(
    sitename="OdeMHPlanner.jl",
    modules=[OdeMHPlanner],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Inference and Sampler Tuning" => "examples/sampling.md",
            "Optimal Control" => "examples/control.md",
        ],
        "Experiments" => "experiments.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo="github.com/TUM-ITR/ode-mh-planner.git",
)