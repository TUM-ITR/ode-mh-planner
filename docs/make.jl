using Documenter, OdeMHPlanner

makedocs(
    sitename="OdeMHPlanner.jl",
    modules=[OdeMHPlanner],
    checkdocs=:exports,
    warnonly=[:missing_docs],
    pages=[
        "Home" => "index.md",
        "Experiments" => [
            "Overview" => "experiments/overview.md",
            "Inference and Sampler Tuning" => "experiments/sampling.md",
            "Optimal Control" => "experiments/control.md",
            "Monte Carlo Study" => "experiments/monte_carlo.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(
    repo="github.com/TUM-ITR/ode-mh-planner.git",
)