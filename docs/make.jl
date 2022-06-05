push!(LOAD_PATH,"../src/")

using Documenter, Spike;

makedocs(sitename = "Spike.jl", pages = [
    "Home" => "index.md",
    "Manual" => [
        "Guide" => "guide.md",
        "Examples" => "examples.md"
    ],
    "Library" => "library.md"
])