using Documenter
using Kiwi

makedocs(
    sitename = "Kiwi.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://srossd.github.io/Kiwi.jl",
        assets = String[],
        repolink = "https://github.com/srossd/Kiwi.jl",
    ),
    modules = [Kiwi],
    checkdocs = :none,  # Don't require all docstrings to be documented
    warnonly = [:cross_references, :missing_docs, :example_block, :autodocs_block],  # Allow warnings
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Getting Started" => "guide/getting_started.md",
            "Lie Algebras" => "guide/lie_algebras.md",
            "Representations" => "guide/representations.md",
            "Characters" => "guide/characters.md",
            "Tensor Products" => "guide/tensor_products.md",
            "Weyl Groups" => "guide/weyl_groups.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Lie Algebras" => "api/algebras.md",
            "Representations" => "api/representations.md",
            "Characters" => "api/characters.md",
            "Tensor Products" => "api/tensor_products.md",
            "Weyl Groups" => "api/weyl.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/srossd/Kiwi.jl.git",
    devbranch = "main",
    push_preview = true,
)
