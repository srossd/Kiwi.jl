using Documenter
using Kiwi

makedocs(
    sitename = "Kiwi.jl",
    format = Documenter.HTML(
        prettyurls = false,
        assets = String[],
    ),
    modules = [Kiwi],
    remotes = nothing,  # Disable remote source links for local build
    checkdocs = :none,  # Don't require all docstrings to be documented
    warnonly = [:cross_references, :missing_docs, :example_block],  # Allow warnings
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

# Remove deploydocs for local building
# deploydocs(
#     repo = "github.com/srossd/Kiwi.jl.git",
#     devbranch = "main",
# )
