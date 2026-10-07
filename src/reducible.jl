"""
Reducible representation type and tensor product decomposition using Klimyk's formula.
"""

"""
    Rep

Represents a (possibly reducible) representation as a linear combination of irreps.
Stored as a dictionary mapping irreps to their multiplicities.

The multiplicity type `T` is `Int` unless a computation produces multiplicities
that do not fit in 64 bits, in which case it is `BigInt`.
"""
struct Rep{T<:Integer}
    algebra::LieAlgebra
    components::Dict{Irrep, T}  # Irrep -> multiplicity

    function Rep(algebra::LieAlgebra, components::Dict{Irrep, T}) where {T<:Integer}
        # Verify all irreps are from the same algebra
        for irrep in keys(components)
            if irrep.algebra != algebra
                error("All irreps must be from the same algebra")
            end
        end
        new{T}(algebra, components)
    end
end

# Construct from single irrep
Rep(irrep::Irrep) = Rep(irrep.algebra, Dict{Irrep, Int}(irrep => 1))

# Construct from pairs of irreps and multiplicities
function Rep(algebra::LieAlgebra, pairs::Pair{Irrep, Int}...)
    Rep(algebra, Dict(pairs...))
end

# Display
function Base.show(io::IO, rep::Rep)
    if isempty(rep.components)
        print(io, "0")
        return
    end
    
    terms = String[]
    for (irrep, mult) in sort(collect(rep.components), by=x->x[1].dynkin_labels)
        label_str = "[" * join(irrep.dynkin_labels, ",") * "]"
        if mult == 1
            push!(terms, label_str)
        else
            push!(terms, "$mult × $label_str")
        end
    end
    print(io, join(terms, " ⊕ "))
end

# Addition of representations
function Base.:+(rep1::Rep, rep2::Rep)
    @assert rep1.algebra == rep2.algebra "Representations must be from the same algebra"
    
    T = promote_type(valtype(rep1.components), valtype(rep2.components))
    components = Dict{Irrep, T}(rep1.components)
    for (irrep, mult) in rep2.components
        if haskey(components, irrep)
            components[irrep] += mult
        else
            components[irrep] = mult
        end
    end
    
    return Rep(rep1.algebra, components)
end

# Scalar multiplication
function Base.:*(n::Int, rep::Rep)
    components = Dict(irrep => n * mult for (irrep, mult) in rep.components)
    return Rep(rep.algebra, components)
end

Base.:*(rep::Rep, n::Int) = n * rep

# Get irreps with their multiplicities
function irreps(rep::Rep)
    return collect(rep.components)
end

# Indexing: rep[irrep] returns multiplicity
function Base.getindex(rep::Rep, irrep::Irrep)
    return get(rep.components, irrep, 0)
end

"""
    dimension(rep::Rep)

Compute the dimension of a (possibly reducible) representation.
This is the sum of dimensions of all component irreps, weighted by their multiplicities.
"""
function dimension(rep::Rep)
    return sum(dimension(irrep) * mult for (irrep, mult) in rep.components)
end

"""
    tensor_product(rep1::Irrep, rep2::Irrep)

Compute the tensor product of two irreducible representations using Klimyk's formula.

Klimyk's formula (also known as the Racah-Speiser algorithm) states that:

    rep1 ⊗ rep2 = Σ_μ sign(w) · [w(λ + μ + ρ) - ρ]

where:
- λ is the highest weight of the larger-dimensional irrep
- The sum is over all weights μ (with multiplicity) of the smaller irrep, whose
  character is obtained from its dominant weights (Freudenthal) and Weyl orbits
- w is the Weyl group element taking λ + μ + ρ to the dominant chamber; terms on a
  chamber wall are dropped
- sign(w) is the determinant of the Weyl group element w

The formula automatically handles the necessary cancellations, and only dominant 
weights contribute to the final result.

# References
- Klimyk, "Decomposition of the direct product of irreducible representations" (1968)
- Racah, "Group Theory and Spectroscopy" (1965)

# Example
```julia
g = A_series(2)
fund1 = Irrep(g, 1, 0)
fund2 = Irrep(g, 0, 1)
result = tensor_product(fund1, fund2)
```
"""
function tensor_product(rep1::Irrep, rep2::Irrep)
    @assert rep1.algebra == rep2.algebra "Representations must be from the same algebra"

    # Choose the order so we compute the character of the smaller representation
    dim1 = dimension(rep1)
    dim2 = dimension(rep2)

    if dim1 <= dim2
        return _tensor_product_klimyk(rep1, rep2)
    else
        return _tensor_product_klimyk(rep2, rep1)
    end
end

"""
    tensor_product(rep1::Rep, rep2::Rep)

Compute the tensor product of two (possibly reducible) representations.

Uses bilinearity: (⊕ᵢ nᵢRᵢ) ⊗ (⊕ⱼ mⱼSⱼ) = ⊕ᵢⱼ (nᵢmⱼ)(Rᵢ ⊗ Sⱼ)
"""
function tensor_product(rep1::Rep, rep2::Rep)
    @assert rep1.algebra == rep2.algebra "Representations must be from the same algebra"
    g = rep1.algebra
    ad = algebra_data(g)

    # Expand the character of the smaller side once, then apply Klimyk's formula
    # to every irrep of the other side.
    if dimension(rep1) < dimension(rep2)
        rep1, rep2 = rep2, rep1
    end
    W, M = _rep_full_character(rep2)
    T = promote_type(valtype(rep1.components), valtype(rep2.components))
    X = Dict{Vector{Int}, T}(irrep.dynkin_labels => mult for (irrep, mult) in rep1.components)
    out = _virtual_times_character(ad, X, W, M)
    return _dict_to_rep(g, out)
end

# All weights (merged) of a possibly reducible / virtual representation
function _rep_full_character(rep::Rep)
    g = rep.algebra
    W = Vector{Vector{Int}}()
    M = Vector{valtype(rep.components)}()
    for (irrep, mult) in rep.components
        mult == 0 && continue
        w, m = _full_character(g, irrep.dynkin_labels)
        append!(W, w)
        append!(M, m .* mult)
    end
    length(rep.components) > 1 ? _merge_weights(W, M) : (W, M)
end

function _dict_to_rep(g::LieAlgebra, d::Dict{Vector{Int}, T}) where {T}
    S = all(c -> typemin(Int) <= c <= typemax(Int), values(d)) ? Int : BigInt
    comps = Dict{Irrep, S}()
    for (λ, c) in d
        c == 0 && continue
        comps[Irrep(g, λ)] = S(c)
    end
    return Rep(g, comps)
end

"""
    tensor_product(rep1::Irrep, rep2::Rep)
    tensor_product(rep1::Rep, rep2::Irrep)

Tensor product between irreducible and reducible representations.
"""
tensor_product(rep1::Irrep, rep2::Rep) = tensor_product(Rep(rep1), rep2)
tensor_product(rep1::Rep, rep2::Irrep) = tensor_product(rep1, Rep(rep2))

"""
Internal implementation of Klimyk's formula.

Algorithm:
1. Take the character of rep1, the smaller irrep (all weights μ with multiplicities)
2. For each weight μ, compute λ₂ + μ + ρ (λ₂ = highest weight of rep2, ρ = Weyl vector)
3. Reflect to the dominant chamber: (λ₂ + μ + ρ)' with sign; drop it if it lies on a wall
4. Subtract ρ: ν = (λ₂ + μ + ρ)' - ρ and add [ν] with multiplicity sign × mult_μ
6. Cancellations from negative signs give the correct decomposition
"""
function _tensor_product_klimyk(rep1::Irrep, rep2::Irrep)
    g = rep1.algebra
    ad = algebra_data(g)
    W, M = _full_character(g, rep1.dynkin_labels)
    out = Dict{Vector{Int}, Int}()
    _klimyk!(out, ad, rep2.dynkin_labels .+ 1, 1, W, M)
    filter!(p -> p.second != 0, out)
    return _dict_to_rep(g, out)
end

# Unicode aliases (defined after tensor_product)
const ⊕ = +
const ⊗ = tensor_product
