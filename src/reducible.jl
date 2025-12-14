"""
Reducible representation type and tensor product decomposition using Klimyk's formula.
"""

"""
    Rep

Represents a (possibly reducible) representation as a linear combination of irreps.
Stored as a dictionary mapping irreps to their multiplicities.
"""
struct Rep
    algebra::LieAlgebra
    components::Dict{Irrep, Int}  # Irrep -> multiplicity
    
    function Rep(algebra::LieAlgebra, components::Dict{Irrep, Int})
        # Verify all irreps are from the same algebra
        for irrep in keys(components)
            if irrep.algebra != algebra
                error("All irreps must be from the same algebra")
            end
        end
        new(algebra, components)
    end
end

# Construct from single irrep
Rep(irrep::Irrep) = Rep(irrep.algebra, Dict(irrep => 1))

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
    
    components = copy(rep1.components)
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

    rep1 ⊗ rep2 = Σ_μ Σ_{w ∈ W} sign(w) · [w(λ₁ + μ)]

where:
- λ₁ is the highest weight of rep1 (we choose rep1 to be the smaller-dimensional one)
- The first sum is over all weights μ in the character of rep2
- W is the Weyl group
- [w(λ₁ + μ)] means: reflect w(λ₁ + μ) to the dominant chamber and interpret as an irrep
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
    
    result = Rep(rep1.algebra, Dict{Irrep, Int}())
    
    for (irrep1, mult1) in rep1.components
        for (irrep2, mult2) in rep2.components
            # Compute irrep1 ⊗ irrep2
            prod = tensor_product(irrep1, irrep2)
            
            # Add to result with combined multiplicity
            for (irrep, mult) in prod.components
                contribution = mult1 * mult2 * mult
                if haskey(result.components, irrep)
                    result.components[irrep] += contribution
                else
                    result.components[irrep] = contribution
                end
            end
        end
    end
    
    # Remove zero entries
    filter!(p -> p.second != 0, result.components)
    
    return result
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
1. Take the character of rep2 (all weights μ with their multiplicities)
2. For each weight μ, compute λ₁ + μ + ρ (where ρ is the Weyl vector)
3. Reflect to the dominant chamber: (λ₁ + μ + ρ)' with sign
4. Subtract ρ: ν = (λ₁ + μ + ρ)' - ρ
5. Check if ν is dominant - if yes, add [ν] with multiplicity = sign × mult_μ
6. Cancellations from negative signs give the correct decomposition
"""
function _tensor_product_klimyk(rep1::Irrep, rep2::Irrep)
    g = rep1.algebra
    λ₁ = highest_weight(rep1)
    ρ = weyl_vector(g)
    
    # Get the character of rep2 (all weights with multiplicities)
    char2 = character(rep2)
    
    # Dictionary to accumulate multiplicities (can be negative during computation)
    result_components = Dict{Irrep, Int}()
    
    # For each weight μ in the character of rep2
    for wm in char2
        μ = wm.weight
        mult_μ = wm.multiplicity
        
        # Compute λ₁ + μ + ρ
        shifted_weight = λ₁ + μ + ρ
        
        # Reflect to dominant chamber and get the sign
        reflected_weight, sign = reflect_to_dominant(shifted_weight)
        
        # Subtract ρ to get the final weight
        final_weight = reflected_weight - ρ
        
        # Check if final_weight is dominant
        if is_dominant(final_weight)
            # Check if the result has non-negative integer coordinates
            if all(isinteger.(final_weight.coordinates)) && all(final_weight.coordinates .>= 0)
                # The weight coordinates in fundamental weight basis ARE the Dynkin labels
                labels = [Int(c) for c in final_weight.coordinates]
                irrep = Irrep(g, labels)
                
                # Add the signed multiplicity from the character of rep2
                contribution = sign * mult_μ
                if haskey(result_components, irrep)
                    result_components[irrep] += contribution
                else
                    result_components[irrep] = contribution
                end
            end
        end
    end
    
    # Remove zero entries (from cancellations)
    filter!(p -> p.second != 0, result_components)
    
    return Rep(g, result_components)
end

# Unicode aliases (defined after tensor_product)
const ⊕ = +
const ⊗ = tensor_product
