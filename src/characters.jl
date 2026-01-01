"""
Character computation for irreducible representations using Freudenthal's formula.

The character of a representation is the collection of all weights with their multiplicities.
Freudenthal's recursive formula computes weight multiplicities by relating each weight
to higher weights in the representation.
"""

"""
    WeightMultiplicity

A weight with its multiplicity in a representation.
"""
struct WeightMultiplicity
    weight::Weight
    multiplicity::Int
end

# Display
function Base.show(io::IO, wm::WeightMultiplicity)
    print(io, "Weight(", wm.weight.coordinates, ") with multiplicity ", wm.multiplicity)
end

"""
    Character

Represents a character of an irreducible representation.
Can be eagerly computed (all weights precomputed) or lazy (weights computed on demand).

# Fields
- `algebra::LieAlgebra`: The Lie algebra
- `irrep::Union{Irrep, Nothing}`: The irrep (needed for lazy evaluation)
- `weights::Dict{Weight, Int}`: Computed weight multiplicities
- `lazy::Bool`: Whether to compute weights on demand
- `ρ::Union{Weight, Nothing}`: Weyl vector (cached for lazy evaluation)
- `positive_roots::Union{Vector{Weight}, Nothing}`: Positive roots (cached for lazy evaluation)
"""
mutable struct Character
    algebra::LieAlgebra
    irrep::Union{Irrep, Nothing}
    weights::Dict{Weight, Int}
    lazy::Bool
    ρ::Union{Weight, Nothing}
    positive_roots::Union{Vector{Weight}, Nothing}
    
    # Constructor for eager evaluation (backward compatibility)
    function Character(g::LieAlgebra, weights::Dict{Weight, Int})
        new(g, nothing, weights, false, nothing, nothing)
    end
    
    # Constructor for lazy or eager evaluation
    function Character(irrep::Irrep; lazy::Bool=false)
        g = irrep.algebra
        λ = highest_weight(irrep)
        weights = Dict{Weight, Int}(λ => 1)
        
        if lazy
            ρ = weyl_vector(g)
            positive_roots_list = positive_roots(g)
            new(g, irrep, weights, true, ρ, positive_roots_list)
        else
            new(g, irrep, weights, false, nothing, nothing)
        end
    end
end

# Display
function Base.show(io::IO, char::Character)
    if char.lazy
        n_computed = length(char.weights)
        print(io, "Character for ", char.irrep, " (", n_computed, " weights computed)")
    else
        n_weights = length(char.weights)
        dim = sum(values(char.weights))
        println(io, "Character with ", n_weights, " distinct weights (dimension ", dim, ")")
        
        # Sort weights for consistent display
        sorted_weights = sort(collect(char.weights), by = x -> x[1].coordinates, lt = lexicographic_compare)
        
        for (weight, mult) in sorted_weights
            println(io, "  Weight(", weight.coordinates, ") → multiplicity ", mult)
        end
    end
end

# Indexing: character[weight] returns multiplicity
function Base.getindex(char::Character, w::Weight)
    # Check if already computed
    if haskey(char.weights, w)
        return char.weights[w]
    end
    
    # If not lazy, weight doesn't exist
    if !char.lazy
        return 0
    end
    
    # Lazy evaluation: compute on demand
    w_dom, _ = reflect_to_dominant(w)
    
    # Compute multiplicity using Freudenthal's formula
    λ = highest_weight(char.irrep)
    mult = compute_multiplicity_freudenthal_lazy(w_dom, λ, char.ρ, char.positive_roots, char.weights)
    
    # Cache the result for entire Weyl orbit
    # All weights in an orbit have the same multiplicity
    orbit_dict = weyl_orbit(w_dom)
    orbit_weights = vcat(orbit_dict[1], orbit_dict[-1])
    for w_orbit in orbit_weights
        char.weights[w_orbit] = mult
    end
    
    return mult
end

# Iterate over weight-multiplicity pairs
function Base.iterate(char::Character, state=1)
    weights_vec = collect(char.weights)
    if state > length(weights_vec)
        return nothing
    end
    weight, mult = weights_vec[state]
    wm = WeightMultiplicity(weight, mult)
    return (wm, state + 1)
end

Base.length(char::Character) = length(char.weights)

"""
    character(rep::Irrep; lazy=false)

Compute the complete character of an irreducible representation.
Returns a Character object containing all weights with their multiplicities.

Uses Freudenthal's multiplicity formula, which recursively computes the multiplicity
of each weight μ from the multiplicities of weights μ + kα for positive roots α.

# Arguments
- `rep::Irrep`: The irreducible representation
- `lazy::Bool=false`: If true, computes weights on demand.
                      If false, all weights precomputed.

# Algorithm
Starting from the highest weight λ (with multiplicity 1), we work our way down
by subtracting positive roots. For each new weight μ, Freudenthal's formula gives:

    2(λ+ρ, μ) m(μ) = Σ_{α>0} Σ_{k≥1} 2(μ+kα, α) m(μ+kα)

where the sum is over all positive roots α and all k such that μ+kα is a weight.

# Returns
A `Character` object.

# Examples
```julia
# Eager evaluation (default): precomputes all weights
rep = Irrep(SU(3), [1, 0])
char = character(rep)

# Lazy evaluation: computes weights on demand
char_lazy = character(rep; lazy=true)
```

# References
- Humphreys, "Introduction to Lie Algebras and Representation Theory" (1972), §22.3
- Fulton & Harris, "Representation Theory" (1991), §15.2
"""
function character(rep::Irrep; lazy::Bool=false)
    # Create Character object (lazy or eager)
    char = Character(rep; lazy=lazy)
    
    # If not lazy, compute all weights now
    if !lazy
        g = rep.algebra
        λ = highest_weight(rep)
        ρ = weyl_vector(g)
        positive_roots_list = positive_roots(g)
        simple_roots_list = simple_roots(g)
        
        # Use the weights dict from char
        multiplicities = char.weights
        
        # Queue of weights to process, starting with highest weight
        # We process weights in order from "highest" to "lowest"
        to_process = [λ]
        processed = Set{Weight}()
        
        while !isempty(to_process)
            μ = popfirst!(to_process)
            
            # Skip if already processed
            if μ in processed
                continue
            end
            push!(processed, μ)
            
            # Generate new weights by subtracting simple roots only
            for α in simple_roots_list
                ν = μ - α
                
                # Check if this weight should be included
                if !(ν in processed) && !(ν in to_process)
                    # Check if we already know its multiplicity (possibly from a Weyl-related weight)
                    if haskey(multiplicities, ν)
                        # We know the multiplicity, but still need to process it to explore further
                        if multiplicities[ν] > 0
                            push!(to_process, ν)
                        end
                    else
                        # Try to compute multiplicity using Freudenthal's formula
                        mult = compute_multiplicity_freudenthal(ν, λ, ρ, positive_roots_list, multiplicities)
                        
                        if mult > 0
                            # Cache for entire Weyl orbit to avoid redundant calculations
                            orbit_dict = weyl_orbit(ν)
                            for w_orbit in vcat(orbit_dict[1], orbit_dict[-1])
                                multiplicities[w_orbit] = mult
                            end
                            
                            push!(to_process, ν)
                        end
                    end
                end
            end
        end
    end
    
    return char
end

"""
    compute_multiplicity_freudenthal(μ, λ, ρ, roots, known_multiplicities)

Compute the multiplicity of weight μ in the irrep with highest weight λ
using Freudenthal's recursive formula.

Formula: (|λ+ρ|² - |μ+ρ|²) m(μ) = Σ_{α>0} Σ_{k≥1} 2(μ+kα, α) m(μ+kα)

Returns 0 if the weight is not in the representation.
"""
function compute_multiplicity_freudenthal(
    μ::Weight, 
    λ::Weight, 
    ρ::Weight,
    roots::Vector{Weight},
    known_multiplicities::Dict{Weight, Int}
)
    # Left-hand side coefficient: |λ+ρ|² - |μ+ρ|²
    λ_plus_ρ = λ + ρ
    μ_plus_ρ = μ + ρ
    lhs_coeff = inner_product(λ_plus_ρ, λ_plus_ρ) - inner_product(μ_plus_ρ, μ_plus_ρ)
    
    # If coefficient is zero or negative, this weight is not in the representation
    # (Actually, for weights in the rep, this should always be positive)
    if lhs_coeff <= 0
        return 0
    end
    
    # Right-hand side: sum over positive roots and multiplicities
    rhs = Rational{Int}(0)
    
    for α in roots
        # Sum over k ≥ 1 such that μ + kα is a known weight
        k = 1
        while true
            μ_plus_kα = μ + k * α
            
            # Check if this is a known weight
            if haskey(known_multiplicities, μ_plus_kα)
                mult = known_multiplicities[μ_plus_kα]
                
                # Add contribution: 2(μ+kα, α) * m(μ+kα)
                contribution = 2 * inner_product(μ_plus_kα, α) * mult
                rhs += contribution
                
                k += 1
            else
                # No more weights of the form μ + kα
                break
            end
        end
    end
    
    # Solve for m(μ): m(μ) = rhs / lhs_coeff
    if rhs == 0
        return 0
    end
    
    multiplicity = rhs / lhs_coeff
    
    # Multiplicity should be a non-negative integer
    if denominator(multiplicity) != 1 || numerator(multiplicity) < 0
        return 0
    end
    
    return Int(numerator(multiplicity))
end

"""
    lexicographic_compare(v1, v2)

Compare two vectors lexicographically.
Returns true if v1 < v2 in lexicographic order.
"""
function lexicographic_compare(v1::Vector, v2::Vector)
    for i in 1:min(length(v1), length(v2))
        if v1[i] < v2[i]
            return true
        elseif v1[i] > v2[i]
            return false
        end
    end
    return length(v1) < length(v2)
end

"""
    dimension(char::Character)

Compute the dimension of a character.
For eager characters, sums all weight multiplicities.
For lazy characters, uses the irrep's dimension formula.
"""
function dimension(char::Character)
    if char.lazy && !isnothing(char.irrep)
        return dimension(char.irrep)
    else
        return sum(values(char.weights))
    end
end

"""
    compute_multiplicity_freudenthal_lazy(μ, λ, ρ, roots, known_multiplicities)

Lazy version of Freudenthal's formula that recursively computes needed multiplicities.
This version will recursively compute any missing multiplicities as needed.
"""
function compute_multiplicity_freudenthal_lazy(
    μ::Weight,
    λ::Weight,
    ρ::Weight,
    roots::Vector{Weight},
    known_multiplicities::Dict{Weight, Int}
)
    # Check if already computed
    if haskey(known_multiplicities, μ)
        return known_multiplicities[μ]
    end
    
    # Left-hand side coefficient: |λ+ρ|² - |μ+ρ|²
    λ_plus_ρ = λ + ρ
    μ_plus_ρ = μ + ρ
    lhs_coeff = inner_product(λ_plus_ρ, λ_plus_ρ) - inner_product(μ_plus_ρ, μ_plus_ρ)
    
    # If coefficient is zero or negative, this weight is not in the representation
    if lhs_coeff <= 0
        known_multiplicities[μ] = 0
        return 0
    end
    
    # Right-hand side: sum over positive roots and multiplicities
    rhs = Rational{Int}(0)
    
    for α in roots
        # Sum over k ≥ 1 such that μ + kα is a weight
        k = 1
        while true
            μ_plus_kα = μ + k * α
            
            # Recursively compute if not known
            if !haskey(known_multiplicities, μ_plus_kα)
                mult = compute_multiplicity_freudenthal_lazy(μ_plus_kα, λ, ρ, roots, known_multiplicities)
                if mult == 0
                    break
                end
            end
            
            # Check if this is a known weight with non-zero multiplicity
            if haskey(known_multiplicities, μ_plus_kα)
                mult = known_multiplicities[μ_plus_kα]
                if mult == 0
                    break
                end
                
                # Add contribution: 2(μ+kα, α) * m(μ+kα)
                contribution = 2 * inner_product(μ_plus_kα, α) * mult
                rhs += contribution
                
                k += 1
            else
                break
            end
        end
    end
    
    # Solve for m(μ): m(μ) = rhs / lhs_coeff
    if rhs == 0
        known_multiplicities[μ] = 0
        return 0
    end
    
    multiplicity = rhs / lhs_coeff
    
    # Multiplicity should be a non-negative integer
    if denominator(multiplicity) != 1 || numerator(multiplicity) < 0
        known_multiplicities[μ] = 0
        return 0
    end
    
    result = Int(numerator(multiplicity))
    known_multiplicities[μ] = result
    return result
end

"""
    dominant_weights(rep::Irrep)

Compute the dominant weights of an irreducible representation.

Returns a dictionary mapping dominant weights (weights with all non-negative Dynkin labels)
to their multiplicities in the character.

# Arguments
- `rep`: The irreducible representation

# Returns
A `Dict{Weight, Int}` containing only the dominant weights and their multiplicities.

# Example
```julia
g = A_series(2)
rep = Irrep(g, [2, 1])
dom_weights = dominant_weights(rep)
```
"""
function dominant_weights(rep::Irrep)
    # Compute the full character
    char = character(rep)
    
    # Filter to only dominant weights
    result = Dict{Weight, Int}()
    for (weight, mult) in char.weights
        if is_dominant(weight)
            result[weight] = mult
        end
    end
    
    return result
end
