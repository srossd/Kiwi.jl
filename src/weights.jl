"""
Weight vector in the weight space of a Lie algebra.
Weights are represented as linear combinations of fundamental weights: w = Σ c_i ω_i.
"""
struct Weight
    algebra::LieAlgebra
    coordinates::Vector{Rational{Int}}  # Coefficients in fundamental weight basis
    
    function Weight(algebra::LieAlgebra, coords::Vector{<:Real})
        new(algebra, Rational{Int}.(coords))
    end
end

# Arithmetic operations
Base.:+(w1::Weight, w2::Weight) = begin
    @assert w1.algebra == w2.algebra "Weights must be from the same algebra"
    Weight(w1.algebra, w1.coordinates .+ w2.coordinates)
end

Base.:-(w1::Weight, w2::Weight) = begin
    @assert w1.algebra == w2.algebra "Weights must be from the same algebra"
    Weight(w1.algebra, w1.coordinates .- w2.coordinates)
end

Base.:*(scalar::Real, w::Weight) = Weight(w.algebra, scalar .* w.coordinates)
Base.:*(w::Weight, scalar::Real) = scalar * w

Base.:(==)(w1::Weight, w2::Weight) = w1.algebra == w2.algebra && w1.coordinates == w2.coordinates
Base.hash(w::Weight, h::UInt) = hash((w.algebra, w.coordinates), h)

"""
    simple_root_squared_lengths(g::LieAlgebra)

Compute the squared lengths of simple roots from the Cartan matrix.
For a Cartan matrix C where C_ij = 2(α_i, α_j)/(α_j, α_j),
we can deduce the relative lengths of roots.
Normalized so that shortest root has squared length 1 and longest has squared length 2.
"""
function simple_root_squared_lengths(g::LieAlgebra)
    C = cartan_matrix(g)
    n = g.rank
    
    # For simply-laced algebras (A, D, E), all roots have the same length
    # For non-simply-laced, we determine lengths from C_ij * C_ji
    # C_ij * C_ji = 4(α_i,α_j)²/((α_i,α_i)(α_j,α_j))
    
    lengths_sq = ones(Rational{Int}, n)
    
    # Use the constraint that C_ij / C_ji determines the ratio of squared lengths
    # C_ij / C_ji = (α_j,α_j) / (α_i,α_i)
    
    for i in 1:n
        for j in i+1:n
            if C[i,j] != 0
                # (α_j,α_j)/(α_i,α_i) = C_ji/C_ij
                ratio = Rational{Int}(C[j,i]) // Rational{Int}(C[i,j])
                lengths_sq[j] = ratio * lengths_sq[i]
            end
        end
    end
    
    # Then scale so longest root has length² = 2  
    max_length = maximum(lengths_sq)
    lengths_sq = 2 * lengths_sq / max_length
    
    return lengths_sq
end

# Inner product
"""
    inner_product(w1::Weight, w2::Weight)

Compute the inner product of two weights using the Killing form.
For weights λ = Σ λ_i ω_i and μ = Σ μ_i ω_i in the fundamental weight basis,
with simple roots α_i = Σ_j C_ij ω_j (rows of Cartan matrix),
the metric is determined by the condition that (α_i, α_j) follows from the Cartan matrix.
"""
function inner_product(w1::Weight, w2::Weight)
    @assert w1.algebra == w2.algebra "Weights must be from the same algebra"
    
    g = w1.algebra
    C = cartan_matrix(g)
    n = length(w1.coordinates)
    
    # Get squared lengths of simple roots
    lengths_sq = simple_root_squared_lengths(g)
    
    # Build the metric tensor
    # We know: C_ij = 2(α_i, α_j)/(α_j, α_j)
    # So: (α_i, α_j) = C_ij * (α_j, α_j) / 2
    # And: α_i = Σ_k C_ik ω_k (rows of C)
    # So: (α_i, α_j) = Σ_{k,l} C_ik C_jl G_kl
    # This gives: C G C^T = target where target_ij = C_ij * (α_j, α_j) / 2
    
    # Build target matrix: T_ij = C_ij * (α_j, α_j) / 2
    D = diagm(lengths_sq)
    target = Rational{Int}.(C) * D / 2
    
    # Solve: C G C^T = target => G = C^{-1} target (C^T)^{-1}
    G = inv(Rational{Int}.(C)) * target * inv(Rational{Int}.(transpose(C)))
    
    result = Rational{Int}(0)
    for i in 1:n
        for j in 1:n
            result += w1.coordinates[i] * w2.coordinates[j] * G[i, j]
        end
    end
    return result
end

"""
    is_positive(w::Weight)

Check if a weight is positive (first non-zero coordinate is positive).
"""
function is_positive(w::Weight)
    for c in w.coordinates
        if c > 0
            return true
        elseif c < 0
            return false
        end
    end
    return false  # zero weight
end

"""
    fundamental_weights(g::LieAlgebra)

Return the fundamental weights of the Lie algebra.
In the fundamental weight basis, ω_i is simply the i-th standard basis vector.
"""
function fundamental_weights(g::LieAlgebra)
    n = g.rank
    weights = Vector{Weight}(undef, n)
    
    for i in 1:n
        coords = zeros(Rational{Int}, n)
        coords[i] = 1
        weights[i] = Weight(g, coords)
    end
    
    return weights
end

"""
    positive_roots(g::LieAlgebra)

Return all positive roots of the Lie algebra as Weight objects.
Follows the procedure in Georgi section 8.8: positive roots are built by 
systematically adding simple roots and checking the root string condition.
"""
function positive_roots(g::LieAlgebra)
    simple = simple_roots(g)
    n = length(simple)
    C = cartan_matrix(g)
    
    roots = copy(simple)
    seen = Set([r.coordinates for r in roots])
    processed = Set{Vector{Rational{Int}}}()
    
    # Build roots by adding simple roots
    # For each root α, check each component (α)_i = C[i,:]·α in the fundamental weight basis
    # If (α)_i < 0, then α is at the bottom of an α_i-string
    # and we can add α_i up to |(α)_i| times to get new roots
    queue = copy(roots)
    
    while !isempty(queue)
        α = popfirst!(queue)
        
        # Skip if already processed
        α.coordinates in processed && continue
        push!(processed, α.coordinates)
        
        for i in 1:n
            # Compute (α)_i = α_i · α where both are in fundamental weight basis
            # Since α_i = C[i,:] (i-th row of Cartan matrix), we have (α)_i = dot(C[i,:], α.coordinates)
            coeff_i = α.coordinates[i]
            
            # If coeff_i < 0, then α is at the bottom of an α_i-string
            # and we can add α_i up to |coeff_i| times
            if coeff_i < 0
                k_max = Int(-coeff_i)
                for k in 1:k_max
                    candidate = α + k * simple[i]
                    if !(candidate.coordinates in seen)
                        push!(roots, candidate)
                        push!(queue, candidate)
                        push!(seen, candidate.coordinates)
                    end
                end
            end
        end
    end
    
    return roots
end

"""
    simple_roots(g::LieAlgebra)

Return the simple roots of the Lie algebra in the fundamental weight basis.
The simple root α_i in fundamental weight coordinates: α_i = Σ_j C_ij ω_j (rows of C)
"""
function simple_roots(g::LieAlgebra)
    C = cartan_matrix(g)
    n = g.rank
    
    roots = Vector{Weight}(undef, n)
    for i in 1:n
        # α_i in fundamental weight basis: α_i = Σ_j C_ij ω_j (i-th row of C)
        coords = zeros(Rational{Int}, n)
        for j in 1:n
            coords[j] = C[i, j]
        end
        roots[i] = Weight(g, coords)
    end
    
    return roots
end

"""
    weyl_vector(g::LieAlgebra)

Compute the Weyl vector ρ as half the sum of positive roots.
"""
function weyl_vector(g::LieAlgebra)
    roots = positive_roots(g)
    if isempty(roots)
        return Weight(g, zeros(Rational{Int}, g.rank))
    end
    
    # Sum all positive roots
    sum_roots = roots[1]
    for i in 2:length(roots)
        sum_roots = sum_roots + roots[i]
    end
    
    # Return half the sum
    return (1//2) * sum_roots
end
