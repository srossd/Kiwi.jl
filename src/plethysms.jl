"""
Plethysm computations for Lie algebras.

This module implements the LiE software approach to computing plethysms,
which involves alternating sums and dominant weight projections.
"""

"""
    alternating_dominant(weights::Dict{Weight, Int}, w::WeylWord=longest_weyl_word(g))

Apply the alternating dominant transformation to a weight dictionary using a Weyl word.

This is a key step in plethysm computations. For each simple reflection k in the Weyl word
(read left to right), the weight dictionary is transformed as follows:
- For each weight λ with multiplicity n:
  - If λ_k ≥ 0 (the k-th Dynkin label is non-negative): leave it unchanged
  - If λ_k = -1: delete the weight
  - If λ_k ≤ -2: replace λ with λ - (λ_k + 1)α_k and flip the sign of n

where λ_k = 2⟨λ, α_k⟩/⟨α_k, α_k⟩ is the k-th Dynkin label and α_k is the k-th simple root.

If no Weyl word is provided, defaults to the longest element of the Weyl group.

# Arguments
- `weights`: Dictionary mapping weights to their multiplicities
- `w`: The Weyl word specifying the sequence of transformations (default: longest_weyl_word)

# Returns
A new `Dict{Weight, Int}` with the transformations applied.

# Example
```julia
g = A_series(2)
rep = Irrep(g, [1, 0])
char = Character(rep)
weights_alt = alternating_dominant(char.weights)  # Uses longest word by default
```
"""
function alternating_dominant(weights::Dict{Weight, Int}, w::Union{WeylWord, Nothing}=nothing)
    # Get algebra from the first weight (all weights should have the same algebra)
    if isempty(weights)
        return Dict{Weight, Int}()
    end
    
    g = first(keys(weights)).algebra
    
    # Default to longest Weyl word if not provided
    if w === nothing
        w = longest_weyl_word(g)
    end
    
    @assert w.algebra == g "Weights and WeylWord must be from the same algebra"
    
    roots = simple_roots(g)
    
    # Start with a copy of the weights
    result_weights = copy(weights)
    
    # Process each reflection in the word from left to right
    for k in w.word
        α_k = roots[k]
        new_weights = Dict{Weight, Int}()
        
        for (λ, n) in result_weights
            if λ[k] >= 0
                # Leave it unchanged
                new_weights[λ] = get(new_weights, λ, 0) + n
            elseif λ[k] == -1
                # Delete it (do nothing)
                continue
            else  # λ[k] <= -2
                # Replace λ with λ - (λ[k] + 1)α_k and flip sign of n
                # Create new weight: λ - (λ[k] + 1)α_k
                new_λ_coords = λ.coordinates - (λ[k] + 1) * α_k.coordinates
                new_λ = Weight(g, new_λ_coords)
                new_weights[new_λ] = get(new_weights, new_λ, 0) - n
            end
        end
        
        result_weights = new_weights
    end

    result_weights = Dict(filter(t -> t[2] != 0, result_weights))
    
    return result_weights
end

"""
    virtual_decomposition(weight::Weight)

Compute the virtual decomposition of a weight by taking its Weyl orbit and applying
the alternating dominant projection with the longest Weyl word.

The Weyl orbit includes alternating signs for reflected weights, and the alternating
dominant operation projects these to the dominant chamber, resulting in a virtual
character (a formal sum with integer coefficients, which may be negative).

# Arguments
- `weight`: The weight to decompose

# Returns
A `Dict{Weight, Int}` representing the virtual character as a sum of weights with multiplicities.

# Example
```julia
g = A_series(2)
λ = Weight(g, [2, 1])
virt = virtual_decomposition(λ)
```
"""
function virtual_decomposition(weight::Weight)
    # Get the Weyl orbit with alternating signs
    orbit = weyl_orbit(weight)
    
    # Convert to a weight dictionary
    weights = Dict{Weight, Int}()
    for (_, weight_list) in orbit
        for w in weight_list
            weights[w] = 1
        end
    end
    
    # Apply alternating_dominant with longest Weyl word (default)
    return alternating_dominant(weights)
end

"""
    adams(n::Int, irrep::Irrep)

Compute the n-th Adams operation on an irreducible representation.

This operation rescales all dominant weights by n and computes the virtual decomposition
of each rescaled weight, weighted by its multiplicity.

# Arguments
- `n`: The scaling factor (positive integer)
- `irrep`: The irreducible representation

# Returns
A `Dict{Weight, Int}` representing the result as a virtual character.

# Example
```julia
g = A_series(2)
rep = Irrep(g, [1, 0])
result = adams(2, rep)
```
"""
function adams(n::Int, irrep::Irrep)
    g = irrep.algebra
    ad = algebra_data(g)
    W, M = _full_character(g, irrep.dynkin_labels)
    # ψⁿ(χ) = Σ_μ m_μ e^{nμ}; decompose with the Brauer/Racah–Speiser rule
    out = Dict{Vector{Int}, Int}()
    _klimyk!(out, ad, ones(Int, ad.n), 1, W, M; scale = n)
    result = Dict{Weight, Int}()
    for (λ, c) in out
        c == 0 && continue
        result[Weight(g, λ)] = c
    end
    return result
end

"""
    plethysm(lambda::Irrep, rho::SymmetricIrrep)

Compute the plethysm of a Lie algebra irrep λ with a symmetric group irrep ρ.

This computes the character of the composition λ[ρ], which describes how the 
representation λ transforms under the symmetric group action.

# Algorithm
For each conjugacy class of S_n with partition (p₁, ..., pₖ):
1. Compute adams(pᵢ, λ) for each part pᵢ
2. Take the tensor product of all these (treating as reducible representations)
3. Weight by (conjugacy class size) × character(ρ, class) / n!
4. Sum over all conjugacy classes

# Arguments
- `lambda`: The Lie algebra irreducible representation
- `rho`: The symmetric group irreducible representation

# Returns
A `Rep` object representing the decomposition into Lie algebra irreps.

# Example
```julia
g = A_series(2)
lambda = Irrep(g, [1, 0])  # Fundamental rep of SU(3)
rho = SymmetricIrrep([2])   # Symmetric square representation of S_2
result = plethysm(lambda, rho)
```
"""
function plethysm(lambda::Irrep, rho::SymmetricIrrep)
    return plethysm(Rep(lambda), rho)
end

"""
    plethysm(V::Rep, rho::SymmetricIrrep)

Plethysm of a (possibly reducible) representation `V` with the Schur functor
labelled by the partition of `rho`.
"""
function plethysm(V::Rep, rho::SymmetricIrrep)
    g = V.algebra
    ad = algebra_data(g)
    parts = rho.partition.parts
    n = isempty(parts) ? 0 : sum(parts)
    if n == 0
        return Rep(g, Dict(Irrep(g, zeros(Int, g.rank)) => 1))
    end
    W, M = _rep_full_character(V)
    d = isempty(M) ? 0 : sum(big, M)
    T = _int_type_for(big(d + 1)^n * factorial(big(n)) * 4)
    memo = Dict{Vector{Int}, Dict{Vector{Int}, T}}()
    memo[Int[]] = Dict{Vector{Int}, T}(zeros(Int, ad.n) => one(T))
    MT = T.(M)
    result = _schur_plethysm!(memo, ad, copy(parts), W, MT)
    return _dict_to_rep(g, result)
end

"""
    _schur_plethysm!(memo, ad, λ, W, M)

Decomposition of the Schur functor s_λ applied to the representation with
weights `W` and multiplicities `M`.  Uses the identity
`Σ_k p_k p_k^⊥ = deg` together with the Murnaghan–Nakayama rule
`p_k^⊥ s_λ = Σ_{λ/μ border strip of size k} (-1)^{ht} s_μ`, which gives

    n · s_λ[V] = Σ_{k=1}^{n} ψᵏ(V) ⊗ Σ_{λ/μ = k-strip} (-1)^{ht} s_μ[V].

Each step is a Klimyk product of a genuine representation with an Adams
operation of `V`; the results for smaller partitions are memoised.
"""
function _schur_plethysm!(memo::Dict{Vector{Int}, Dict{Vector{Int}, T}}, ad::AlgebraData,
                          λ::Vector{Int}, W::Vector{Vector{Int}}, M::Vector{T}) where {T}
    r = get(memo, λ, nothing)
    r === nothing || return r
    n = sum(λ)
    total = Dict{Vector{Int}, T}()
    for k in 1:n
        Y = Dict{Vector{Int}, T}()
        for (μ, sgn) in _rim_hook_removals(λ, k)
            for (ν, c) in _schur_plethysm!(memo, ad, μ, W, M)
                Y[ν] = get(Y, ν, zero(T)) + sgn * c
            end
        end
        for (ν, c) in Y
            c == 0 && continue
            _klimyk!(total, ad, ν .+ 1, c, W, M; scale = k)
        end
    end
    result = Dict{Vector{Int}, T}()
    for (ν, c) in total
        c == 0 && continue
        q, rem = divrem(c, n)
        rem == 0 || error("Non-integer multiplicity in plethysm (internal error)")
        result[ν] = q
    end
    memo[copy(λ)] = result
    return result
end
"""
    symmetric_power(n::Int, irrep::Irrep)

Compute the n-th symmetric power of an irreducible representation.

This is a convenience function that computes the plethysm with the trivial
(all-symmetric) representation of Sₙ.

# Arguments
- `n`: The power (positive integer)
- `irrep`: The irreducible representation

# Returns
A `Rep` object representing the decomposition into irreps.

# Example
```julia
g = SU(3)
fund = Irrep(g, [1, 0])
sym2 = symmetric_power(2, fund)  # S²(3) = [2,0]
```
"""
function symmetric_power(n::Int, irrep::Irrep)
    # The n-th symmetric power corresponds to the partition [n]
    rho = SymmetricIrrep([n])
    return plethysm(irrep, rho)
end

"""
    antisymmetric_power(n::Int, irrep::Irrep)

Compute the n-th antisymmetric (exterior) power of an irreducible representation.

This is a convenience function that computes the plethysm with the alternating
representation of Sₙ.

# Arguments
- `n`: The power (positive integer)
- `irrep`: The irreducible representation

# Returns
A `Rep` object representing the decomposition into irreps.

# Example
```julia
g = SU(3)
fund = Irrep(g, [1, 0])
alt2 = antisymmetric_power(2, fund)  # Λ²(3) = [0,1]
```
"""
function antisymmetric_power(n::Int, irrep::Irrep)
    # The n-th antisymmetric power corresponds to the partition [1,1,...,1] (n times)
    rho = SymmetricIrrep(fill(1, n))
    return plethysm(irrep, rho)
end

