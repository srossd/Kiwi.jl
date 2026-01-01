"""
    Irrep

Represents an irreducible representation of a Lie algebra by its highest weight.
The highest weight is specified by Dynkin labels.
"""
struct Irrep
    algebra::LieAlgebra
    dynkin_labels::Vector{Int}
    
    function Irrep(algebra::LieAlgebra, labels::Vector{Int})
        if length(labels) != lie_rank(algebra)
            error("Number of Dynkin labels ($(length(labels))) must equal rank of algebra ($(lie_rank(algebra)))")
        end
        if any(labels .< 0)
            error("Dynkin labels must be non-negative")
        end
        new(algebra, labels)
    end
end

# Convenient constructor with varargs
Irrep(algebra::LieAlgebra, labels::Int...) = Irrep(algebra, collect(labels))

# Display
function Base.show(io::IO, rep::Irrep)
    print(io, "Irrep($(rep.algebra), [$(join(rep.dynkin_labels, ", "))])")
end

"""
    dynkin_labels(rep::Irrep)

Get the Dynkin labels of an irreducible representation.
"""
dynkin_labels(rep::Irrep) = rep.dynkin_labels

"""
    highest_weight(rep::Irrep)

Convert Dynkin labels to highest weight as a Weight object.
In the fundamental weight basis, the highest weight coordinates are just the Dynkin labels.
"""
function highest_weight(rep::Irrep)
    # In fundamental weight basis, highest weight is λ = Σ λ_i ω_i
    # So the coordinates are just the Dynkin labels
    return Weight(rep.algebra, rep.dynkin_labels)
end

"""
    is_trivial(rep::Irrep)

Check if the representation is the trivial representation.
"""
is_trivial(rep::Irrep) = all(rep.dynkin_labels .== 0)

"""
    conjugate(rep::Irrep)

Return the conjugate (dual) representation.
"""
function conjugate(rep::Irrep)
    g = rep.algebra
    
    if g.series == :A
        # For A_n, conjugate reverses Dynkin labels
        return Irrep(g, reverse(rep.dynkin_labels))
    elseif g.series == :B || g.series == :C || g.series == :F || g.series == :G
        # For B_n, C_n, F_4, and G_2, representations are self-conjugate
        return rep
    elseif g.series == :D
        # For D_n, swap last two Dynkin labels
        labels = copy(rep.dynkin_labels)
        n = length(labels)
        labels[n-1], labels[n] = labels[n], labels[n-1]
        return Irrep(g, labels)
    elseif g.series == :E
        # For E_6, conjugate swaps nodes 1↔6, 3↔5
        # For E_7 and E_8, representations are self-conjugate
        if g.rank == 6
            labels = copy(rep.dynkin_labels)
            labels[1], labels[6] = labels[6], labels[1]
            labels[3], labels[5] = labels[5], labels[3]
            return Irrep(g, labels)
        else
            return rep
        end
    end
    
    error("Conjugate not implemented for this algebra type")
end

"""
    is_self_conjugate(rep::Irrep)

Check if the representation is self-conjugate.
"""
is_self_conjugate(rep::Irrep) = conjugate(rep).dynkin_labels == rep.dynkin_labels

# Equality
Base.:(==)(rep1::Irrep, rep2::Irrep) = rep1.algebra == rep2.algebra && rep1.dynkin_labels == rep2.dynkin_labels

# Hash for using as dictionary keys
Base.hash(rep::Irrep, h::UInt) = hash((rep.algebra, rep.dynkin_labels), h)

"""
    CongruencyClass

Represents the congruency class (character of the center) of an irreducible representation.
Different Lie algebras have different center structures:
- A_n: Z_{n+1}
- B_n: Z_2
- C_n: Z_2
- D_{2n}: Z_2 × Z_2
- D_{2n+1}: Z_4
- E_6: Z_3
- E_7: Z_2
- E_8, F_4, G_2: trivial
"""
struct CongruencyClass
    algebra::LieAlgebra
    value::Union{Int, Tuple{Int,Int}, Nothing}
    
    function CongruencyClass(algebra::LieAlgebra, value::Union{Int, Tuple{Int,Int}, Nothing})
        # Validate the value based on algebra type
        if algebra.series == :A
            # Z_{n+1}
            if !(value isa Int)
                error("A_n congruency class must be an integer")
            end
        elseif algebra.series == :B || algebra.series == :C || algebra.series == :E && algebra.rank == 7
            # Z_2
            if !(value isa Int)
                error("Congruency class must be an integer")
            end
        elseif algebra.series == :D
            if iseven(algebra.rank)
                # D_{2n}: Z_2 × Z_2
                if !(value isa Tuple{Int,Int})
                    error("D_{2n} congruency class must be a pair of integers")
                end
            else
                # D_{2n+1}: Z_4
                if !(value isa Int)
                    error("D_{2n+1} congruency class must be an integer")
                end
            end
        elseif algebra.series == :E && algebra.rank == 6
            # Z_3
            if !(value isa Int)
                error("E_6 congruency class must be an integer")
            end
        elseif (algebra.series == :E && algebra.rank == 8) || algebra.series == :F || algebra.series == :G
            # Trivial center
            if value !== nothing
                error("E_8, F_4, and G_2 have trivial center")
            end
        end
        
        new(algebra, value)
    end
end

# Display
function Base.show(io::IO, cc::CongruencyClass)
    if cc.value === nothing
        print(io, "CongruencyClass(trivial)")
    elseif cc.value isa Tuple
        print(io, "CongruencyClass($(cc.value[1]), $(cc.value[2]))")
    else
        print(io, "CongruencyClass($(cc.value))")
    end
end

# Equality
Base.:(==)(cc1::CongruencyClass, cc2::CongruencyClass) = cc1.algebra == cc2.algebra && cc1.value == cc2.value

# Hash
Base.hash(cc::CongruencyClass, h::UInt) = hash((cc.algebra, cc.value), h)

"""
    congruency_class(rep::Irrep)

Compute the congruency class of an irreducible representation from its Dynkin labels.

The formulas for each algebra type are:
- A_n: p_1 + 2p_2 + ... + n*p_n (mod n+1)
- B_n: p_n (mod 2)
- C_n: p_1 + p_3 + p_5 + ... (mod 2)
- D_{2n}: (p_{2n-1} + p_{2n}, p_1 + p_3 + ... + p_{2n-3} + (n-1)*p_{2n-1} + n*p_{2n}) (mod 2 for each)
- D_{2n+1}: 2(p_1 + p_3 + ... + p_{2n-1}) + (2n-1)*p_{2n} + (2n+1)*p_{2n+1} (mod 4)
- E_6: p_1 - p_2 + p_4 - p_5 (mod 3)
- E_7: p_4 + p_6 + p_7 (mod 2)
- E_8, F_4, G_2: trivial
"""
function congruency_class(rep::Irrep)
    g = rep.algebra
    p = rep.dynkin_labels
    n = g.rank
    
    if g.series == :A
        # p_1 + 2*p_2 + ... + n*p_n (mod n+1)
        value = sum(i * p[i] for i in 1:n) % (n + 1)
        return CongruencyClass(g, value)
        
    elseif g.series == :B
        # p_n (mod 2)
        value = p[n] % 2
        return CongruencyClass(g, value)
        
    elseif g.series == :C
        # p_1 + p_3 + p_5 + ... (mod 2)
        value = sum(p[i] for i in 1:2:n) % 2
        return CongruencyClass(g, value)
        
    elseif g.series == :D
        if iseven(n)
            # D_{2n}: (p_{n-1} + p_n, p_1 + p_3 + ... + p_{n-3} + (n/2-1)*p_{n-1} + (n/2)*p_n)
            # First component
            val1 = (p[n-1] + p[n]) % 2
            # Second component
            half_n = div(n, 2)
            sum_odd = sum(p[i] for i in 1:2:(n-3)) # p_1 + p_3 + ... + p_{n-3}
            val2 = (sum_odd + (half_n - 1) * p[n-1] + half_n * p[n]) % 2
            return CongruencyClass(g, (val1, val2))
        else
            # D_{2n+1}: 2(p_1 + p_3 + ... + p_{n-2}) + (n-1)*p_{n-1} + n*p_n
            sum_odd = sum(p[i] for i in 1:2:(n-2)) # p_1 + p_3 + ... + p_{n-2}
            value = (2 * sum_odd + (n - 2) * p[n-1] + n * p[n]) % 4
            return CongruencyClass(g, value)
        end
        
    elseif g.series == :E
        if n == 6
            # E_6: p_1 - p_2 + p_4 - p_5 (mod 3)
            value = (p[1] - p[2] + p[4] - p[5]) % 3
            # Ensure positive
            value = (value + 3) % 3
            return CongruencyClass(g, value)
        elseif n == 7
            # E_7: p_4 + p_6 + p_7 (mod 2)
            value = (p[4] + p[6] + p[7]) % 2
            return CongruencyClass(g, value)
        elseif n == 8
            # E_8: trivial
            return CongruencyClass(g, nothing)
        end
        
    elseif g.series == :F || g.series == :G
        # F_4, G_2: trivial
        return CongruencyClass(g, nothing)
    end
    
    error("Congruency class not implemented for this algebra type")
end

"""
    adjoint_irrep(g::LieAlgebra)

Return the adjoint representation of the Lie algebra as an Irrep.

The adjoint representation corresponds to the Lie algebra acting on itself via the
commutator. For any simple Lie algebra, the adjoint has highest weight equal to the
highest root, which can be expressed in Dynkin labels.

# Dynkin labels for adjoint representations
- A_n: [1, 0, ..., 0, 1] (highest root)
- B_n: [0, 1, 0, ..., 0] (second fundamental weight)
- C_n: [2, 0, ..., 0] (twice the first fundamental weight)
- D_n: [0, 1, 0, ..., 0] (second fundamental weight)
- E_6: [0, 0, 0, 0, 0, 1] (sixth fundamental weight)
- E_7: [1, 0, 0, 0, 0, 0, 0] (first fundamental weight)
- E_8: [0, 0, 0, 0, 0, 0, 1, 0] (seventh fundamental weight)
- F_4: [1, 0, 0, 0] (first fundamental weight)
- G_2: [0, 1] (second fundamental weight)

# References
- Fulton & Harris, "Representation Theory" (1991), Section 21.1
- Humphreys, "Introduction to Lie Algebras and Representation Theory" (1972)
- Online database: http://www.liegroups.org/
"""
function adjoint_irrep(g::LieAlgebra)
    n = g.rank
    
    if g.series == :A
        # A_n: adjoint is [1, 0, ..., 0, 1] for n >= 2
        # A_1 is an exception: adjoint is [2]
        labels = zeros(Int, n)
        if n == 1
            labels[1] = 2
        else
            labels[1] = 1
            labels[n] = 1
        end
        return Irrep(g, labels)
        
    elseif g.series == :B
        # B_n: adjoint is [0, 1, 0, ..., 0] for n >= 3
        # B_2 is an exception: adjoint is [0, 2]
        labels = zeros(Int, n)
        if n == 2
            labels[2] = 2
        else
            labels[2] = 1
        end
        return Irrep(g, labels)
        
    elseif g.series == :C
        # C_n: adjoint is [2, 0, ..., 0]
        labels = zeros(Int, n)
        labels[1] = 2
        return Irrep(g, labels)
        
    elseif g.series == :D
        # D_n: adjoint is [0, 1, 0, ..., 0] for n >= 4
        # D_2 is an exception: adjoint is [0, 2]
        # D_3 is an exception: adjoint is [0, 1, 1]
        labels = zeros(Int, n)
        if n == 2
            labels[2] = 2
        elseif n == 3
            labels[2] = 1
            labels[3] = 1
        else
            labels[2] = 1
        end
        return Irrep(g, labels)
        
    elseif g.series == :E
        if n == 6
            # E_6: adjoint is [0, 0, 0, 0, 0, 1]
            return Irrep(g, [0, 0, 0, 0, 0, 1])
        elseif n == 7
            # E_7: adjoint is [1, 0, 0, 0, 0, 0, 0]
            return Irrep(g, [1, 0, 0, 0, 0, 0, 0])
        elseif n == 8
            # E_8: adjoint is [0, 0, 0, 0, 0, 0, 1, 0]
            return Irrep(g, [0, 0, 0, 0, 0, 0, 1, 0])
        end
        
    elseif g.series == :F && n == 4
        # F_4: adjoint is [1, 0, 0, 0]
        return Irrep(g, [1, 0, 0, 0])
        
    elseif g.series == :G && n == 2
        # G_2: adjoint is [0, 1]
        return Irrep(g, [0, 1])
    end
    
    error("Adjoint not implemented for $(g)")
end

# Global cache for irreps_up_to_dim results
# Maps (algebra) => (maxDim, result_vector)
const _IRREPS_CACHE = Dict{LieAlgebra, Tuple{Int, Vector{Irrep}}}()

"""
    irreps_up_to_dim(algebra::LieAlgebra, maxDim::Int; use_cache::Bool=true)

Return all irreducible representations of the given algebra with dimension up to `maxDim`.

Uses breadth-first search starting from the trivial representation, incrementing
Dynkin labels one at a time until the dimension exceeds `maxDim`.

Results are cached, so if a larger `maxDim` has been computed for the same algebra,
the cached result will be filtered instead of recomputing from scratch.

Set `use_cache=false` to disable caching and force recomputation.
"""
function irreps_up_to_dim(algebra::LieAlgebra, maxDim::Int; use_cache::Bool=true)
    if maxDim < 1
        return Irrep[]
    end
    
    # Check cache if enabled
    if use_cache && haskey(_IRREPS_CACHE, algebra)
        cached_maxDim, cached_result = _IRREPS_CACHE[algebra]
        
        # If we have a cached result with maxDim >= requested maxDim, filter it
        if cached_maxDim >= maxDim
            return filter(rep -> dimension(rep) <= maxDim, cached_result)
        end
        # If we have a smaller cached result, we'll extend it below
    end
    
    rank = lie_rank(algebra)
    result = Irrep[]
    visited = Set{Vector{Int}}()
    queue = Vector{Vector{Int}}()
    
    # Start with the trivial representation [0, 0, ..., 0]
    trivial_labels = zeros(Int, rank)
    push!(queue, trivial_labels)
    push!(visited, trivial_labels)
    
    # BFS
    while !isempty(queue)
        current_labels = popfirst!(queue)
        current_rep = Irrep(algebra, current_labels)
        
        # Add to result if dimension is within bounds
        if dimension(current_rep) <= maxDim
            push!(result, current_rep)
            
            # Explore neighbors: increment each Dynkin label by 1
            for i in 1:rank
                new_labels = copy(current_labels)
                new_labels[i] += 1
                
                # Only explore if not visited
                if !(new_labels in visited)
                    new_rep = Irrep(algebra, new_labels)
                    
                    # Only add to queue if dimension is still within bounds
                    if dimension(new_rep) <= maxDim
                        push!(queue, new_labels)
                        push!(visited, new_labels)
                    end
                end
            end
        end
    end
    
    # Update cache if enabled
    if use_cache
        # Only update cache if this is a larger maxDim than previously cached
        if !haskey(_IRREPS_CACHE, algebra) || _IRREPS_CACHE[algebra][1] < maxDim
            _IRREPS_CACHE[algebra] = (maxDim, result)
        end
    end
    
    return result
end

"""
    clear_irreps_cache!()

Clear the cache for irreps_up_to_dim.
"""
function clear_irreps_cache!()
    empty!(_IRREPS_CACHE)
    return nothing
end

"""
    dimension_label(rep::Irrep)

Return a string label for the representation combining dimension with distinguishing marks.

The label consists of:
- The dimension
- A bar (̄) if the representation appears after its conjugate in the canonical ordering
- Prime marks (′) to distinguish between representations with the same dimension and conjugation

For SO(8) (D_4), uses special triality labels: 8_v (vector), 8_s (spinor), 8_c (conjugate spinor).

The canonical ordering is by: increasing Dynkin index, increasing congruency class value,
then decreasing lexicographic order of Dynkin labels.
"""
function dimension_label(rep::Irrep)
    dim = dimension(rep)
    g = rep.algebra
    
    # Get all irreps with the same dimension
    same_dim_irreps = filter(r -> dimension(r) == dim, irreps_up_to_dim(g, dim))
    
    # Define ordering for congruency classes
    function cc_order(cc::CongruencyClass)
        if cc.value === nothing
            return 0
        elseif cc.value isa Tuple
            # For Z_2 × Z_2, use (a,b) -> 2a + b for ordering
            return 2 * cc.value[1] + cc.value[2]
        else
            return cc.value
        end
    end
    
    # Sort by: increasing Dynkin index, increasing congruency class, decreasing lexicographic Dynkin labels
    sorted_irreps = sort(same_dim_irreps, by = r -> (
        dynkin_index(r),
        cc_order(congruency_class(r)),
        tuple((-l for l in dynkin_labels(r))...)  # Negative for decreasing order
    ))
    
    # Check if rep appears after its conjugate
    conj_rep = conjugate(rep)
    rep_idx = findfirst(==(rep), sorted_irreps)
    conj_idx = findfirst(==(conj_rep), sorted_irreps)
    has_bar = rep_idx !== nothing && conj_idx !== nothing && rep_idx > conj_idx
    
    # Filter to remove irreps whose conjugates appear earlier
    filtered_irreps = []
    for r in sorted_irreps
        conj_r = conjugate(r)
        conj_r_idx = findfirst(==(conj_r), sorted_irreps)
        r_idx = findfirst(==(r), sorted_irreps)
        # Keep if conjugate doesn't appear earlier
        if conj_r_idx === nothing || conj_r_idx >= r_idx
            push!(filtered_irreps, r)
        end
    end
    
    # Find where rep or its conjugate appears in filtered list
    prime_count = 0
    for (i, r) in enumerate(filtered_irreps)
        if r == rep || r == conj_rep
            prime_count = i - 1
            break
        end
    end
    
    # SO(8) triality subscript logic
    subscript = ""
    is_so8_with_triality = false
    if g.series == :D && g.rank == 4
        labels = rep.dynkin_labels
        # Extract positions {1, 3, 4} (0-indexed: 0, 2, 3)
        triality_coords = [labels[1], labels[3], labels[4]]
        n_distinct = length(unique(triality_coords))
        
        if n_distinct > 1
            is_so8_with_triality = true
            p1, p3, p4 = labels[1], labels[3], labels[4]
            
            if n_distinct == 2
                # Two equal, one different
                if p3 == p4
                    subscript = "v"
                elseif p1 == p4
                    subscript = "c"
                elseif p1 == p3
                    subscript = "s"
                end
            elseif n_distinct == 3
                # All three different
                if p1 > p4 > p3
                    subscript = "vs"
                elseif p1 > p3 > p4
                    subscript = "vc"
                elseif p4 > p1 > p3
                    subscript = "sv"
                elseif p4 > p3 > p1
                    subscript = "sc"
                elseif p3 > p1 > p4
                    subscript = "cv"
                elseif p3 > p4 > p1
                    subscript = "cs"
                end
            end
        end
    end
    
    # Build the label
    label = string(dim)
    
    # SO(8) irreps with triality subscripts don't get bars or primes
    if is_so8_with_triality
        label *= "_" * subscript
    else
        if has_bar
            label *= "b"
        end
        label *= "′"^prime_count  # Unicode prime U+2032
    end
    
    return label
end
