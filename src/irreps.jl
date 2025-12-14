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
