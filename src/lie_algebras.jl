"""
    LieAlgebra

Represents a simple Lie algebra by its type and rank.
Supported types: A, B, C, D (classical) and E, F, G (exceptional).
"""
struct LieAlgebra
    series::Symbol  # :A, :B, :C, :D, :E, :F, :G
    rank::Int
    
    function LieAlgebra(series::Symbol, rank::Int)
        # Validate the series and rank
        if series == :A && rank >= 1
            new(:A, rank)
        elseif series == :B && rank >= 2
            new(:B, rank)
        elseif series == :C && rank >= 2
            new(:C, rank)
        elseif series == :D && rank >= 3
            new(:D, rank)
        elseif series == :E && rank in [6, 7, 8]
            new(:E, rank)
        elseif series == :F && rank == 4
            new(:F, rank)
        elseif series == :G && rank == 2
            new(:G, rank)
        else
            error("Invalid Lie algebra: $(series)_$(rank)")
        end
    end
end

# Convenient constructors
"""
    A_series(n::Int)

Construct the A_n Lie algebra (su(n+1)).
"""
A_series(n::Int) = LieAlgebra(:A, n)

"""
    B_series(n::Int)

Construct the B_n Lie algebra (so(2n+1)).
"""
B_series(n::Int) = LieAlgebra(:B, n)

"""
    C_series(n::Int)

Construct the C_n Lie algebra (sp(2n)).
"""
C_series(n::Int) = LieAlgebra(:C, n)

"""
    D_series(n::Int)

Construct the D_n Lie algebra (so(2n)).
"""
D_series(n::Int) = LieAlgebra(:D, n)

"""
    E_series(n::Int)

Construct the E_n Lie algebra (n = 6, 7, or 8).
"""
E_series(n::Int) = LieAlgebra(:E, n)

"""
    F_series(n::Int)

Construct the F_4 Lie algebra (n must be 4).
"""
F_series(n::Int) = LieAlgebra(:F, n)

"""
    G_series(n::Int)

Construct the G_2 Lie algebra (n must be 2).
"""
G_series(n::Int) = LieAlgebra(:G, n)

# Standard notation constructors
"""
    SU(n::Int)

Construct the su(n) Lie algebra, which is isomorphic to A_{n-1}.
For n ≥ 2, this is the Lie algebra of n×n traceless Hermitian matrices.
"""
SU(n::Int) = n >= 2 ? LieAlgebra(:A, n-1) : error("SU(n) requires n ≥ 2")

"""
    SO(n::Int)

Construct the so(n) Lie algebra.
- For n = 2k+1 ≥ 3 (odd), this is B_k
- For n = 2k ≥ 4 (even), this is D_k
so(n) is the Lie algebra of n×n skew-symmetric matrices.
"""
function SO(n::Int)
    if n < 3
        error("SO(n) requires n ≥ 3")
    elseif isodd(n)
        k = (n - 1) ÷ 2
        LieAlgebra(:B, k)
    else
        k = n ÷ 2
        LieAlgebra(:D, k)
    end
end

"""
    Sp(n::Int)

Construct the sp(n) Lie algebra.
For even n = 2k, this is C_k, the Lie algebra of n×n matrices preserving a symplectic form.
Requires even n ≥ 2.
"""
function Sp(n::Int)
    if n < 2 || isodd(n)
        error("Sp(n) requires even n ≥ 2")
    end
    k = n ÷ 2
    LieAlgebra(:C, k)
end

# Display
Base.show(io::IO, g::LieAlgebra) = print(io, "$(g.series)_$(g.rank)")

# Equality
Base.:(==)(g1::LieAlgebra, g2::LieAlgebra) = g1.series == g2.series && g1.rank == g2.rank

# Hash for using as dictionary keys
Base.hash(g::LieAlgebra, h::UInt) = hash((g.series, g.rank), h)

"""
    lie_rank(g::LieAlgebra)

Return the rank of the Lie algebra.
"""
lie_rank(g::LieAlgebra) = g.rank

"""
    dimension(g::LieAlgebra)

Return the dimension of the Lie algebra (number of generators).
"""
function dimension(g::LieAlgebra)
    adj = adjoint_irrep(g)
    return dimension(adj)
end

"""
    dual_coxeter_number(g::LieAlgebra)

Compute the dual Coxeter number hⱽ of the Lie algebra.
This is defined as the Dynkin index of the adjoint representation.
"""
function dual_coxeter_number(g::LieAlgebra)
    # The dual Coxeter number equals the Dynkin index of the adjoint representation
    adj = adjoint_irrep(g)
    return Int(round(dynkin_index(adj)))
end

const hⱽ(g::LieAlgebra) = dual_coxeter_number(g)

"""
    cartan_matrix(g::LieAlgebra)

Return the Cartan matrix for the given Lie algebra.
"""
function cartan_matrix(g::LieAlgebra)
    n = g.rank
    
    if g.series == :A
        # A_n Cartan matrix
        C = zeros(Int, n, n)
        for i in 1:n
            C[i, i] = 2
            if i < n
                C[i, i+1] = -1
                C[i+1, i] = -1
            end
        end
        return C
        
    elseif g.series == :B
        # B_n Cartan matrix
        C = zeros(Int, n, n)
        for i in 1:n
            C[i, i] = 2
            if i < n-1
                C[i, i+1] = -1
                C[i+1, i] = -1
            elseif i == n-1
                C[i, i+1] = -2
                C[i+1, i] = -1
            end
        end
        return C
        
    elseif g.series == :C
        # C_n Cartan matrix
        C = zeros(Int, n, n)
        for i in 1:n
            C[i, i] = 2
            if i < n-1
                C[i, i+1] = -1
                C[i+1, i] = -1
            elseif i == n-1
                C[i, i+1] = -1
                C[i+1, i] = -2
            end
        end
        return C
        
    elseif g.series == :D
        # D_n Cartan matrix
        C = zeros(Int, n, n)
        for i in 1:n
            C[i, i] = 2
            if i < n-1
                C[i, i+1] = -1
                C[i+1, i] = -1
            end
        end
        # Special branching: α_{n-2} connects to both α_{n-1} and α_n
        C[n-2, n] = -1
        C[n, n-2] = -1
        return C
        
    elseif g.series == :E
        # E_n Cartan matrices
        if n == 6
            return [
                2  -1   0   0   0   0;
               -1   2  -1   0   0   0;
                0  -1   2  -1   0  -1;
                0   0  -1   2  -1   0;
                0   0   0  -1   2   0;
                0   0  -1   0   0   2
            ]
        elseif n == 7
            return [
                2  -1   0   0   0   0   0;
               -1   2  -1   0   0   0   0;
                0  -1   2  -1   0   0  -1;
                0   0  -1   2  -1   0   0;
                0   0   0  -1   2  -1   0;
                0   0   0   0  -1   2   0;
                0   0  -1   0   0   0   2
            ]
        elseif n == 8
            return [
                2  -1   0   0   0   0   0   0;
               -1   2  -1   0   0   0   0   0;
                0  -1   2  -1   0   0   0  -1;
                0   0  -1   2  -1   0   0   0;
                0   0   0  -1   2  -1   0   0;
                0   0   0   0  -1   2  -1   0;
                0   0   0   0   0  -1   2   0;
                0   0  -1   0   0   0   0   2
            ]
        end
        
    elseif g.series == :F
        # F_4 Cartan matrix
        return [
            2  -1   0   0;
           -1   2  -2   0;
            0  -1   2  -1;
            0   0  -1   2
        ]
        
    elseif g.series == :G
        # G_2 Cartan matrix
        # Convention: α_1 is short, α_2 is long
        return [
            2  -1;
           -3   2
        ]
        
    end
    
    error("Unknown Lie algebra: $(g)")
end
