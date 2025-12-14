"""
    dimension(rep::Irrep)

Compute the dimension of an irreducible representation using Weyl's dimension formula:
dim(λ) = ∏_{α > 0} (λ+ρ, α) / (ρ, α)
where the product is over all positive roots α, and ρ is the Weyl vector.
"""
function dimension(rep::Irrep)
    g = rep.algebra
    
    # Get highest weight, Weyl vector, and positive roots
    λ = highest_weight(rep)
    ρ = weyl_vector(g)
    roots = positive_roots(g)
    
    # Shifted weight λ + ρ
    λ_shifted = λ + ρ
    
    # Compute product using Weyl's formula
    numerator = big(1)
    denominator = big(1)
    
    for α in roots
        # Compute (λ+ρ, α) and (ρ, α)
        num_val = inner_product(λ_shifted, α)
        den_val = inner_product(ρ, α)
        
        # Multiply by rational numbers
        numerator *= Base.numerator(num_val) * Base.denominator(den_val)
        denominator *= Base.denominator(num_val) * Base.numerator(den_val)
    end
    
    result = div(numerator, denominator)
    return Int(result)
end

"""
    quadratic_casimir(rep::Irrep)

Compute the quadratic Casimir invariant C_2(rep).
The eigenvalue is given by (λ, λ + 2ρ) / 2 where ρ is the Weyl vector.
This normalization is standard in physics literature.
"""
function quadratic_casimir(rep::Irrep)
    g = rep.algebra
    λ = highest_weight(rep)
    ρ = weyl_vector(g)
    
    # Compute (λ, λ + 2ρ) / 2
    λ_shifted = λ + 2 * ρ
    
    return inner_product(λ, λ_shifted) / 2
end

"""
    dynkin_index(rep::Irrep)

Compute the Dynkin index (also known as the second index) of the representation.
For the adjoint representation, this equals 2h^∨ (twice the dual Coxeter number).
"""
function dynkin_index(rep::Irrep)
    g = rep.algebra
    
    d = dimension(rep)
    c2 = quadratic_casimir(rep)
    dim_g = dimension(g)
    
    # Dynkin index formula
    return d * c2 / dim_g
end