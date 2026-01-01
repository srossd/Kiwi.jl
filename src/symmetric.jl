"""
Symmetric group representations and character theory.

This module implements:
- Irreducible representations labeled by partitions (Young diagrams)
- Conjugacy classes labeled by cycle types (partitions)
- Character computation via the Murnaghan-Nakayama rule
"""

"""
    Partition

A partition of a positive integer n, representing a Young diagram.
Stored in decreasing order: λ₁ ≥ λ₂ ≥ ... ≥ λₖ > 0

Partitions label both irreducible representations and conjugacy classes
of the symmetric group Sₙ where n = Σᵢ λᵢ.
"""
struct Partition
    parts::Vector{Int}
    
    function Partition(parts::Vector{Int})
        # Remove zeros and sort in decreasing order
        filtered = filter(x -> x > 0, parts)
        sorted = sort(filtered, rev=true)
        new(sorted)
    end
end

# Convenience constructors
Partition(parts::Int...) = Partition(collect(parts))
Partition(n::Int) = Partition([n])  # Single row partition

# Conversion from Vector{Int}
Base.convert(::Type{Partition},v::Vector{T}) where {T<:Int} = Partition(v)

"""
    SymmetricIrrep

An irreducible representation of the symmetric group Sₙ, labeled by a partition λ.

The representation corresponds to the Specht module associated with the Young diagram λ.
"""
struct SymmetricIrrep
    partition::Partition
end

# Convenience constructors
SymmetricIrrep(parts::Vector{Int}) = SymmetricIrrep(Partition(parts))
SymmetricIrrep(parts::Int...) = SymmetricIrrep(Partition(parts...))

# Conversion from Vector{Int}
Base.convert(::Type{SymmetricIrrep}, v::Vector{T}) where {T<:Int} = SymmetricIrrep(v)

"""
    ConjugacyClass

A conjugacy class of the symmetric group Sₙ, labeled by a partition (cycle type).

A partition [λ₁, λ₂, ..., λₖ] represents elements with cycle structure:
- λ₁-cycles, λ₂-cycles, etc.

For example, [3, 2, 1] in S₆ represents permutations with one 3-cycle, one 2-cycle, 
and one fixed point.
"""
struct ConjugacyClass
    partition::Partition
end

# Convenience constructors
ConjugacyClass(parts::Vector{Int}) = ConjugacyClass(Partition(parts))
ConjugacyClass(parts::Int...) = ConjugacyClass(Partition(parts...))

# Conversion from Vector{Int}
Base.convert(::Type{ConjugacyClass}, v::Vector{T}) where {T<:Int} = ConjugacyClass(v)

# Display methods
function Base.show(io::IO, p::Partition)
    print(io, "[", join(p.parts, ", "), "]")
end

function Base.show(io::IO, irrep::SymmetricIrrep)
    n = sum(irrep.partition.parts)
    print(io, "S", n, " irrep ", irrep.partition)
end

function Base.show(io::IO, cc::ConjugacyClass)
    n = sum(cc.partition.parts)
    print(io, "S", n, " conjugacy class ", cc.partition)
end

# Equality
Base.:(==)(p1::Partition, p2::Partition) = p1.parts == p2.parts
Base.:(==)(r1::SymmetricIrrep, r2::SymmetricIrrep) = r1.partition == r2.partition
Base.:(==)(c1::ConjugacyClass, c2::ConjugacyClass) = c1.partition == c2.partition

# Hashing for use in dictionaries
Base.hash(p::Partition, h::UInt) = hash(p.parts, h)
Base.hash(r::SymmetricIrrep, h::UInt) = hash(r.partition, h)
Base.hash(c::ConjugacyClass, h::UInt) = hash(c.partition, h)

"""
    partition_size(p::Partition)

Compute n = Σᵢ λᵢ, the integer that the partition divides.
This gives the degree of the symmetric group Sₙ.
"""
function partition_size(p::Partition)
    return sum(p.parts)
end

partition_size(irrep::SymmetricIrrep) = partition_size(irrep.partition)
partition_size(cc::ConjugacyClass) = partition_size(cc.partition)

"""
    conjugacy_class_size(cc::ConjugacyClass)

Compute the size of a conjugacy class in Sₙ.

For a partition λ = [λ₁^m₁, λ₂^m₂, ..., λₖ^mₖ] (where λᵢ^mᵢ means λᵢ appears mᵢ times),
the size is:

    |C_λ| = n! / (∏ᵢ λᵢ^mᵢ · mᵢ!)

where n = Σᵢ mᵢ · λᵢ.
"""
function conjugacy_class_size(cc::ConjugacyClass)
    λ = cc.partition.parts
    n = sum(λ)
    
    # Count multiplicities: how many times each part appears
    multiplicities = Dict{Int, Int}()
    for part in λ
        multiplicities[part] = get(multiplicities, part, 0) + 1
    end
    
    # Compute denominator: ∏ᵢ λᵢ^mᵢ · mᵢ!
    denominator = 1
    for (part, mult) in multiplicities
        denominator *= part^mult * factorial(mult)
    end
    
    return factorial(n) ÷ denominator
end

"""
    hook_length(partition::Partition, i::Int, j::Int)

Compute the hook length at position (i,j) in a Young diagram.

The hook length h(i,j) is:
    h(i,j) = λᵢ - j + λ'ⱼ - i + 1

where λ'ⱼ is the length of the j-th column (i.e., the number of parts ≥ j).
"""
function hook_length(partition::Partition, i::Int, j::Int)
    λ = partition.parts
    
    # Check if (i,j) is in the diagram
    if i > length(λ) || j > λ[i]
        return 0
    end
    
    # Arm length: cells to the right in row i
    arm = λ[i] - j
    
    # Leg length: cells below in column j
    leg = count(k -> k >= j, λ[i+1:end])
    
    return arm + leg + 1
end

"""
    dimension(irrep::SymmetricIrrep)

Compute the dimension of a symmetric group irrep using the hook length formula.

    dim(λ) = n! / ∏_{(i,j)∈λ} h(i,j)

where h(i,j) is the hook length at position (i,j).
"""
function dimension(irrep::SymmetricIrrep)
    λ = irrep.partition.parts
    n = sum(λ)
    
    # Compute product of hook lengths
    hook_product = 1
    for i in 1:length(λ)
        for j in 1:λ[i]
            hook_product *= hook_length(irrep.partition, i, j)
        end
    end
    
    return factorial(n) ÷ hook_product
end

"""
    remove_rim_hook(partition::Partition, length::Int)

Find all ways to remove a rim hook (border strip) of given length from a partition.

A rim hook is a connected set of boundary cells that doesn't contain a 2×2 block.
Returns a vector of tuples (new_partition, sign) where sign = (-1)^(height-1)
and height is the number of rows the rim hook spans.

This is used in the Murnaghan-Nakayama rule.
"""
function remove_rim_hook(partition::Partition, target_length::Int)
    λ = partition.parts
    if isempty(λ) || target_length <= 0
        return Tuple{Partition, Int}[]
    end
    
    if target_length > sum(λ)
        return Tuple{Partition, Int}[]
    end
    
    results = Tuple{Partition, Int}[]
    
    # Build the maximal rim (boundary path from top-right to bottom-left)
    # Starting at (1, λ[1]), we move down if possible, otherwise left
    rim_cells = Tuple{Int, Int}[]
    i, j = 1, λ[1]
    
    while i <= length(λ) && j >= 1
        push!(rim_cells, (i, j))
        if j == 1 && i == length(λ)
            break
        end
        
        # Try to move down (increase row)
        if i < length(λ) && λ[i+1] >= j
            i += 1
        else
            # Move left (decrease column)
            j -= 1
        end
    end
    
    # Now search for contiguous strips of target_length along the rim
    # Each strip must be a valid rim hook (connected and no 2×2 blocks)
    for start_idx in 1:length(rim_cells)
        if start_idx + target_length - 1 > length(rim_cells)
            break
        end
        
        # Extract the potential rim hook
        hook_cells = rim_cells[start_idx:(start_idx + target_length - 1)]

        # Try to remove these cells and see if we get a valid partition
        new_parts = copy(λ)
        
        # Group cells by row
        cells_by_row = Dict{Int, Vector{Int}}()
        for (ri, rj) in hook_cells
            if !haskey(cells_by_row, ri)
                cells_by_row[ri] = Int[]
            end
            push!(cells_by_row[ri], rj)
        end
        
        # For each row, we need to remove cells from the right end
        # The cells must be contiguous and at the end of the row
        valid = true
        for (row_idx, cols) in cells_by_row
            if maximum(cols) != new_parts[row_idx]
                valid = false
                break
            end
            # Remove the cells
            new_parts[row_idx] -= length(cols)
        end

        # Check that the result is still a valid partition (decreasing)
        for k in 1:(length(new_parts)-1)
            if new_parts[k] < new_parts[k+1]
                valid = false
                break
            end
        end
        
        if !valid
            continue
        end
        
        # Calculate height (number of distinct rows in the hook)
        rows_in_hook = Set(cell[1] for cell in hook_cells)
        height = length(rows_in_hook)
        sign = (-1)^(height - 1)
        
        # Filter out zero parts and create new partition
        filtered_parts = filter(x -> x > 0, new_parts)
        new_partition = Partition(filtered_parts)
        
        push!(results, (new_partition, sign))
    end
    
    return unique(results)
end

"""
    murnaghan_nakayama(irrep::SymmetricIrrep, cc::ConjugacyClass)

Compute the character value χ^λ(μ) using the Murnaghan-Nakayama rule.

The character of irrep λ evaluated on conjugacy class μ = [μ₁, μ₂, ..., μₖ] is:

    χ^λ(μ) = Σ sign(T) 

summed over all ways to sequentially remove rim hooks of lengths μ₁, μ₂, ..., μₖ
from the Young diagram λ, where sign(T) is the product of signs from each removal.

Returns 0 if the irrep and conjugacy class correspond to different values of n.
"""
function murnaghan_nakayama(irrep::SymmetricIrrep, cc::ConjugacyClass)
    λ = irrep.partition
    μ = cc.partition.parts
    
    # Check that both are for the same Sₙ
    if partition_size(λ) != partition_size(cc.partition)
        return 0
    end
    
    # Base case: empty partition
    if isempty(μ) || all(x -> x == 0, μ)
        return isempty(λ.parts) || all(x -> x == 0, λ.parts) ? 1 : 0
    end
    
    # Recursive case: remove rim hooks of length μ[1]
    # then compute character for remaining parts
    total = 0
    hook_length = μ[1]
    remaining_cycle_type = length(μ) > 1 ? μ[2:end] : Int[]
    
    removals = remove_rim_hook(λ, hook_length)
    
    for (new_partition, sign) in removals
        # Recursively compute character for the remaining diagram and cycle type
        remaining_cc = ConjugacyClass(Partition(remaining_cycle_type))
        remaining_irrep = SymmetricIrrep(new_partition)
        
        sub_char = murnaghan_nakayama(remaining_irrep, remaining_cc)
        total += sign * sub_char
    end
    
    return total
end

"""
    character(irrep::SymmetricIrrep, cc::ConjugacyClass)

Compute the character value χ^λ(μ) of irrep λ on conjugacy class μ.
Alias for `murnaghan_nakayama`.
"""
character(irrep::SymmetricIrrep, cc::ConjugacyClass) = murnaghan_nakayama(irrep, cc)

"""
    character(irrep::SymmetricIrrep; lazy=false)

Create a SymmetricCharacter for the given irrep.

# Arguments
- `irrep::SymmetricIrrep`: The irreducible representation
- `lazy::Bool=false`: If true, computes character values on demand.
                      If false (default), precomputes all values.

# Returns
A `SymmetricCharacter` object.

# Examples
```julia
# Eager evaluation (default): precomputes all values
irrep = SymmetricIrrep([2, 1])
char = character(irrep)

# Lazy evaluation: computes on demand
char_lazy = character(irrep; lazy=true)
```
"""
character(irrep::SymmetricIrrep; lazy::Bool=false) = SymmetricCharacter(irrep; lazy=lazy)

"""
    character_table(n::Int)

Compute the full character table of Sₙ.

Returns a dictionary mapping (irrep, conjugacy_class) pairs to character values.
"""
function character_table(n::Int)
    # Generate all partitions of n
    partitions = all_partitions(n)
    
    # Build character table
    table = Dict{Tuple{SymmetricIrrep, ConjugacyClass}, Int}()
    
    for λ in partitions
        irrep = SymmetricIrrep(λ)
        for μ in partitions
            cc = ConjugacyClass(μ)
            table[(irrep, cc)] = character(irrep, cc)
        end
    end
    
    return table
end

"""
    all_partitions(n::Int)

Generate all partitions of n in decreasing lexicographic order.
"""
function all_partitions(n::Int)
    if n <= 0
        return [Partition(Int[])]
    end
    
    partitions = Partition[]
    
    function generate(remaining::Int, max_part::Int, current::Vector{Int})
        if remaining == 0
            push!(partitions, Partition(copy(current)))
            return
        end
        
        for part in min(remaining, max_part):-1:1
            push!(current, part)
            generate(remaining - part, part, current)
            pop!(current)
        end
    end
    
    generate(n, n, Int[])
    return partitions
end

"""
    SymmetricCharacter

A character of a symmetric group irreducible representation.
Stores character values for conjugacy classes, with eager or lazy evaluation.

# Fields
- `irrep::SymmetricIrrep`: The irreducible representation
- `values::Dict{ConjugacyClass, Int}`: Memoized character values
- `lazy::Bool`: Whether to compute values on demand (true) or precompute all (false)

# Constructors
    SymmetricCharacter(irrep::SymmetricIrrep; lazy=false)

Create a character for the given irrep. By default (`lazy=false`), precomputes
all character values for all conjugacy classes of Sₙ. With `lazy=true`, computes
values on demand when accessed.

# Indexing
    char[cc::ConjugacyClass]  # Get character value for conjugacy class
    char[cycle_type::Vector{Int}]  # Convenience: convert vector to ConjugacyClass

# Example
```julia
# Eager evaluation (default): precomputes all values
irrep = SymmetricIrrep([2, 1])  # Standard representation of S₃
char = SymmetricCharacter(irrep)
char[ConjugacyClass([3])]      # => 0
char[[2, 1]]                    # => -1

# Lazy evaluation: computes on demand
char_lazy = SymmetricCharacter(irrep; lazy=true)
char_lazy[[1, 1, 1]]            # => 2 (computed and cached)
```
"""
mutable struct SymmetricCharacter
    irrep::SymmetricIrrep
    values::Dict{ConjugacyClass, Int}
    lazy::Bool
    
    function SymmetricCharacter(irrep::SymmetricIrrep; lazy::Bool=false)
        values = Dict{ConjugacyClass, Int}()
        
        if !lazy
            # Precompute all character values for this Sₙ
            n = partition_size(irrep)
            partitions = all_partitions(n)
            
            for μ in partitions
                cc = ConjugacyClass(μ)
                values[cc] = murnaghan_nakayama(irrep, cc)
            end
        end
        
        new(irrep, values, lazy)
    end
end

# Display
function Base.show(io::IO, char::SymmetricCharacter)
    n = partition_size(char.irrep)
    if char.lazy
        n_computed = length(char.values)
        print(io, "SymmetricCharacter for ", char.irrep, " (", n_computed, " values computed)")
    else
        n_classes = length(all_partitions(n))
        print(io, "SymmetricCharacter for ", char.irrep, " (", n_classes, " values precomputed)")
    end
end

# Indexing by ConjugacyClass
function Base.getindex(char::SymmetricCharacter, cc::ConjugacyClass)
    # Check if already computed
    if haskey(char.values, cc)
        return char.values[cc]
    end
    
    # Compute and cache
    value = murnaghan_nakayama(char.irrep, cc)
    char.values[cc] = value
    return value
end

# Convenience indexing by Vector{Int}
function Base.getindex(char::SymmetricCharacter, cycle_type::Vector{<:Integer})
    cc = ConjugacyClass(cycle_type)
    return char[cc]
end

"""
    dimension(char::SymmetricCharacter)

Compute the dimension of the representation.
"""
function dimension(char::SymmetricCharacter)
    return dimension(char.irrep)
end
