# Recoupling Coefficients

Kiwi computes explicit representation matrices, Clebsch–Gordan coefficients and
recoupling symbols (F/6j, R, A, B, 1j, 2j, 3j) for every compact simple Lie
algebra, including tensor products with outer multiplicity.  This page states
the conventions precisely: individual coefficients are only meaningful together
with them.

## Quick tour

```@example rc
using Kiwi
su3 = SU(3)
a = Irrep(su3, 1, 1)                     # 8
nsymbol(a, a, a)                         # 8 appears twice in 8 ⊗ 8
```

```@example rc
C = fusion_tensor(a, a, a)               # ⟨8 i; 8 j | 8 k; μ⟩, size (8, 8, 8, 2)
size(C)
```

```@example rc
F = fsymbol(a, a, a, a, a, a)            # F^{888}_8[8, 8], size (2, 2, 2, 2)
size(F)
```

```@example rc
E, Fm, H = generators(Irrep(SU(2), 2))  # spin 1: E = J₊, F = J₋, H = 2J_z
Matrix(Fm[1])
```

## Conventions

### Bases of irreducible modules

* Weights of ``V(λ)`` are ordered by depth (the height of ``λ − μ``) and, within a
  depth, by decreasing Dynkin labels.  The highest weight comes first, and
  [`basis_weights`](@ref) lists the weight of every basis vector.
* **Lowering basis** (`basis = :lowering`, the default, all algebras).  The
  candidates spanning ``V(λ)_μ`` are the vectors ``F_i b`` for ``i = 1, …, n`` in
  turn and ``b`` running over the basis of ``V(λ)_{μ+α_i}``.  The basis is the first
  independent subfamily of candidates, chosen greedily in this order, and then
  Gram–Schmidt orthonormalised in the same order with respect to the
  contravariant form (``⟨v_λ, v_λ⟩ = 1``, ``F_i^† = E_i``).  Each basis vector
  therefore has positive overlap with the candidate it came from.  For ``\mathfrak{su}(2)``
  this is the Condon–Shortley basis, with ``F = J_-``, ``E = J_+`` and ``H = 2J_z``.
* **Gelfand–Tsetlin basis** (`basis = :gelfand_tsetlin`, ``\mathfrak{su}(N)`` only): the
  standard GT basis with nonnegative raising-operator matrix elements, the same basis
  as in SUNRepresentations.jl.  Within a weight space, patterns are ordered by
  decreasing rows ``(m_{N-1}, …, m_1)``.
* Generators are the Chevalley generators ``E_i, F_i, H_i``, with
  ``[E_i, F_j] = δ_{ij} H_i`` and ``[H_i, E_j] = ⟨α_j, α_i^∨⟩ E_j``.  In both bases
  ``E_i = F_i^T``.
* With `T = Rational{BigInt}` the module is built in exact arithmetic, in the
  unnormalised lowering basis (with its Gram matrices).  Its Gram–Schmidt
  orthonormalisation is exactly the floating-point lowering basis.

### Clebsch–Gordan coefficients and the multiplicity convention

`fusion_tensor(a, b, c)[i, j, k, μ] = ⟨a i; b j | c k; μ⟩` uses TensorKit's layout.
For each ``μ`` the map ``e_k ↦ \sum_{ij} C[i,j,k,μ]\, e_i ⊗ e_j`` is an isometric
intertwiner ``V_c → V_a ⊗ V_b``.

The intertwiners are fixed by their highest-weight vectors ``X``.  Let ``K`` be the
space of highest-weight vectors of weight ``c`` in ``V_a ⊗ V_b``, and
``f_1, f_2, …`` the basis of ``V_b(c − a)``.  The basis ``X_1, …, X_N`` is the
Gram–Schmidt orthonormalisation of the projections ``Π_K(v_a ⊗ f_j)``, taken for
``j = 1, 2, …`` in order and skipping those that depend on earlier ones.  So
``⟨X_μ, v_a ⊗ f_{j_μ}⟩ > 0``.  In the multiplicity-free case, the first nonzero
coefficient of ``X`` along ``v_a ⊗ V_b(c − a)`` is positive.  For ``\mathfrak{su}(2)`` this is
the Condon–Shortley convention ``⟨j_1 j_1; j_2\, (J − j_1) | J J⟩ > 0``, and Kiwi
reproduces the Racah formula.

`fusion_tensor(a, b, c; sparse = true)` returns the same coefficients as a sparse
``(\dim a \cdot \dim b) × (\dim c \cdot N)`` matrix.

### F- and R-symbols

Both follow TensorKitSectors.jl:

* ``F^{abc}_d[e, f]`` is an array of size ``(N^{ab}_e, N^{ec}_d, N^{bc}_f, N^{af}_d)``
  relating ``((a ⊗ b → e) ⊗ c → d)`` and ``(a ⊗ (b ⊗ c → f) → d)``:
  ```math
  F[μ,ν,κ,λ] = \sum C^{af}_d[i,n,1,λ]\, C^{bc}_f[j,l,n,κ]\, C^{ab}_e[i,j,m,μ]\, C^{ec}_d[m,l,1,ν].
  ```
* ``R^{ab}_c[μ, ν] = \sum_{ij} C^{ba}_c[j,i,1,ν]\, C^{ab}_c[i,j,1,μ]``.

Kiwi evaluates these contractions on highest-weight vectors (index 1), without
forming the full Clebsch–Gordan tensors.  For ``\mathfrak{su}(2)``,
```math
F^{j_1 j_2 j_3}_j[j_{12}, j_{23}] = (-1)^{j_1+j_2+j_3+j}\sqrt{(2j_{12}+1)(2j_{23}+1)}
\begin{Bmatrix} j_1 & j_2 & j_{12} \\ j_3 & j & j_{23}\end{Bmatrix},
\qquad R^{j_1 j_2}_J = (-1)^{j_1+j_2-J}.
```

### 1j, 2j and 3j symbols

* **2j symbol.** ``Z_a = \sqrt{\dim a}\; C^{a\bar a}_1[:, :, 1, 1]`` ([`twoj`](@ref))
  is a real orthogonal matrix (rows: basis of ``a``, columns: basis of ``\bar a``).
  The multiplicity convention gives ``Z_a[\mathrm{hw}(a), \mathrm{lw}(\bar a)] > 0``.
  For ``\mathfrak{su}(2)``, ``Z_j[m, m'] = (-1)^{j-m} δ_{m', -m}``.
* **1j symbol.** By Schur's lemma ``Z_{\bar a} = ϰ_a Z_a^T`` with ``ϰ_a = ±1``
  ([`onej`](@ref)).  For self-conjugate ``a``, ``ϰ_a`` is the Frobenius–Schur
  indicator (``+1`` real, ``-1`` pseudoreal).  For ``a ≠ \bar a``, ``ϰ_a = ϰ_{\bar a}``
  depends on the convention.  In all cases ``ϰ_a`` equals TensorKit's
  Frobenius–Schur phase ``\operatorname{sign} F^{a\bar a a}_a[1, 1]``.  This is what
  makes TensorKit's line-bending A- and B-symbols ([`asymbol`](@ref),
  [`bsymbol`](@ref)) agree whether they are computed from F-symbols or from fusion
  tensors.
* **3j symbol.**
  ``(a\,b\,c)[i,j,k,μ] = \sum_m C^{ab}_{\bar c}[i,j,m,μ]\, Z_c[k,m] / \sqrt{\dim c}``
  ([`threej`](@ref)) is an invariant of ``V_a ⊗ V_b ⊗ V_c`` normalised to
  ``\sum |·|^2 = 1``.  For ``\mathfrak{su}(2)`` it is related to Wigner's symbol by
  ``\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}_{\text{Wigner}} = (-1)^{j_1-j_2-j_3}\,(j_1\,j_2\,j_3)_{\text{Kiwi}}``.

## Exact values

With `T = Rational{BigInt}` everything is computed in exact arithmetic.  F- and
R-symbols of multiplicity-free vertices are returned as [`SqrtRational`](@ref)
values ``s\sqrt{q}``:

```@example rc
q, qb = Irrep(su3, 1, 0), Irrep(su3, 0, 1)        # 3, 3̄
fsymbol(q, q, qb, q, qb, a; T = Rational{BigInt})  # F^{3 3 3̄}_3[3̄, 8]
```

## TensorKit

After `using TensorKitSectors` (or TensorKit), `tensorkit_sector(g)` is a sector type
whose topological data are Kiwi's, so tensors with the symmetry of any compact
simple Lie group can be used in TensorKit:

```julia
using Kiwi, TensorKitSectors
I = tensorkit_sector(Kiwi.SU(3))
s = tensorkit_sector(Irrep(Kiwi.SU(3), 1, 1))
Fsymbol(s, s, s, s, s, s)
```

Kiwi's sectors pass TensorKitSectors' own test suite: pentagon, hexagon, triangle,
unitarity, and the fusion-tensor/F-symbol consistency of the A, B and
Frobenius–Schur data.
