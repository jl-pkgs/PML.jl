"""
    tridiagonal_solver(a, b, c, d)

Tridiagonal solver

# Description

Converted into a R code from the original code of Gordon Bonan: Bonan, G.
(2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge:
Cambridge University Press. doi:10.1017/9781107339217

Solve for U given the set of equations R * U = D, where U is a vector
of length N, D is a vector of length N, and R is an N x N tridiagonal
matrix defined by the vectors A, B, C each of length N. A(1) and
C(N) are undefined and are not referenced.

    |B(1) C(1) ...  ...  ...                     |
    |A(2) B(2) C(2) ...  ...                     |
R = |     A(3) B(3) C(3) ...                     |
    |                    ... A(N-1) B(N-1) C(N-1)|
    |                    ... ...    A(N)   B(N)  |

The system of equations is written as:

  A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i

for i = 1 to N. The solution is found by rewriting the
equations so that:

  U_i = F_i - E_i * U_i+1

# Return
- `Solution: U`
"""
function tridiagonal_solver!(
  u::V, e::V, f::V,
  a::V, b::V, c::V, d::V) where {T<:Real,V<:AbstractVector{T}}

  n = length(a)
  e[1] = c[1] / b[1]
  @inbounds for i in 2:(n-1)
    e[i] = c[i] / (b[i] - a[i] * e[i-1])
  end

  f[1] = d[1] / b[1]
  @inbounds for i in 2:(n)
    f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
  end

  u[n] = f[n]
  @inbounds for i in n-1:-1:1
    u[i] = f[i] - e[i] * u[i+1]
  end
  return u
end


function tridiagonal_solver(a::V, b::V, c::V, d::V) where {T<:Real,V<:AbstractVector{T}}
  n = length(b)
  u = fill(NaN, n)
  e = fill(NaN, n) # n - 1
  f = fill(NaN, n)
  tridiagonal_solver!(u, e, f, a, b, c, d)
end


# 三角矩阵求解的临时变量
@with_kw struct TriDiagonal{FT}
  n::Int = 10
  u::Vector{FT} = zeros(FT, n) # [I↑₀, I↓₀, I↑₁, I↓₁, ..., I↑ₙ, I↓ₙ]
  a::Vector{FT} = zeros(FT, n)
  b::Vector{FT} = zeros(FT, n)
  c::Vector{FT} = zeros(FT, n)
  d::Vector{FT} = zeros(FT, n)
  e::Vector{FT} = zeros(FT, n)
  f::Vector{FT} = zeros(FT, n)
end

function tridiagonal_solver!(tri::TriDiagonal{FT}) where {FT}
  (; u, e, f, a, b, c, d) = tri
  tridiagonal_solver!(u, e, f, a, b, c, d)
end
