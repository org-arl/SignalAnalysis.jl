using Random

export RedGaussian, PinkGaussian

# We deliberately avoid subtyping Distributions.ContinuousMultivariateDistribution
# and define our own NoiseDistribution instead: Distributions pulls in a large
# transitive dependency tree (StatsFuns, Rmath, PDMats, ...), and libRmath-julia's
# TLS segment breaks library loading on old-glibc systems by exhausting the static
# TLS surplus (https://github.com/JuliaStats/Rmath-julia/issues/56, see also
# https://github.com/org-arl/SignalAnalysis.jl/issues/50). All we need from
# Distributions is the family of rand()/rand!() fallbacks, which we define
# ourselves below.
abstract type NoiseDistribution end

"""
$(TYPEDEF)
Red Gaussian noise distribution for random variate generation.

# Example:
```julia-repl
julia> rand(RedGaussian(1000))
1000-element Array{Float64,1}:
   ⋮

julia> rand(RedGaussian(1000, 2.0))
1000-element Array{Float64,1}:
   ⋮
```
"""
Base.@kwdef struct RedGaussian{T} <: NoiseDistribution
  n::Int
  σ::T = 1.0
end

RedGaussian(n) = RedGaussian(n=n)

"""
$(TYPEDEF)
Pink Gaussian noise distribution for random variate generation.

# Example:
```julia-repl
julia> rand(PinkGaussian(1000))
1000-element Array{Float64,1}:
   ⋮

julia> rand(PinkGaussian(1000, 2.0))
1000-element Array{Float64,1}:
   ⋮
```
"""
Base.@kwdef struct PinkGaussian{T} <: NoiseDistribution
  n::Int
  σ::T = 1.0
end

PinkGaussian(n) = PinkGaussian(n=n)

function Random.rand!(rng::AbstractRNG, d::RedGaussian{T}, x::AbstractVector{T}) where {T<:Real}
  length(x) >= d.n || throw(ArgumentError("length of x must be at least n"))
  extra = 100
  v = d.σ .* randn(rng, T, length(x)+extra)
  for j = 2:length(v)
    v[j] += v[j-1]
  end
  removedc!(v)
  x .= v[extra+1:end]/3.20377
  return x
end

function Random.rand!(rng::AbstractRNG, d::PinkGaussian{T}, x::AbstractVector{T}) where {T<:Real}
  length(x) >= d.n || throw(ArgumentError("length of x must be at least n"))
  # based on https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html
  hb = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
  ha = [1.0, -2.494956002, 2.017265875, -0.522189400]
  extra = 1430
  v = d.σ .* randn(rng, T, length(x)+extra)
  v = filt(hb, ha, v)/0.08680587859687908
  x .= v[extra+1:end]
  return x
end

Base.length(d::NoiseDistribution) = d.n

Base.rand(rng::AbstractRNG, d::NoiseDistribution) = rand!(rng, d, zeros(typeof(d.σ), d.n))
Base.rand(d::NoiseDistribution) = rand(Random.default_rng(), d)
Random.rand!(d::NoiseDistribution, x::AbstractVecOrMat) = rand!(Random.default_rng(), d, x)
Base.rand(rng::AbstractRNG, d::NoiseDistribution, m::Integer) = rand!(rng, d, zeros(typeof(d.σ), d.n, m))
Base.rand(d::NoiseDistribution, m::Integer) = rand(Random.default_rng(), d, m)

function Random.rand!(rng::AbstractRNG, d::NoiseDistribution, A::AbstractMatrix)
  foreach(x -> rand!(rng, d, x), eachcol(A))
  return A
end
