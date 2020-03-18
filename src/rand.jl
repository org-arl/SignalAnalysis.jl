export RedGaussian, PinkGaussian

# TODO: doesn't handle rand(::RedGaussian, m) for m samples of n-vectors yet (see _rand!)

"""
  $(TYPEDEF)
Red Gaussian noise distribution for random variate generation.

# Example:
```jldoctest
julia> using SignalAnalysis
julia> rand(RedGaussian(1000))
1000-element Array{Float64,1}:
[...]

julia> rand(RedGaussian(1000, 2.0))
1000-element Array{Float64,1}:
[...]
"""
Base.@kwdef struct RedGaussian{T} <: Distributions.ContinuousMultivariateDistribution
  n::Int
  σ::T = 1.0
end

RedGaussian(n) = RedGaussian(n=n)

"""
  $(TYPEDEF)
Pink Gaussian noise distribution for random variate generation.

# Example:
```jldoctest
julia> using SignalAnalysis
julia> rand(PinkGaussian(1000))
1000-element Array{Float64,1}:
[...]

julia> rand(PinkGaussian(1000, 2.0))
1000-element Array{Float64,1}:
[...]
"""
Base.@kwdef struct PinkGaussian{T} <: Distributions.ContinuousMultivariateDistribution
  n::Int
  σ::T = 1.0
end

PinkGaussian(n) = PinkGaussian(n=n)

function Random.rand!(rng::AbstractRNG, d::RedGaussian{T}, x::AbstractVector{T}) where {T<:Real}
  length(x) >= d.n || throw(ArgumentError("length of x must be at least n"))
  extra = 100
  v = rand(rng, Normal{T}(0.0,d.σ), length(x)+extra)
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
  v = rand(rng, Normal{T}(0.0,d.σ), length(x)+extra)
  v = filt(hb, ha, v)/0.08680587859687908
  x .= v[extra+1:end]
  return x
end

Base.length(d::RedGaussian) = d.n
Base.length(d::PinkGaussian) = d.n

Base.rand(rng::AbstractRNG, d::RedGaussian) = Random.rand!(rng, d, zeros(d.n))
Base.rand(rng::AbstractRNG, d::PinkGaussian) = Random.rand!(rng, d, zeros(d.n))
