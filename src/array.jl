using LinearAlgebra
import Base.@kwdef

export steering, beamform
export Bartlett, Capon, Music

abstract type Beamformer end

"""
Frequency-domain Bartlett beamformer.
"""
struct Bartlett <: Beamformer end

"""
Frequency-domain Capon beamformer with diagonal loading factor `γ`.
"""
@kwdef struct Capon <: Beamformer
  γ = 0.0
end

"""
Frequency-domain MUSIC beamformer with `nsignals` signals.
"""
@kwdef struct Music <: Beamformer
  nsignals = 1
end

"""
$(SIGNATURES)
Computes steering delays for specified receiver positions `rxpos`, signal
propagation speed `c`, and angles `θ`. The delays are computed with a far-field
assumption, i.e., for plane incoming waves.

# Examples:
```julia-repl
julia> steering(0.0:1.0:5.0, 1500.0, range(0.0, π; length=181))
6×181 Array{Float64,2}:
  0.00166667    0.00166641   …  -0.00166641   -0.00166667
  0.001         0.000999848     -0.000999848  -0.001
  0.000333333   0.000333283     -0.000333283  -0.000333333
 -0.000333333  -0.000333283      0.000333283   0.000333333
 -0.001        -0.000999848      0.000999848   0.001
 -0.00166667   -0.00166641   …   0.00166641    0.00166667

julia> rxpos = [  # can be 2D or 3D coordinates
  0.0  1.0  2.0  3.0  4.0  5.0
  0.0  0.0  0.0  0.0  0.0  0.0
];
julia> steering(rxpos, 1500.0, range(0.0, π; length=181))
6×181 Array{Float64,2}:
  0.00166667    0.00166641   …  -0.00166641   -0.00166667
  0.001         0.000999848     -0.000999848  -0.001
  0.000333333   0.000333283     -0.000333283  -0.000333333
 -0.000333333  -0.000333283      0.000333283   0.000333333
 -0.001        -0.000999848      0.000999848   0.001
 -0.00166667   -0.00166641   …   0.00166641    0.00166667
```
"""
function steering(rxpos::AbstractMatrix, c, θ)
  # TODO: add support for 2D steering
  nsensors = size(rxpos, 2)
  ndir = length(θ)
  pos0 = sum(rxpos; dims=2) / nsensors
  rxpos = rxpos .- pos0
  size(rxpos, 1) == 1 && (rxpos = vcat(rxpos, zeros(1, nsensors)))
  size(rxpos, 1) > 2 && (rxpos = rxpos[1:2,:])
  sd = zeros(nsensors, ndir)
  for k ∈ 1:ndir
    for j ∈ 1:nsensors
      sd[j,k] = -rxpos[:,j]' * [cos(θ[k]), sin(θ[k])] / c
    end
  end
  sd
end

steering(rxpos::AbstractVector, c, θ) = steering(collect(rxpos'), c, θ)

"""
    beamform(s, sd; fs=framerate(s))

Broadband time-domain delay-and-sum beamformer. Takes in passband or baseband
signals `s` and produces beamformer output for all directions specified by the
steering distances `sd`. The beamformer output is a timeseries signal for each
steering direction.

# Example:
```julia-repl
julia> x = cw(100.0, 1.0, 44100.0);
julia> sd = steering(0.0:1.0:3.0, 1500.0, range(0.0, π; length=181));
julia> bfo = beamform([x x x x], sd)
44100-frame, 181-channel SampleBuf{Complex{Float64}, 2}
1.0s sampled at 44100.0Hz
  ⋮
```
"""
function beamform(s, sd; fs=framerate(s))
  ndir = size(sd, 2)
  bfo = zeros(eltype(s), (nframes(s), ndir))
  for k ∈ 1:ndir
    for j ∈ 1:nchannels(s)
      bfo[:,k] .+= padded(s[:,j], 0; delay=round(Int, -sd[j,k]*fs))
    end
  end
  signal(bfo, fs)
end

"""
    beamform(s, f, n, sd; fs=framerate(s), method=Bartlett())
    beamform(s, f, sd; fs=framerate(s), method=Bartlett())

Narrowband frequency-domain beamformer. Takes in passband signals `s` and
produces beamformer output for all directions specified by the steering
distances `sd`. The beamformer output is an energy estimate (or equivalent)
for each steering direction. The beamforming only uses a narrowband signal
cenetered at frequency `f` with a bandwidth of about `fs/n`.

If `n` is not specified, or is 1, then the input signal is assumed to be
narrowband, and centered at frequency `f`.

The narrowband assumption requires that the bandwidth be no greater than
about 5/T, where T is the maximum time taken for a signal to propagate through
the array.

Several beamforming methods are available:
- [`Bartlett`](@ref)
- [`Capon`](@ref)
- [`Music`](@ref)

Custom beamformers can be implemented by creating a subtype of
`SignalAnalysis.Beamformer` and implementing the `SignalAnalysis.beamformer()`
method dispatched on that type.

# Example:
```julia-repl
julia> x = cw(100.0, 1.0, 44100.0);
julia> sd = steering(0.0:1.0:3.0, 1500.0, range(0.0, π; length=181));
julia> bfo = beamform([x x x x], 100.0, 4096, sd; method=Capon(0.1))
181-element Array{Float64,1}:
 0.12406290296318974
 0.1240975045516605
 ⋮
 0.12406290296318974
```
"""
function beamform(s, f, n, sd; fs=framerate(s), method=Bartlett())
  s̃ = n > 1 ? goertzel(s, f, n; fs=fs) : signal(s, fs)
  diff([extrema(abs.(sd))...])[1] > 5/framerate(s̃) && @warn "Narrowband assumption not met, try increasing n"
  R = samples(s̃)' * samples(s̃)
  sv = cis.(2π * f * sd)
  beamformer(method, R, sv)
end

beamform(s, f, sd; fs=framerate(s), method=Bartlett()) = beamform(s, f, 1, sd; fs=fs, method=method)

# TODO: broadband beamformer

### various frequency-domain beamformer implementations
#
# beamformer(opt::type, R, sv)
#   R: narrowband covariance matrix
#   sv: narrowband steering vector
#   opt: beamformer-specific options

function beamformer(::Bartlett, R, sv)
  # TODO: add support for shading
  [abs.(sv[:,j]' * R * sv[:,j]) for j ∈ 1:size(sv,2)]
end

function beamformer(opt::Capon, R, sv)
  opt.γ != 0 && (R .+= opt.γ * I(size(R, 1)))
  [1 ./ abs.(sv[:,j]' * inv(R) * sv[:,j]) for j ∈ 1:size(sv,2)]
end

function beamformer(opt::Music, R, sv)
  F = eigen(R)
  ndx = sortperm(F.values; rev=true)
  Q = F.vectors[:, ndx[opt.nsignals+1:end]]
  Q = Q * Q'
  [1 ./ abs.(sv[:,j]' * Q * sv[:,j]) for j ∈ 1:size(sv,2)]
end
