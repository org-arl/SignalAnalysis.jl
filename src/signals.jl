using MetaArrays
using Base.Iterators: partition, flatten

using DSP, DSP.Filters
using PaddedViews
using ProgressMeter

export signal, analytic, isanalytic, samples
export padded, toframe, domain
export partition, flatten

struct SamplingInfo
  fs::Float32
end

const SampledSignal = MetaArray{<:Any,SamplingInfo,T} where T
const SampledSignalVector = MetaArray{<:AbstractVector,SamplingInfo,T} where T
const SampledSignalMatrix = MetaArray{<:AbstractMatrix,SamplingInfo,T} where T

SignalBase.nframes(x::AbstractArray) = size(x, 1)
SignalBase.framerate(x::SampledSignal) = x.fs
SignalBase.framerate(x::AbstractArray) = 1.0
SignalBase.nchannels(x::AbstractArray) = size(x, 2)
SignalBase.sampletype(x::AbstractArray) = eltype(x)

"""
$(SIGNATURES)
Returns the domain of the signal.
"""
domain(x) = (0:nframes(x)-1) ./ framerate(x)

function Base.show(io::IO, mime::MIME"text/plain", s::SampledSignal{T}) where T
  print(io, "SampledSignal @ ", framerate(s), " Hz, ")
  show(io, mime, samples(s))
end

"""
$(SIGNATURES)
Creates a signal with frame rate `fs`.
"""
signal(x::AbstractArray, fs) = MetaArray(SamplingInfo(inHz(fs)), x)
signal(x::Base.Iterators.Flatten, fs) = signal(collect(x), fs)

"""
$(SIGNATURES)
Creates a signal with frame rate `fs`. If the original signal's frame rate is
the same as `fs`, this method simply returns the original signal. Otherwise, it
creates a new signal with the specified frame rate and data from the original
signal. Do note that this method does not resample the signal.
"""
signal(x::SampledSignal, fs) = fs == framerate(x) ? x : signal(samples(x), fs)

"""
$(SIGNATURES)
Creates an empty signal of length `n` samples, and frame rate `fs`.
"""
signal(n::Int, fs) = signal(Array{Float64}(undef, n), fs)

"""
$(SIGNATURES)
Creates an empty signal of length `n` samples, `ch` channels, and frame rate `fs`.
"""
signal(n::Int, ch::Int, fs) = signal(Array{Float64}(undef, n, ch), fs)

"""
$(SIGNATURES)
Creates an empty signal of type `T`, length `n` samples, and frame rate `fs`.
"""
signal(T::Type, n::Int, fs) = signal(Array{T}(undef, n), fs)

"""
$(SIGNATURES)
Creates an empty signal of type `T`, length `n` samples, `ch` channels, and frame rate `fs`.
"""
signal(T::Type, n::Int, ch::Int, fs) = signal(Array{T}(undef, n, ch), fs)

"""
$(SIGNATURES)
Creates a curried function that takes in an array and creates a signal with sampling rate `fs`.
"""
signal(fs) = x -> signal(x, fs)

"""
$(SIGNATURES)
Converts a signal to analytic representation.
"""
analytic(s::SampledSignal) = isanalytic(s) ? s : signal(hilbert(s)/√2.0, framerate(s))
analytic(s) = isanalytic(s) ? s : hilbert(s)/√2.0

"""
$(SIGNATURES)
Checks if a signal is analytic.
"""
isanalytic(s) = eltype(s) <: Complex

"""
$(SIGNATURES)
Gets the underlying samples in the signal.
"""
samples(s::SampledSignal) = getfield(s, :data)
samples(s) = s

"""
$(SIGNATURES)
Generates a padded view of a signal with optional delay/advance.
"""
function padded(s::AbstractVector{T}, padding; delay=0, fill=zero(T)) where {T}
  if length(padding) == 1
    left = padding
    right = padding
  else
    left = padding[1]
    right = padding[2]
  end
  PaddedView(fill, s, (1-left:length(s)+right,), (1+delay:delay+length(s),))
end

function padded(s::SampledSignal{T}, padding; delay=0, fill=zero(T)) where T
  signal(padded(samples(s), padding; delay=delay, fill=fill), framerate(s))
end

"""
    partition(x, n, step)

Iterates over the signal `x`, `n` samples at a time, with a step size of `step`.
"""
function Base.Iterators.partition(s::SampledSignalVector, n::Integer, step::Integer)
    n < 1 && throw(ArgumentError("cannot create partitions of length $n"))
    step < 1 && throw(ArgumentError("cannot create partitions with step size $step"))
    v = samples(s)
    return StepPartitionIterator{typeof(v)}(v, Int(n), Int(step))
end

## modified from Julia Iterators implementation

struct StepPartitionIterator{T <: AbstractVector}
  c::T
  n::Int
  step::Int
end

Base.eltype(::Type{StepPartitionIterator{T}}) where {T<:AbstractVector} = AbstractVector{eltype(T)}
Base.eltype(::Type{StepPartitionIterator{T}}) where {T<:Vector} = SubArray{eltype(T), 1, T, Tuple{UnitRange{Int}}, true}

Base.IteratorEltype(::Type{<:StepPartitionIterator{T}}) where {T<:AbstractVector} = Base.EltypeUnknown()
Base.IteratorEltype(::Type{<:StepPartitionIterator{T}}) where {T<:Vector} = Base.IteratorEltype(T)

partition_iteratorsize(::Base.HasShape) = Base.HasLength()
partition_iteratorsize(isz) = isz

function Base.IteratorSize(::Type{StepPartitionIterator{T}}) where {T}
  partition_iteratorsize(Base.IteratorSize(T))
end

function Base.length(itr::StepPartitionIterator)
  l = length(itr.c)
  return div(l, itr.step) + ((mod(l, itr.step) > 0) ? 1 : 0)
end

function Base.iterate(itr::StepPartitionIterator{<:AbstractRange}, state=1)
  state > length(itr.c) && return nothing
  r = min(state + itr.n - 1, length(itr.c))
  return @inbounds itr.c[state:r], state + itr.step
end

function Base.iterate(itr::StepPartitionIterator{<:AbstractVector}, state=1)
  state > length(itr.c) && return nothing
  r = min(state + itr.n - 1, length(itr.c))
  return @inbounds view(itr.c, state:r), state + itr.step
end

## end of Julia iterators implementation

"""
$(SIGNATURES)
Converts time to signal frame number.

# Examples:
```julia-repl
julia> x = signal(randn(2000), 8kHz);
julia> toframe(0.2s, x)
1601

julia> toframe([0.2s, 201ms], x)
2-element Array{Int64,1}:
 1601
 1609

julia> toframe(0.2:0.01:0.3, x)
 11-element Array{Int64,1}:
  1601
  1681
  1761
   ⋮
```
"""
toframe(t, s::SampledSignal) = 1 .+ round.(Int, inseconds.(t)*framerate(s))

"""
    (:)(start::Unitful.Time, stop::Unitful.Time)

Generates a time range index for a signal.

# Examples:
```julia-repl
julia> x = signal(randn(2000), 8kHz)
julia> x[0.2𝓈:0.201𝓈]
SampledSignal @ 8000.0 Hz, 9-element Array{Float64,1}:
 -0.08671384898800058
 -0.665143340284631
 -0.3955367460364236
  1.2386430598616671
 -0.4882254309443194
 -1.080437097803303
  0.8209785486953832
  1.3477512734963886
 -0.27722340584395494
```
"""
(::Colon)(t1::Units.Unitful.Time, t2::Units.Unitful.Time) = tuple(inseconds(t1), inseconds(t2))
(::Colon)(t1::Real, t2::Units.Unitful.Time) = tuple(inseconds(t1), inseconds(t2))

Base.getindex(s::SampledSignal, t::NTuple{2}) = Base.getindex(s, toframe(t[1], s):toframe(t[2], s))
Base.getindex(s::SampledSignal, t::NTuple{2}, ndx...) = Base.getindex(s, toframe(t[1], s):toframe(t[2], s), ndx...)
Base.setindex!(s::SampledSignal, v, t::NTuple{2}) = Base.setindex!(s, v, toframe(t[1], s):toframe(t[2], s))
Base.setindex!(s::SampledSignal, v, t::NTuple{2}, ndx...) = Base.setindex!(s, v, toframe(t[1], s):toframe(t[2], s), ndx...)
