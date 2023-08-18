using MetaArrays: MetaArray
using Base.Iterators: partition
using WAV: wavread
using DSP: hilbert
using PaddedViews: PaddedView

export signal, analytic, isanalytic, samples
export padded, toframe, domain
export partition, samerateas

struct SamplingInfo
  fs::Float32
end

Base.keys(x::SamplingInfo) = fieldnames(typeof(x))

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
$(TYPEDSIGNATURES)
Creates a signal with frame rate `fs`.
"""
signal(x::AbstractArray, fs) = MetaArray(SamplingInfo(inHz(fs)), x)
signal(x::Base.Iterators.Flatten, fs) = signal(collect(x), fs)

"""
$(TYPEDSIGNATURES)
Creates a signal with frame rate `fs`. If the original signal's frame rate is
the same as `fs`, this method simply returns the original signal. Otherwise, it
creates a new signal with the specified frame rate and data from the original
signal. Do note that this method does not resample the signal.
"""
signal(x::SampledSignal, fs) = fs == framerate(x) ? x : signal(samples(x), fs)

"""
$(TYPEDSIGNATURES)
Loads a signal from a WAV file.
"""
function signal(filename::AbstractString; start=1, nsamples=missing)
  if start == 1 && nsamples === missing
    data, fs = wavread(filename)
  elseif start == 1
    data, fs = wavread(filename; subrange=nsamples)
  else
    if nsamples === missing
      n, ch = wavread(filename; format="size")
      nsamples = n - start + 1
    end
    data, fs = wavread(filename; subrange=start:(start+nsamples-1))
  end
  signal(data, fs)
end

"""
$(TYPEDSIGNATURES)
Creates an empty signal of length `n` samples, and frame rate `fs`.
"""
signal(n::Int, fs) = signal(Array{Float64}(undef, n), fs)

"""
$(TYPEDSIGNATURES)
Creates an empty signal of length `n` samples, `ch` channels, and frame rate `fs`.
"""
signal(n::Int, ch::Int, fs) = signal(Array{Float64}(undef, n, ch), fs)

"""
$(TYPEDSIGNATURES)
Creates an empty signal of type `T`, length `n` samples, and frame rate `fs`.
"""
signal(T::Type, n::Int, fs) = signal(Array{T}(undef, n), fs)

"""
$(TYPEDSIGNATURES)
Creates an empty signal of type `T`, length `n` samples, `ch` channels, and frame rate `fs`.
"""
signal(T::Type, n::Int, ch::Int, fs) = signal(Array{T}(undef, n, ch), fs)

"""
$(TYPEDSIGNATURES)
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
$(TYPEDSIGNATURES)
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
  PaddedView(fill, s, (firstindex(s,1)-left:lastindex(s,1)+right,), (firstindex(s,1)+delay:delay+lastindex(s,1),))
end

"""
$(TYPEDSIGNATURES)
Generates a padded view of a signal with optional delay/advance.
"""
function padded(s::AbstractMatrix{T}, padding; delay=0, fill=zero(T)) where {T}
  if length(padding) == 1
    left = padding
    right = padding
  else
    left = padding[1]
    right = padding[2]
  end
  channelindices = firstindex(s,2):lastindex(s,2)
  PaddedView(fill, s, (firstindex(s,1)-left:lastindex(s,1)+right,channelindices), (firstindex(s,1)+delay:delay+lastindex(s,1),channelindices))
end

function padded(s::SampledSignal{T}, padding; delay=0, fill=zero(T)) where T
  samerateas(s, padded(samples(s), padding; delay=delay, fill=fill))
end

"""
    partition(x::SampledSignal, n; step=n, flush=true)

Iterates over the signal `x`, `n` samples at a time, with a step size of `step`. If `flush` is
enabled, the last partition may be smaller than `n` samples.

When applied to a multichannel signal `x`, each partition contains samples from all channels.

# Examples:
```julia-repl
julia> x = signal(collect(1:10), 1.0);
julia> collect(partition(x, 5))
2-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:
 [1, 2, 3, 4, 5]
 [6, 7, 8, 9, 10]

julia> collect(partition(x, 5; step=2))
5-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:
 [1, 2, 3, 4, 5]
 [3, 4, 5, 6, 7]
 [5, 6, 7, 8, 9]
 [7, 8, 9, 10]
 [9, 10]

julia> collect(partition(x, 5; step=2, flush=false))
3-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:
 [1, 2, 3, 4, 5]
 [3, 4, 5, 6, 7]
 [5, 6, 7, 8, 9]

julia> x = signal(hcat(collect(1:10), collect(11:20)), 1.0);
julia> collect(partition(x, 5))[1]
5×2 view(::Array{Int64,2}, 1:5, :) with eltype Int64:
 1  11
 2  12
 3  13
 4  14
 5  15
```
"""
function Base.Iterators.partition(s::SampledSignal, n::Integer; step::Integer=n, flush::Bool=true)
  n < 1 && throw(ArgumentError("cannot create partitions of length $n"))
  step < 1 && throw(ArgumentError("cannot create partitions with step size $step"))
  SignalPartitionIterator{typeof(s)}(s, Int(n), Int(step), flush)
end

struct SignalPartitionIterator{T <: AbstractArray}
  c::T
  n::Int
  step::Int
  flush::Bool
end

Base.IteratorEltype(::Type{SignalPartitionIterator}) = Base.HasEltype()
Base.IteratorSize(::Type{SignalPartitionIterator}) = Base.HasLength()

Base.eltype(::Type{SignalPartitionIterator{T}}) where {T<:AbstractVector} = AbstractVector{eltype(T)}
Base.eltype(::Type{SignalPartitionIterator{T}}) where {T<:AbstractMatrix} = AbstractMatrix{eltype(T)}
Base.eltype(::Type{SignalPartitionIterator{T}}) where {T<:Vector} = SubArray{eltype(T), 1, T, Tuple{UnitRange{Int}}, true}
Base.eltype(::Type{SignalPartitionIterator{T}}) where {T<:Matrix} = SubArray{eltype(T), 2, T, Tuple{UnitRange{Int},Base.Slice{Base.OneTo{Int64}}}, false}

function Base.length(itr::SignalPartitionIterator)
  l = size(itr.c, 1)
  itr.flush || (l -= itr.n-1)
  div(l, itr.step) + ((mod(l, itr.step) > 0) ? 1 : 0)
end

function Base.iterate(itr::SignalPartitionIterator{<:AbstractRange}, state=firstindex(itr.c))
  l = lastindex(itr.c)
  state > l && return nothing
  itr.flush || state + itr.n - 1 <= l || return nothing
  r = min(state + itr.n - 1, l)
  @inbounds itr.c[state:r], state + itr.step
end

function Base.iterate(itr::SignalPartitionIterator{<:AbstractVector}, state=firstindex(itr.c))
  l = lastindex(itr.c)
  state > l && return nothing
  itr.flush || state + itr.n - 1 <= l || return nothing
  r = min(state + itr.n - 1, l)
  @inbounds view(itr.c, state:r), state + itr.step
end

function Base.iterate(itr::SignalPartitionIterator{<:AbstractMatrix}, state=firstindex(itr.c))
  l = lastindex(itr.c, 1)
  state > l && return nothing
  itr.flush || state + itr.n - 1 <= l || return nothing
  r = min(state + itr.n - 1, l)
  @inbounds view(itr.c, state:r, :), state + itr.step
end

"""
$(SIGNATURES)
Converts time to signal frame number.

# Examples:
```julia-repl
julia> x = signal(randn(2000), 8kHz);
julia> toframe(0.2𝓈, x)
1601

julia> toframe([0.2𝓈, 0.201𝓈], x)
2-element Array{Int64,1}:
  1601
  1609

julia> toframe(0.2:0.201𝓈, x)
1601:1609

julia> toframe((0.2, 0.201), x)
1601:1609

julia> toframe(0.2:0.01:0.3, x)
11-element Array{Int64,1}:
  1601
  1681
  1761
   ⋮
```
"""
toframe(t, s::SampledSignal) = 1 .+ round.(Int, inseconds.(t) .* framerate(s))
toframe(t::NTuple{2,<:Real}, s::SampledSignal) = UnitRange((1 .+ round.(Int, inseconds.(t) .* framerate(s)))...)

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

function Base.vec(s::SampledSignal)
  if ndims(s) < 3 && isone(size(s, 2))
    signal(vec(samples(s)), framerate(s))
  else
    throw(ArgumentError("reshape a multi-channel signal as a single-channel signal is undefined."))
  end
end

"""
$(SIGNATURES)
Checks if two signals have the same sampling rate. If the sampling rate of a
signal is unknown (because it is not a `SampledSignal`), it is assumed to have
the same rate as the other signal.
"""
issamerate(x::SampledSignal, y::SampledSignal) = framerate(x) == framerate(y)
issamerate(x::SampledSignal, y::AbstractArray) = true
issamerate(x::AbstractArray, y::SampledSignal) = true
issamerate(x::AbstractArray, y::AbstractArray) = true

"""
$(SIGNATURES)
Create a signal with the same sampling rate as signal `x`.

# Examples:
```julia-repl
julia> x = signal(randn(100), 8kHz)
julia> y = samerateas(x)(randn(5))
SampledSignal @ 8000.0 Hz, 5-element Array{Float64,1}:
 -0.08671384898800058
 -0.665143340284631
 -0.3955367460364236
 0.8209785486953832
 1.3477512734963886
```
"""
samerateas(x) = y -> signal(samples(y), framerate(x))

"""
$(SIGNATURES)
Create a signal with samples `y` and sampling rate same as signal `x`.

# Examples:
```julia-repl
julia> x = signal(randn(100), 8kHz)
julia> y = samerateas(x, randn(5))
SampledSignal @ 8000.0 Hz, 5-element Vector{Float64}:
 -0.3053704876108388
 -0.5474123820044299
 -0.6916442204609657
 -0.5185296405433826
 -0.4598263144701988
```
"""
samerateas(x, y) = signal(samples(y), framerate(x))
