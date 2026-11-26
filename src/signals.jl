using MetaArrays: MetaArray
using Base.Iterators: partition
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

function Base.show(io::IO, mime::MIME"text/plain", s::SampledSignal{T}) where T
  print(io, "SampledSignal @ ", framerate(s), " Hz, ")
  show(io, mime, samples(s))
end

"""
$(SIGNATURES)
Returns the domain of the signal.
"""
domain(x) = (0:nframes(x)-1) ./ framerate(x)

"""
    time(i, x::SampledSignal)

Gets the time of the `i`th sample in the signal. The index `i` can be a range
or an array of indices.
"""
Base.time(i::Real, x::SampledSignal) = (i - 1) / framerate(x)
Base.time(i::AbstractArray{<:Real}, x::SampledSignal) = (i .- 1) ./ framerate(x)

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
  ext = Base.get_extension(@__MODULE__, :WAVExt)
  ext === nothing && error("Loading signals from WAV files requires the WAV package to be loaded; run `using WAV` first")
  ext.loadsignal(filename; start, nsamples)
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
$(TYPEDSIGNATURES)
Creates a signal by sampling function `f` at frame rate `fs`, from time `t0` up
to, but not including, time `t1`. Function `f` is called with time in seconds.
The time axis of the resulting signal starts at 0.

# Examples:
```julia-repl
julia> x = signal(t -> sin(2π * 100 * t), 1𝓈, 2𝓈, 8kHz)
julia> x = signal(sinpi, 0, 2, 100)
```
"""
function signal(f::Function, t0, t1, fs)
  t = (inseconds(t0):1/inHz(fs):inseconds(t1))[1:end-1]
  signal(f.(t), fs)
end

"""
$(TYPEDSIGNATURES)
Creates a signal by sampling function `f` at frame rate `fs`, from time 0 up
to, but not including, time `T`. Function `f` is called with time in seconds.

# Examples:
```julia-repl
julia> x = signal(t -> sin(2π * 100 * t), 500ms, 8kHz)
```
"""
signal(f::Function, T, fs) = signal(f, 0, T, fs)

"""
$(SIGNATURES)
Converts a signal to analytic representation. The conversion preserves energy,
i.e., to convert back to a real signal while conserving energy, multiply by √2.
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
samples(s::SampledSignal) = samples(getfield(s, :data))
samples(s::PaddedView) = collect(s)
samples(s) = s

const PaddedIndex = Union{Integer,AbstractRange{Int},Vector{Int}}

struct InfinitePaddedView{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  data::A
  fill::T
  delay::Int
end

Base.size(p::InfinitePaddedView) = size(p.data)
Base.axes(p::InfinitePaddedView) = axes(p.data)

Base.checkbounds(::Type{Bool}, p::InfinitePaddedView, i::PaddedIndex, j::PaddedIndex...) =
  Base.checkbounds_indices(Bool, Base.tail(axes(p.data)), j)

function Base.getindex(p::InfinitePaddedView{T,1}, i::Int) where T
  k = i - p.delay
  checkbounds(Bool, p.data, k) ? p.data[k] : p.fill
end

function Base.getindex(p::InfinitePaddedView{T,2}, i::Int, j::Int) where T
  Base.checkindex(Bool, axes(p.data, 2), j) || throw(BoundsError(p, (i, j)))
  k = i - p.delay
  checkbounds(Bool, p.data, k, j) ? p.data[k, j] : p.fill
end

samples(p::InfinitePaddedView) = collect(p)

"""
    padded(s, padding; delay=0, fill=0)
    padded(s; delay=0, fill=0)

Pads a signal with `fill`, with optional delay/advance. A positive `delay` shifts
the signal later, so that `padded(s; delay=k)[i] == s[i - k]`.

Without `padding`, the signal is padded indefinitely in both directions. The result
is a view with the same indices and length as `s`, but any index is valid. This is
the simplest way to access a signal beyond its ends, e.g. `y[i-k]` in a loop.

With `padding` (a number, or a tuple of left and right amounts), the signal is
extended by that many samples:
- For an array, the result is a view that keeps the original indexing, so its
  indices start at `1 - left`.
- For a `SampledSignal`, the result is a new signal indexed from 1, with the
  padding prepended and appended. Padding on the left therefore delays the
  signal: its first sample moves to index `left + 1`, at time `left / fs`.
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

padded(s::AbstractVecOrMat{T}; delay=0, fill=zero(T)) where {T} = InfinitePaddedView(s, convert(T, fill), Int(delay))

padded(s::SampledSignal{T}, padding; delay=0, fill=zero(T)) where T =
  samerateas(s, padded(samples(s), padding; delay=delay, fill=fill))

padded(s::SampledSignal{T}; delay=0, fill=zero(T)) where T =
  MetaArray(getfield(s, :meta), padded(getfield(s, :data); delay=delay, fill=fill))

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
    toframe(t, s::SampledSignal)
    toframe(t, fs)

Converts time to signal frame number.

# Examples:
```julia-repl
julia> x = signal(randn(2000), 8kHz);
julia> toframe(0.2𝓈, x)
1601

julia> toframe(0.2𝓈, 8kHz)
1601

julia> toframe([0.2𝓈, 0.201𝓈], x)
2-element Array{Int64,1}:
  1601
  1609

julia> toframe(0.2:0.201𝓈, x)
1601:1609

julia> toframe((0.2, 0.201), x)
1601:1609

julia> toframe((0.2, 0.201), 8000)
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
toframe(t, fs) = 1 .+ round.(Int, inseconds.(t) .* inHz(fs))
toframe(t::Tuple{Any,Any}, s::SampledSignal) = UnitRange((1 .+ round.(Int, inseconds.(t) .* framerate(s)))...)
toframe(t::Tuple{Any,Any}, fs) = UnitRange((1 .+ round.(Int, inseconds.(t) .* inHz(fs)))...)

"""
    (:)(start::Unitful.Time, stop::Unitful.Time)
    (:)(start::Unitful.Frequency, stop::Unitful.Frequency)

Generates a time range index for a signal, or a frequency range index for a
spectrum. If `start` has no units, it takes the units of `stop`, so `200:300ms`
is 200 ms to 300 ms, and `1:2kHz` is 1 kHz to 2 kHz.

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
(::Colon)(t1::Units.Unitful.Time, t2::Units.Unitful.Time) = (t1, t2)
(::Colon)(t1::Real, t2::Units.Unitful.Time) = (t1 * Units.Unitful.unit(t2), t2)
(::Colon)(f1::Units.Unitful.Frequency, f2::Units.Unitful.Frequency) = (f1, f2)
(::Colon)(f1::Real, f2::Units.Unitful.Frequency) = (f1 * Units.Unitful.unit(f2), f2)

const TimeRange = Tuple{Units.Unitful.Time,Units.Unitful.Time}

Base.getindex(s::SampledSignal, t::TimeRange, ndx...) = Base.getindex(s, toframe(t[1], s):toframe(t[2], s), ndx...)
Base.setindex!(s::SampledSignal, v, t::TimeRange, ndx...) = Base.setindex!(s, v, toframe(t[1], s):toframe(t[2], s), ndx...)

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
