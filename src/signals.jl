import SampledSignals
import SampledSignals: SampleBuf, domain

export signal, @rate, @samerateas, domain, analytic, isanalytic, samples
export padded, slide, toframe

SignalBase.nframes(x::SampleBuf) = SampledSignals.nframes(x)
SignalBase.framerate(x::SampleBuf) = SampledSignals.samplerate(x)
SignalBase.nchannels(x::SampleBuf) = SampledSignals.nchannels(x)
SignalBase.sampletype(x::SampleBuf) = eltype(x)

SignalBase.nframes(x::AbstractArray) = size(x,1)
SignalBase.framerate(x::AbstractArray) = 1.0
SignalBase.nchannels(x::AbstractArray) = size(x,2)
SignalBase.sampletype(x::AbstractArray) = eltype(x)

"""
$(SIGNATURES)
Creates a signal with frame rate `fs`.
"""
signal(x::AbstractArray, fs) = SampleBuf(x, inHz(fs))

"""
$(SIGNATURES)
Creates a signal with frame rate `fs`. If the original signal's frame rate is
the same as `fs`, this method simply returns the original signal. Otherwise, it
creates a new signal with the specified frame rate and data from the original
signal. Do note that this method does not resample the signal.
"""
signal(x::SampleBuf, fs) = fs == framerate(x) ? x : SampleBuf(x.data, inHz(fs))

"""
    @rate fs expr
Creates a signal from `expr` with frame rate `fs`. It provides syntactic sugar
on the [`signal`](@ref) method.

# Example:
```julia-repl
julia> x = @rate 44.1kHz randn(44100)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
```
"""
macro rate(fs, expr)
  :( signal($(esc(expr)), $(esc(fs))) )
end

"""
    @samerateas x expr
Creates a signal from `expr` with the same frame rate as signal `x`. This is
useful to preserve frame rate metadata across functions that do not return a
signal.

# Examples:
```julia-repl
julia> x = @rate 44.1kHz randn(44100)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz

julia> y = filt([1.0,0.5], x)     # frame rate stripped by DSP.filt
44100-element Array{Float64,1}:
   ⋮

julia> y = @samerateas x filt([1.0,0.5], x)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
```
"""
macro samerateas(x, expr)
  :( signal($(esc(expr)), framerate($(esc(x)))))
end

"""
    @samerateas n x expr
Creates a signal from `expr` with the a frame rate `n` times that of signal `x`.

# Examples:
```julia-repl
julia> x = @rate 44.1kHz randn(44100)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz

julia> y = @samerateas 1//2 x x[1:2:end]
22050-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 22050.0Hz

julia> y = @samerateas 2//3 x resample(x, 2//3)
29400-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 29400.0Hz
```
"""
macro samerateas(n, x, expr)
  :( signal($(esc(expr)), $(esc(n))*framerate($(esc(x)))))
end

"""
$(SIGNATURES)
Converts a signal to analytic representation.
"""
analytic(s::SampleBuf) = isanalytic(s) ? s : signal(hilbert(s), framerate(s))
analytic(s) = isanalytic(s) ? s : hilbert(s)

"""
$(SIGNATURES)
Checks if a signal is analytic.
"""
isanalytic(s) = eltype(s) <: Complex

"""
$(SIGNATURES)
Gets the underlying samples in the signal.
"""
samples(s::SampleBuf) = s.data
samples(s) = s

"""
$(SIGNATURES)
Generates a padded view of a signal with optional delay/advance.
"""
function padded(s::AbstractVector{T}, padding; delay=0, fill=zero(T)) where {T, N}
  if length(padding) == 1
    left = padding
    right = padding
  else
    left = padding[1]
    right = padding[2]
  end
  PaddedView(fill, s, (1-left:length(s)+right,), (1+delay:delay+length(s),))
end

"""
    slide(f::Function, s::AbstractVector, nframes, overlap=0, args...; showprogress=true)

Slides a window over a signal, processing each window. If the total number of frames
in the signal is not an integral multiple of `nframes`, the last incomplete
block of samples remains unprocessed.

The processing function receives a view on the original signal, and therefore
may modify the signal if desired.

# Examples:
```julia-repl
julia> x = signal(ones(1000), 8kHz);
julia> slide(x, 250) do x1, blknum, firstframe
         println(size(x1), ", ", blknum, ", ", firstframe)
       end
(250,), 1, 1
(250,), 2, 251
(250,), 3, 501
(250,), 4, 751

julia> slide(x, 250) do x1, blknum, firstframe
         x1 .= blknum
       end

julia> x[1], x[251], x[501], x[751]
(1.0, 2.0, 3.0, 4.0)
```
"""
function slide(f::Function, s::AbstractVector, nframes, overlap=0, args...; showprogress=true)
  @assert overlap < nframes "overlap must be less than nframes"
  n = size(s,1)
  m = nframes - overlap
  mmax = (n-nframes)÷m
  showprogress && (p = Progress(mmax+1, 1, "Processing: "))
  for j = 0:mmax
    s1 = @view s[j*m+1:j*m+nframes]
    f(s1, j+1, j*m+1, args...)
    showprogress && next!(p)
  end
end

"""
    slide(f::Function, ::Type{T}, s::AbstractVector, nframes, overlap=0, args...; showprogress=true) where T

Slides a window over a signal, processing each window, and collecting the results of type `T`.
If the total number of frames in the signal is not an integral multiple of
`nframes`, the last incomplete block of samples remains unprocessed.

# Examples:
```julia-repl
julia> x = signal(ones(1000), 8kHz);
julia> slide(Float32, x, 250) do x1, blknum, firstframe
         sum(x1)*blknum
       end
4-element Array{Float32,1}:
  250.0
  500.0
  750.0
 1000.0

julia> slide(Tuple{Int,Float64}, x, 250) do x1, blknum, firstframe
          (blknum, sum(x1)*blknum)
        end
4-element Array{Tuple{Int64,Float64},1}:
  (1, 250.0)
  (2, 500.0)
  (3, 750.0)
  (4, 1000.0)
```
"""
function slide(f::Function, ::Type{T}, s::AbstractVector, nframes, overlap=0, args...; showprogress=true) where {T}
  @assert overlap < nframes "overlap must be less than nframes"
  n = size(s,1)
  m = nframes - overlap
  mmax = (n-nframes)÷m
  out = Array{T,1}(undef, 1+mmax)
  showprogress && (p = Progress(mmax+1, 1, "Processing: "))
  for j = 0:mmax
    s1 = @view s[j*m+1:j*m+nframes]
    out[j+1] = f(s1, j+1, j*m+1, args...)
    showprogress && next!(p)
  end
  return out
end

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
toframe(t, s::SampleBuf) = 1 .+ round.(Int, inseconds.(t)*framerate(s))
