import SampledSignals
import SampledSignals: SampleBuf, domain

using DSP, DSP.Filters
using PaddedViews
using ProgressMeter

export signal, @rate, @samerateas, domain, analytic, isanalytic, samples
export padded, slide, toframe
export upconvert, downconvert, rrcosfir, rcosfir

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

julia> using DSP
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

julia> using DSP
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
analytic(s::SampleBuf) = isanalytic(s) ? s : signal(hilbert(s)/√2.0, framerate(s))
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

"""
$(SIGNATURES)
Converts baseband signal with `sps` symbols per passband sample to a real
passband signal centered around carrier frequency `fc`.
"""
function upconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  pad = cld(length(pulseshape), 2*sps) - 1
  s = vcat(zeros(pad), analytic(s), zeros(pad))
  s = signal(resample(s, sps, pulseshape), sps*fs)
  fc == 0 && (return s)
  √2 * real.(s .* cis.(2π * fc * domain(s)))
end

function upconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  out = Any[]
  for j in 1:nchannels(s)
    push!(out, upconvert(s[:,j], sps, fc, pulseshape; fs=fs))
  end
  hcat(out...)
end

"""
$(SIGNATURES)
Converts passband signal centered around carrier frequency `fc` to baseband,
and downsamples it by a factor of `sps`. If the `pulseshape` is specified to
be `nothing`, downsampling is performed without filtering.
"""
function downconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  s = signal(analytic(s), fs)
  s .*= cis.(-2π * fc * domain(s))
  sps == 1 && return signal(s, fs)
  pulseshape == nothing && return signal(s[1:sps:end,:], fs/sps)
  signal(resample(s, 1//sps, pulseshape), fs/sps)
end

function downconvert(s::AbstractMatrix, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  out = Any[]
  for j in 1:nchannels(s)
    push!(out, downconvert(s[:,j], sps, fc, pulseshape; fs=fs))
  end
  hcat(out...)
end

"""
$(SIGNATURES)
Root-raised cosine filter.
"""
function rrcosfir(β, sps, span = β < 0.68 ? 33-floor(Int, 44β) : 4)
  # default span based on http://www.commsys.isy.liu.se/TSKS04/lectures/3/MichaelZoltowski_SquareRootRaisedCosine.pdf
  delay = fld(span*sps, 2)
  t = collect(-delay:delay)/sps
  h = Array{Float64}(undef, size(t))
  for i ∈ 1:length(t)
    if t[i] == 0
      h[i] = (1 + β*(4/π - 1))/sps
    elseif abs(t[i]) == 1/(4β)
      h[i] = β/(√2*sps) * ((1+2/pi)*sin(π/(4β)) + (1-2/pi)*cos(π/(4β)))
    else
      h[i] = (sin(π*t[i]*(1-β)) + 4β*t[i]*cos(π*t[i]*(1+β))) / (π*t[i]*(1 - (4β*t[i])^2)) / sps
    end
  end
  h / √sum(h.^2)
end

"""
$(SIGNATURES)
Raised cosine filter.
"""
function rcosfir(β, sps, span = β < 0.68 ? 33-floor(Int, 44β) : 4)
  # default span based on http://www.commsys.isy.liu.se/TSKS04/lectures/3/MichaelZoltowski_SquareRootRaisedCosine.pdf
  # since the span is for rrcosfir, for rcosfir, it is very conservative
  delay = fld(span*sps, 2)
  t = collect(-delay:delay)/sps
  h = Array{Float64}(undef, size(t))
  for i ∈ 1:length(t)
    if abs(t[i]) == 1/(2β)
      h[i] = π/(4sps) * sinc(1/(2β))
    else
      h[i] = sinc(t[i]) * cos(π * β * t[i]) / (1-(2β * t[i])^2) / sps
    end
  end
  h / √sum(h.^2)
end
