export isanalytic, analytic
export padded, slide, toframe
export energy, meantime, rmsduration, meanfrequency, rmsbandwidth, ifrequency

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
$(SIGNATURES)
Slides a window over a signal, processing each window. If the total number of frames
in the signal is not an integral multiple of `nframes`, the last incomplete
block of samples remains unprocessed.

The processing function receives a view on the original signal, and therefore
may modify the signal if desired.

# Examples:
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units
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
$(SIGNATURES)
Slides a window over a signal, processing each window, and collecting the results.
If the total number of frames in the signal is not an integral multiple of
`nframes`, the last incomplete block of samples remains unprocessed.

# Examples:
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units
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
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units
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
[...]
```
"""
toframe(t, s::SampleBuf) = 1 .+ round.(Int, inseconds.(t)*framerate(s))

"""
$(SIGNATURES)
Computes total signal energy.
"""
energy(s::AbstractVector; fs=framerate(s)) = sum(abs2, s)/inHz(fs)
energy(s::AbstractMatrix; fs=framerate(s)) = vec(sum(abs2, s; dims=1))./inHz(fs)

"""
$(SIGNATURES)
Computes mean time of the signal.
"""
meantime(s::SampleBuf) = wmean(domain(s), abs2.(samples(s)))
meantime(s; fs) = wmean((0:size(s,1)-1)/fs, abs2.(s))

"""
$(SIGNATURES)
Computes RMS duration of the signal.
"""
rmsduration(s::SampleBuf) = sqrt.(wmean(domain(s).^2, abs2.(samples(s))) .- meantime(s).^2)
rmsduration(s; fs) = sqrt.(wmean(((0:size(s,1)-1)/fs).^2, abs2.(s)) .- meantime(s; fs=fs).^2)

"""
$(SIGNATURES)
Computes instantaneous frequency of the signal.
"""
function ifrequency(s; fs=framerate(s))
  s1 = analytic(s)
  f1 = inHz(fs)/(2π) * diff(unwrap(angle.(s1); dims=1); dims=1)
  f2 = Array{Float32}(undef, size(s))
  f2[1,:] = f1[1,:]
  f2[2:end-1,:] .= (f1[1:end-1,:] .+ f1[2:end,:])./2
  f2[end,:] = f1[end,:]
  signal(f2, fs)
end

"""
$(SIGNATURES)
Computes mean frequency of a signal.
"""
function meanfrequency(s::AbstractMatrix; fs=framerate(s), nfft=1024, window=nothing)
  fs = inHz(fs)
  mapslices(samples(s); dims=1) do s1
    p = welch_pgram(s1, min(nfft, length(s1)); fs=fs, window=window)
    f = freq(p)
    wmean(f, power(p))
  end
end

function meanfrequency(s::AbstractVector; fs=framerate(s), nfft=1024, window=nothing)
  p = welch_pgram(s, min(nfft, length(s)); fs=inHz(fs), window=window)
  f = freq(p)
  wmean(f, power(p))
end

"""
$(SIGNATURES)
Computes RMS bandwidth of a signal.
"""
function rmsbandwidth(s::AbstractMatrix; fs=framerate(s), nfft=1024, window=nothing) where T
  fs = inHz(fs)
  mapslices(samples(s); dims=1) do s1
    p = welch_pgram(s1, min(nfft, length(s1)); fs=fs, window=window)
    f = freq(p)
    f0 = wmean(f, power(p))
    sqrt.(wmean((f.-f0).^2, power(p)))
  end
end

function rmsbandwidth(s::AbstractVector; fs=framerate(s), nfft=1024, window=nothing) where T
  p = welch_pgram(s, min(nfft, length(s)); fs=inHz(fs), window=window)
  f = freq(p)
  f0 = wmean(f, power(p))
  sqrt.(wmean((f.-f0).^2, power(p)))
end

### utility functions

wmean(x::AbstractVector, w::AbstractVector) = (x'w) / sum(w)
wmean(x, w::AbstractVector) = sum(x.*w) ./ sum(w)
wmean(x, w::AbstractMatrix) = sum(x.*w; dims=1) ./ sum(w; dims=1)
