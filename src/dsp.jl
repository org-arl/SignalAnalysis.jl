export fir, removedc, removedc!, demon
export upconvert, downconvert, rrcosfir, rcosfir

"""
$(SIGNATURES)
Designs a `n`-tap FIR filter with a passband from `f1` to `f2` using the
specified `method`. If frame rate `fs` is not specified, `f1` and `f2` are given
in normalized units (1.0 being Nyquist). If `f1` is 0, the designed filter is
a lowpass filter, and if `f2` is `nothing` then it is a highpass filter.

This method is a convenience wrapper around `DSP.digitalfilter`.

# Examples:
```julia-repl
julia> lpf = fir(127, 0, 10kHz; fs=44.1kHz)   # design a lowpass filter
127-element Array{Float64,1}:
  ⋮

julia> hpf = fir(127, 10kHz; fs=44.1kHz)      # design a highpass filter
127-element Array{Float64,1}:
  ⋮

julia> bpf = fir(127, 1kHz, 5kHz; fs=44.1kHz) # design a bandpass filter
127-element Array{Float64,1}:
  ⋮
```
"""
function fir(n, f1, f2=nothing; fs=2.0, method=FIRWindow(hanning(n)))
  fs = inHz(fs)
  if f1 == 0
    f = Lowpass(inHz(f2); fs=fs)
  elseif f2 == nothing || inHz(f2) == fs/2
    f = Highpass(inHz(f1); fs=fs)
  else
    f = Bandpass(inHz(f1), inHz(f2); fs=fs)
  end
  return digitalfilter(f, method)
end

"""
$(SIGNATURES)
DC removal filter. Parameter `α` controls the cutoff frequency. Implementation
based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).

See also: [`removedc`](@ref)
"""
function removedc!(s; α=0.95)
  for k = 1:size(s,2)
    for j = 2:size(s,1)
      s[j,k] += α*s[j-1,k]
    end
    s[2:end,k] .+= -s[1:end-1,k]
  end
  s *= sqrt(α)
  return s
end

"""
$(SIGNATURES)
DC removal filter. Parameter `α` controls the cutoff frequency. Implementation
based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).

See also: [`removedc!`](@ref)
"""
removedc(s; α=0.95) = removedc!(copy(s); α=α)

"""
$(SIGNATURES)
Estimates DEMON spectrum. The output is highpass filtered with a `cutoff`
frequency and downsampled. Supported downsampling methods are `:rms` (default),
`:mean` and `:fir`.
"""
function demon(x; fs=framerate(x), downsample=250, method=:rms, cutoff=1.0)
  local y
  fs /= downsample
  for k = 1:size(x,2)
    if downsample == 1
      y1 = abs.(hilbert(x[:,k]))
    elseif method == :rms
      y1 = sqrt.(mean.(Periodograms.arraysplit(abs2.(hilbert(x[:,k])), downsample, 0)))
    elseif method == :mean
      y1 = mean.(Periodograms.arraysplit(abs.(hilbert(x[:,k])), downsample, 0))
    elseif method == :fir
      aaf = fir(127, 0, 0.48fs; fs=fs)
      y1 = filtfilt(aaf, abs.(hilbert(x[:,k])))[1:downsample:end]
    else
      throw(ArgumentError("Unknown method"))
    end
    if k == 1
      y = zeros(length(y1), size(x,2))
    end
    y[:,k] .= y1
  end
  maxflen = length(y)÷4
  mod(maxflen, 2) == 0 && (maxflen += 1)
  hpf = fir(min(127, maxflen), cutoff; fs=fs)
  signal(filtfilt(hpf, y), fs)
end

function demon(x::AbstractVector{T}; fs=framerate(x), downsample=250, method=:rms, cutoff=1.0) where T
  y = @view samples(x)[:,1:1]
  z = demon(y, fs=fs, downsample=downsample, method=method, cutoff=cutoff)
  @samerateas z dropdims(samples(z), dims=2)
end

"""
$(SIGNATURES)
Converts baseband signal with `sps` symbols per passband sample to a real
passband signal centered around carrier frequency `fc`.
"""
function upconvert(s::AbstractVector, sps, fc, pulseshape=rrcosfir(0.25, sps); fs=framerate(s))
  pad = cld(length(pulseshape), 2*sps) - 1
  s = vcat(zeros(pad), analytic(s), zeros(pad))
  s = signal(resample(s, sps, pulseshape), sps*fs)
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
