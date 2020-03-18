export fir, removedc, removedc!, demon

"""
$(SIGNATURES)
Designs a `n`-tap FIR filter with a passband from `f1` to `f2` using the
specified `method`. If frame rate `fs` is not specified, `f1` and `f2` are given
in normalized units (1.0 being Nyquist). If `f1` is 0, the designed filter is
a lowpass filter, and if `f2` is `nothing` then it is a highpass filter.

This method is a convenience wrapper around [`DSP.digitalfilter`](@ref).

# Examples:
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units
julia> lpf = fir(127, 0, 10kHz; fs=44.1kHz)   # design a lowpass filter
127-element Array{Float64,1}:
[...]

julia> hpf = fir(127, 10kHz; fs=44.1kHz)      # design a highpass filter
127-element Array{Float64,1}:
[...]

julia> bpf = fir(127, 1kHz, 5kHz; fs=44.1kHz) # design a bandpass filter
127-element Array{Float64,1}:
[...]
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
  hpf = fir(127, cutoff; fs=fs)
  signal(filtfilt(hpf, y), fs)
end
