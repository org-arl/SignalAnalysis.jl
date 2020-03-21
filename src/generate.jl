export cw, chirp
export hanning, hamming, tukey, cosine, lanczos, triang, bartlett, gaussian, bartlett_hann, blackman, kaiser, dpss

"""
$(SIGNATURES)
Generates a sinusoidal signal with specified `freq` and `duration` at frame
rate `fs`. The starting `phase` and `window` type may be optionally specified.

# Examples:
```julia-repl
julia> x = cw(5kHz, 200ms, 44.1kHz)
8821-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.20002267573696145s sampled at 44100.0Hz

julia> x = cw(5kHz, 200ms, 44.1kHz; window=hamming)
8821-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.20002267573696145s sampled at 44100.0Hz

julia> x = cw(-5kHz, 200ms, 44.1kHz; phase=45°, window=(tukey, 0.05))
8821-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.20002267573696145s sampled at 44100.0Hz
```
"""
function cw(freq, duration, fs; phase=0.0, window=nothing)
    abs(freq) > fs/2 && throw(ArgumentError("Frequency is beyond Nyquist"))
    t = 0:1/inHz(fs):inseconds(duration)
    x = cis.(2π*inHz(freq).*t .+ phase)
    if window != nothing
      x .*= getwindow(window, length(t))
    end
    signal(x, fs)
end

"""
$(SIGNATURES)
Generates a frequency modulated chirp signal from `freq1` to `freq2` and
specified `duration` at frame rate `fs`. The type of frequency modulation may
be controlled using `shape` (`:linear` (default) or `:hyperbolic`). The
starting `phase` and `window` type may be optionally specified.

# Examples:
```julia-repl
julia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz)
4411-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.10002267573696146s sampled at 44100.0Hz

julia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz; phase=45°, window=hamming)
4411-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.10002267573696146s sampled at 44100.0Hz

julia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz; shape=:hyperbolic, window=(tukey,0.05))
4411-frame, 1-channel SampleBuf{Complex{Float64}, 1}
0.10002267573696146s sampled at 44100.0Hz
```
"""
function chirp(freq1, freq2, duration, fs=2.0; shape=:linear, phase=0.0, window=nothing)
  f1 = inHz(freq1)
  f2 = inHz(freq2)
  d = inseconds(duration)
  t = 0:1/inHz(fs):d
  if shape == :linear
    cby2 = (f2-f1)/d/2.0
    x = cis.(2π .* (cby2.*t.^2 .+ f1.*t) .+ phase)
  elseif shape == :hyperbolic
    s = f2*d/(f2-f1)
    x = cis.(-2π*s*f1 .* log.(abs.(1.0 .- t./s)))
  else
    throw(ArgumentError("Unknown chirp shape: $shape"))
  end
  if window !== nothing
    x .*= getwindow(window, length(t))
  end
    signal(x, fs)
end

### utility functions

function getwindow(window, n)
  window === nothing && return ones(n)
  window isa Tuple && return window[1](n, window[2])
  window(n)
end
