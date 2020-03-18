import SampledSignals
import SampledSignals: SampleBuf, domain

export signal, @rate, @samerateas, domain, analytic, isanalytic, samples

SignalBase.nframes(x::SampleBuf) = SampledSignals.nframes(x)
SignalBase.framerate(x::SampleBuf) = SampledSignals.samplerate(x)
SignalBase.nchannels(x::SampleBuf) = SampledSignals.nchannels(x)
SignalBase.sampletype(x::SampleBuf) = eltype(x)

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
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units
julia> x = @rate 44.1kHz randn(44100)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
[...]
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
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units, DSP
julia> x = @rate 44.1kHz randn(44100);
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
[...]

julia> y = filt([1.0,0.5], x)     # frame rate stripped by filt
44100-element Array{Float64,1}:
[...]

julia> y = @samerateas x filt([1.0,0.5], x)
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
[...]
```
"""
macro samerateas(x, expr)
  :( signal($(esc(expr)), framerate($(esc(x)))))
end

"""
    @samerateas n x expr
Creates a signal from `expr` with the a frame rate `n` times that of signal `x`.

# Examples:
```jldoctest
julia> using SignalAnalysis, SignalAnalysis.Units, DSP
julia> x = @rate 44.1kHz randn(44100);
44100-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 44100.0Hz
[...]

julia> y = @samerateas 1//2 x x[1:2:end]
22050-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 22050.0Hz
[...]

julia> y = @samerateas 2//3 x resample(x, 2//3)
29400-frame, 1-channel SampleBuf{Float64, 1}
1.0s sampled at 29400.0Hz
[...]
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
