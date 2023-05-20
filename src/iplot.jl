using .InteractiveViz
using Statistics
using DSP.Periodograms

export iplot, iplot!, ispecgram

"""
$(SIGNATURES)
Plots interactive timeseries of the signal.
"""
function iplot(s::SampledSignal, args...; kwargs...)
  s1 = samples(s)
  if isanalytic(s1)
    @warn "Plotting only real part of complex signal"
    s1 = real.(s1)
  end
  t = domain(s)
  if maximum(t) <= 1.0
    t *= 1000.0
    xlabel = "Time (ms)"
  else
    xlabel = "Time (s)"
  end
  ilines(t, s1, args...; axis=(; xlabel=xlabel), kwargs...)
end

"""
$(SIGNATURES)
Plots interactive timeseries of the signal over a previous plot.
"""
function iplot!(s::SampledSignal; kwargs...)
  s1 = samples(s)
  if isanalytic(s1)
    @warn "Plotting only real part of complex signal"
    s1 = real.(s1)
  end
  t = domain(s)
  maximum(t) <= 1.0 && (t *= 1000.0)
  ilines!(t, s1; kwargs...)
end

"""
$(SIGNATURES)
Plots interactive spectrogram of the signal.
"""
function ispecgram(s; fs=framerate(s), nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, kwargs...)
  @assert Base.size(s,2) == 1 "ispecgram only works with vectors"
  s1 = vec(samples(s))
  if isanalytic(s1)
    @warn "Using only real part of complex signal"
    s1 = real.(s1)
  end
  p = spectrogram(s1, nfft, noverlap; fs=inHz(fs), window=window)
  t = time(p)
  f = freq(p)
  tunit = "s"
  funit = "Hz"
  if maximum(f) >= 10000
    f /= 1000.0
    funit = "kHz"
  end
  if maximum(t) <= 1.0
    t *= 1000.0
    tunit = "ms"
  end
  xlabel = "Time ("*tunit*")"
  ylabel = "Frequency ("*funit*")"
  z = pow2db.(power(p)')
  cmax = ceil(Int, maximum(z)/5)*5
  cmin = floor(Int, max(cmax-50, minimum(z))/5)*5
  iheatmap(range(first(t), last(t); length=size(z,1)), range(first(f), last(f); length=size(z,2)), z;
    colorrange=(cmin, cmax),
    axis=(; xlabel=xlabel, ylabel=ylabel), kwargs...)
end
