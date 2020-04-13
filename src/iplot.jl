using InteractiveViz
using Statistics
using DSP.Periodograms

export iplot, ispecgram

"""
$(SIGNATURES)
Plots interactive timeseries of the signal.
"""
function InteractiveViz.iplot(s::SampleBuf, args...; kwargs...)
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
    xlabel = "Time (seconds)"
  end
  iplot(t, s1, args...; xlabel=xlabel, kwargs...)
end

"""
$(SIGNATURES)
Plots interactive spectrogram of the signal.
"""
function ispecgram(s; fs=framerate(s), nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, pooling=mean, kwargs...)
  @assert Base.size(s,2) == 1 "ispecgram only works with vectors"
  p = spectrogram(samples(s), nfft, noverlap; fs=inHz(fs), window=window)
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
  cmin = max(cmax-50, floor(Int, minimum(z)/5)*5)
  iheatmap(z, t[1], t[end], f[1], f[end]; clim=(cmin, cmax), xlabel=xlabel, ylabel=ylabel, pooling=pooling, kwargs...)
end
