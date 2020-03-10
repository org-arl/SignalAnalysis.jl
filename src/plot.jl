using Plots

export psd, psd!, specgram, timeseries, timeseries!, filtfreqz, filtfreqz!

"Plot power spectral density of the signal."
function psd(s; fs=deffs[], nfft=1024, plot=plot, window=nothing, legend=false, logfreq=false, kwargs...)
  nfft = nextfastfft(nfft)
  while nfft > size(s,1)
    nfft = div(nfft, 2)
  end
  funit = "Hz"
  for j = 1:size(s,2)
    p = welch_pgram(s[:,j], nfft; fs=freqQ(fs), window=window)
    f = freq(p)
    z = power(p)
    if isanalytic(s)
      ndx = sortperm(f)
      f = f[ndx]
      z = z[ndx]
    end
    if maximum(f) >= 10000
      f /= 1000.0
      funit = "kHz"
    end
    if logfreq
      plot(f[2:end], pow2db.(z[2:end]); leg=legend, xaxis=:log, kwargs...)
    else
      plot(f, pow2db.(z); leg=legend, kwargs...)
    end
    plot = plot!
  end
  xlabel!("Frequency ("*funit*")")
  ylabel!("Power spectral density (dB/Hz)")
end

"Plot power spectral density of the signal."
psd!(s; fs=deffs[], nfft=1024, window=nothing, legend=false, kwargs...) = psd(s; fs=fs, nfft=nfft, plot=plot!, window=window, legend=legend, kwargs...)

"Plot spectrogram of the signal."
function specgram(s; fs=deffs[], nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, kwargs...)
  @assert size(s,2) == 1 "specgram only works with vectors"
  p = spectrogram(s[:,1], nfft, noverlap; fs=freqQ(fs), window=window)
  t = time(p)
  f = freq(p)
  z = power(p)
  if isanalytic(s)
    ndx = sortperm(f)
    f = f[ndx]
    z = z[ndx, :]
  end
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
  heatmap(t, f, pow2db.(z); kwargs...)
  xlabel!("Time ("*tunit*")")
  ylabel!("Frequency ("*funit*")")
end

"Plot timeseries of the signal."
function timeseries(s; fs=deffs[], t0=0.0, downsample=nothing, pooling=nothing, plot=plot, legend=false, kwargs...)
  fs = freqQ(fs)
  n = size(s,1)
  if isanalytic(s)
    s = real(s)
  end
  if downsample == nothing && n > 10000
    downsample = ceil(Int, n/10000)
    @warn "Too many points; downsampling by $downsample"
  end
  if downsample != nothing && downsample != 1
    if pooling == nothing
      phase = floor(Int, downsample/2)
      s = s[1+phase:downsample:end,:]
    elseif ndims(s) == 1
      s = pooling.(Periodograms.arraysplit(s, downsample, 0))
    else
      y = zeros(floor(Int,n/downsample), size(s,2))
      for j in 1:size(s,2)
        y[:,j] = pooling.(Periodograms.arraysplit(s[:,j], downsample, 0))
      end
      s = y
    end
    fs /= downsample
  end
  t = time(s; t0=t0, fs=fs)
  tunit = "s"
  if maximum(t) <= 1.0
    t *= 1000.0
    tunit = "ms"
  end
  plot(t, s; leg=legend, kwargs...)
  xlims!(minimum(t), maximum(t))
  xlabel!("Time ("*tunit*")")
end

"Plot timeseries of the signal."
timeseries!(s; fs=deffs[], downsample=nothing, pooling=nothing, legend=false, kwargs...) = timeseries(s; fs=fs, downsample=downsample, pooling=pooling, plot=plot!, legend=legend, kwargs...)

"Plot frequency response of filter."
function filtfreqz(num::AbstractArray, den::AbstractArray=[1]; fs=deffs[], nfreq=256, plot=plot, legend=false, kwargs...)
  fs = freqQ(fs)
  f = collect(range(0, stop=fs/2, length=nfreq))
  z = Filters.freqz(PolynomialRatio(num, den), f, fs)
  funit = "Hz"
  if maximum(f) >= 10000
    f /= 1000.0
    funit = "kHz"
  end
  p1 = plot(f, pow2db.(abs.(z)); leg=legend, kwargs...)
  xlabel!("Frequency ("*funit*")")
  ylabel!("Magnitude (dB)")
  p2 = plot(f, angle.(z); leg=legend, kwargs...)
  xlabel!("Frequency ("*funit*")")
  ylabel!("Phase (radians)")
  plot(p1, p2, layout=2, size=(600,300))
end

"Plot frequency response of filter."
filtfreqz!(num::AbstractArray, den::AbstractArray=[1]; fs=deffs[], nfreq=256, legend=false, kwargs...) = freqz(num, den; fs=fs, nfreq=nfreq, plot=plot!, legend=legend, kwargs...)
