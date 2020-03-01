using Plots

export psd, psd!, specgram, timeseries, timeseries!

"Plot power spectral density of the signal."
function psd(s; fs=deffs, nfft=1024, plot=plot, window=nothing, label=nothing, kwargs...)
  nfft = nextfastfft(nfft)
  while nfft > size(s,1)
    nfft = div(nfft, 2)
  end
  p = welch_pgram(s, nfft; fs=freqQ(fs), window=window)
  f = freq(p)
  funit = "Hz"
  if maximum(f) >= 10000
    f /= 1000.0
    funit = "kHz"
  end
  plot(f, power(p); label=label, kwargs...)
  xlabel!("Frequency ("*funit*")")
  ylabel!("Power spectral density (dB/Hz)")
end

"Plot power spectral density of the signal."
psd!(s; fs=deffs, nfft=1024, window=nothing, label=nothing, kwargs...) = psd(s; fs=fs, nfft=nfft, plot=plot!, window=window, label=label, kwargs...)

"Plot spectrogram of the signal."
function specgram(s; fs=deffs, nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, kwargs...)
  p = spectrogram(s, nfft, noverlap; fs=freqQ(fs), window=window)
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
  heatmap(t, f, z; kwargs...)
  xlabel!("Time ("*tunit*")")
  ylabel!("Frequency ("*funit*")")
end

"Plot timeseries of the signal."
function timeseries(s; fs=deffs, downsample=nothing, pooling=nothing, plot=plot, label=nothing, kwargs...)
  fs1 = freqQ(fs)
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
    fs1 /= downsample
  end
  t = time(s; fs=fs1)
  tunit = "s"
  if maximum(t) <= 1.0
    t *= 1000.0
    tunit = "ms"
  end
  plot(t, s; label=label, kwargs...)
  xlabel!("Time ("*tunit*")")
end

"Plot timeseries of the signal."
timeseries!(s; fs=deffs, downsample=nothing, pooling=nothing, label=nothing, kwargs...) = timeseries(s; fs=fs, downsample=downsample, pooling=pooling, plot=plot!, label=label, kwargs...)
