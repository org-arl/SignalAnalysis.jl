using Plots

export psd, psd!, specgram, timeseries, timeseries!, filtfreqz, filtfreqz!, orderedextrema

"Plot power spectral density of the signal."
function psd(s; fs=framerate(s), nfft=1024, plot=plot, window=nothing, legend=false, logfreq=false, kwargs...)
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
psd!(s; fs=framerate(s), nfft=1024, window=nothing, legend=false, kwargs...) = psd(s; fs=fs, nfft=nfft, plot=plot!, window=window, legend=legend, kwargs...)

"Plot spectrogram of the signal."
function specgram(s; fs=framerate(s), nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, t0=0.0, downsample=nothing, pooling=mean, kwargs...)
  @assert size(s,2) == 1 "specgram only works with vectors"
  p = spectrogram(s[:,1], nfft, noverlap; fs=freqQ(fs), window=window)
  t = t0 .+ time(p)
  f = freq(p)
  z = power(p)
  if downsample == nothing && length(t) > 5000
    downsample = ceil(Int, length(t)/5000)
    @warn "Too many points; downsampling by $downsample"
  end
  if downsample != nothing && downsample != 1
    if pooling == nothing
      phase = floor(Int, downsample/2)
      t = t[1+phase:downsample:end]
      z = z[:,1+phase:downsample:end]
    else
      y1 = collect(Iterators.flatten(pooling.(Periodograms.arraysplit(z[1,:], downsample, 0))))
      y = zeros(size(z,1), length(y1))
      y[1,:] .= y1
      for j in 2:size(z,1)
        y[j,:] .= collect(Iterators.flatten(pooling.(Periodograms.arraysplit(z[j,:], downsample, 0))))
      end
      z = y
      t = t[1:downsample:end][1:size(z,2)]
    end
  end
  if isanalytic(s)
    ndx = sortperm(f)
    f = f[ndx]
    z = z[ndx,:]
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

"Minmax pooling for perceptual integrity of timeseries."
orderedextrema(x) = argmin(x) < argmax(x) ? (minimum(x), maximum(x)) : (maximum(x), minimum(x))

"Plot timeseries of the signal."
function timeseries(s; fs=framerate(s), t0=0.0, downsample=nothing, pooling=orderedextrema, plot=plot, legend=false, kwargs...)
  fs = freqQ(fs)
  n = size(s,1)
  if isanalytic(s)
    s = real(s)
  end
  if downsample == nothing && n > 5000
    downsample = ceil(Int, n/5000)
    @warn "Too many points; downsampling by $downsample"
  end
  if downsample != nothing && downsample != 1
    if pooling == nothing
      phase = floor(Int, downsample/2)
      s = s[1+phase:downsample:end,:]
    elseif ndims(s) == 1
      s = collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s, downsample, 0))))
    else
      y = collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s[:,1], downsample, 0))))
      if size(s,2) > 1
        y1 = zeros(length(y), size(s,2))
        y1[:,1] .= y
        y = y1
      end
      for j in 2:size(s,2)
        y[:,j] .= collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s[:,j], downsample, 0))))
      end
      s = y
    end
    fs /= n/size(s,1)
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
timeseries!(s; fs=framerate(s), downsample=nothing, pooling=nothing, legend=false, kwargs...) = timeseries(s; fs=fs, downsample=downsample, pooling=pooling, plot=plot!, legend=legend, kwargs...)

"Plot frequency response of filter."
function filtfreqz(num::AbstractArray, den::AbstractArray=[1]; fs, nfreq=256, logfreq=false, plot=plot, legend=false, kwargs...)
  fs = freqQ(fs)
  f = collect(range(0, stop=fs/2, length=nfreq))
  z = Filters.freqz(PolynomialRatio(num, den), f, fs)
  funit = "Hz"
  if maximum(f) >= 10000
    f /= 1000.0
    funit = "kHz"
  end
  if logfreq
    p1 = plot(f[2:end], amp2db.(abs.(z[2:end])); leg=legend, xaxis=:log, kwargs...)
  else
    p1 = plot(f, amp2db.(abs.(z)); leg=legend, kwargs...)
  end
  xlabel!("Frequency ("*funit*")")
  ylabel!("Magnitude (dB)")
  if logfreq
    p2 = plot(f[2:end], angle.(z[2:end]); leg=legend, xaxis=:log, kwargs...)
  else
    p2 = plot(f, angle.(z); leg=legend, kwargs...)
  end
  xlabel!("Frequency ("*funit*")")
  ylabel!("Phase (radians)")
  plot(p1, p2, layout=2, size=(600,300))
end

"Plot frequency response of filter."
filtfreqz!(num::AbstractArray, den::AbstractArray=[1]; fs, nfreq=256, legend=false, kwargs...) = freqz(num, den; fs=fs, nfreq=nfreq, plot=plot!, legend=legend, kwargs...)
