using Plots

"""
    plot(data; kwargs...)

Plots timeseries from a sample buffer.

# Optional keyword arguments
- `t0=0.0`: start time
- `downsample=:auto`: downsampling factor (integer)
- `pooling=:auto`: pooling mode (`:min`, `:max`, `:mean`, `:minmax`, `nothing` or function)
- other `kwargs` are passed on to `plot`
"""
@recipe function plot(s::SampleBuf; t0=0.0, downsample=:auto, pooling=:auto)
  ticks --> :native
  legend --> ndims(s) > 1 && size(s, 2) > 1
  s1 = s.data
  t = domain(s) .+ t0
  if maximum(t) <= 1.0
    t *= 1000.0
    xguide --> "Time (ms)"
  else
    xguide --> "Time (seconds)"
  end
  n = size(s1, 1)
  if downsample == :auto
    downsample = n >= 10000 ? ceil(Int, n/5000) : 1
    #downsample > 1 && @warn "Downsampling by $downsample"
  end
  if downsample !== nothing && downsample != 1
    if pooling === nothing || (pooling == :auto && isanalytic(s1))
      phase = floor(Int, downsample/2)
      s1 = s1[1+phase:downsample:end,:]
      t = t[1+phase:downsample:end]
    else
      if pooling == :minmax || pooling == :auto
        pooling = orderedextrema
      elseif pooling == :min
        pooling = minimum
      elseif pooling == :max
        pooling = maximum
      elseif pooling == :mean
        pooling = mean
      end
      if ndims(s) == 1
        s1 = collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s1, downsample, 0))))
        downsample = round(Int, n/size(s1, 1))
        t = t[1:downsample:end]
      else
        y = collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s1[:,1], downsample, 0))))
        if size(s1, 2) > 1
          y1 = zeros(length(y), size(s1, 2))
          y1[:,1] .= y
          y = y1
        end
        for j in 2:size(s1, 2)
          y[:,j] .= collect(Iterators.flatten(pooling.(Periodograms.arraysplit(s1[:,j], downsample, 0))))
        end
        s1 = y
        downsample = round(Int, n/size(s1, 1))
        t = t[1:downsample:end]
      end
    end
  end
  @series begin
    seriestype := :line
    t[1:length(s1)], s1
  end
end

"""
    psd(data; kwargs...)

Plots the power spectral density of data.

# Optional keyword arguments
- `fs=1.0`: derived from the data if a `SampleBuf` is provided as input
- `nfft=512`: size of FFT window
- `noverlap=nfft÷2`: window overlap size
- `window=nothing`: accepts any window from [`DSP.jl`](https://juliadsp.github.io/DSP.jl/stable/windows/)
- `xscale=:auto`: one of `:auto`, `:identity` or `:log10`
- other `kwargs` are passed on to `plot`
"""
@userplot PSD
@recipe function plot(s::PSD; fs=1.0, nfft=512, noverlap=div(nfft,2),
                      window=nothing, xscale=:auto)
  if length(s.args) != 1 || !(s.args[1] isa AbstractArray)
    error("psd should be provided timeseries data")
  end
  s.args[1] isa SampleBuf && (fs = framerate(s.args[1]))
  s = samples(s.args[1])
  nfft = nextfastfft(nfft)
  while nfft > size(s, 1)
    nfft = nfft ÷ 2
  end
  ticks --> :native
  legend --> ndims(s) > 1 && size(s, 2) > 1
  yguide --> "Power spectral density (dB/Hz)"
  xscale == :auto && (xscale := :identity)
  for j = 1:size(s, 2)
    p = welch_pgram(s[:,j], nfft, noverlap; fs=inHz(fs), window=window)
    f = freq(p)
    if maximum(f) > 10000 && xscale == :auto
      xguide --> "Frequency (kHz)"
      f /= 1000.0
    else
      xguide --> "Frequency (Hz)"
    end
    if xscale == :log10
      xlims --> (f[2], f[end])
    end
    y = pow2db.(power(p))
    ymax = ceil(Int, maximum(y)/5)*5
    ymin = max(ymax-50, floor(Int, minimum(y)/5)*5)
    ylims --> (ymin, ymax)
    @series begin
      seriestype := :line
      f, y
    end
  end
end

"""
    specgram(data; kwargs...)

Plots a spectrogram of the data.

# Optional keyword arguments
- `fs=1.0`: derived from the data if a `SampleBuf` is provided as input
- `nfft=256`: size of FFT window
- `noverlap=nfft÷2`: window overlap size
- `window=nothing`: accepts any window from [`DSP.jl`](https://juliadsp.github.io/DSP.jl/stable/windows/)
- `t0=0.0`: start time
- `downsample=:auto`: downsampling factor (integer) for time axis
- `pooling=:mean`: pooling mode (`:min`, `:max`, `:mean`, `nothing` or function)
- other `kwargs` are passed on to `plot`
"""
@userplot Specgram
@recipe function plot(s::Specgram; fs=1.0, nfft=256, noverlap=div(nfft,2),
                      window=nothing, t0=0.0, downsample=:auto, pooling=:mean)
  if length(s.args) != 1 || !(s.args[1] isa AbstractArray)
    error("specgram should be provided timeseries data")
  end
  s.args[1] isa SampleBuf && (fs = framerate(s.args[1]))
  s1 = samples(s.args[1])
  if ndims(s1) > 1 && size(s1, 2) > 1
    error("specgram does not support multichannel data")
  end
  nfft = nextfastfft(nfft)
  while nfft > size(s1, 1)
    nfft = nfft ÷ 2
  end
  ticks --> :native
  p = spectrogram(s1[:,1], nfft, noverlap; fs=inHz(fs), window=window)
  t = t0 .+ time(p)
  f = freq(p)
  z = power(p)
  if downsample == :auto
    downsample = length(t) > 10000 ? ceil(Int, length(t)/5000) : 1
    #downsample > 1 && @warn "Downsampling by $downsample"
  end
  if downsample !== nothing && downsample != 1
    if pooling === nothing
      phase = floor(Int, downsample/2)
      t = t[1+phase:downsample:end]
      z = z[:,1+phase:downsample:end]
    else
      if pooling == :min
        pooling = minimum
      elseif pooling == :max
        pooling = maximum
      elseif pooling == :mean
        pooling = mean
      end
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
  if isanalytic(s1)
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
  xguide --> "Time ("*tunit*")"
  yguide --> "Frequency ("*funit*")"
  z = pow2db.(z)
  cmax = ceil(Int, maximum(z)/5)*5
  cmin = max(cmax-50, floor(Int, minimum(z)/5)*5)
  clims --> (cmin, cmax)
  @series begin
    seriestype := :heatmap
    t, f, z
  end
end

"""
    freqresp(filter; kwargs...)
    freqresp(num; kwargs...)
    freqresp(num, den; kwargs...)

Plots frequency response of a digital filter.

# Optional keyword arguments
- `fs=1.0`: sampling frequency
- `nfreq=256`: number of frequency points to evaluate filter response at
- `xscale=:auto`: one of `:auto`, `:identity` or `:log10`
- other `kwargs` are passed on to `plot`
"""
@userplot FreqResp
@recipe function plot(s::FreqResp; fs=1.0, nfreq=256, xscale=:auto)
  if length(s.args) ∉ (1, 2)
    error("freqresp should be provided filter details")
  end
  if s.args[1] isa FilterCoefficients
    filt = s.args[1]
  else
    num = s.args[1] isa AbstractArray ? s.args[1] : [s.args[1]]
    if length(s.args) == 2
      den = s.args[2] isa AbstractArray ? s.args[2] : [s.args[2]]
      filt = PolynomialRatio(num, den)
    else
      filt = PolynomialRatio(num, [1.0])
    end
  end
  fs = inHz(fs)
  f = collect(range(0, stop=fs/2, length=nfreq))
  z = Filters.freqz(filt, f, fs)
  funit = "Hz"
  if xscale == :auto
    xscale := :identity
    if maximum(f) >= 10000
      f /= 1000.0
      funit = "kHz"
    end
  end
  if xscale == :log10
    f = f[2:end]
    z = z[2:end]
  end
  layout := @layout [mag; phase]
  @series begin
    seriestype := :line
    subplot := 1
    xguide --> "Frequency ("*funit*")"
    yguide --> "Response (dB)"
    legend --> false
    f, amp2db.(abs.(z))
  end
  @series begin
    seriestype := :line
    subplot := 2
    xguide --> "Frequency ("*funit*")"
    yguide --> "Phase (radians)"
    legend --> false
    f, angle.(z)
  end
end

# helper functions

orderedextrema(x) = argmin(x) < argmax(x) ? (minimum(x), maximum(x)) : (maximum(x), minimum(x))
