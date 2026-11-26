module PlotsExt

using SignalAnalysis
using SignalAnalysis: SampledSignal, SampledSpectrum
using Plots
using DSP.Periodograms

"""
    plot(data::SampledSignal; kwargs...)

Plots timeseries from a sample buffer.

# Optional keyword arguments
- `t0=0.0`: start time
- `downsample=:auto`: downsampling factor (integer)
- `pooling=:auto`: pooling mode (`:min`, `:max`, `:mean`, `:minmax`, `nothing` or function)
- `seriestype=:line`: series type, e.g. `:line`, `:stem` or `:scatter`
- other `kwargs` are passed on to `plot`

Long signals are downsampled for display. By default, line plots use `:minmax`
pooling and other series types use no pooling. Stem plots have circle markers
by default, and `marker=:none` removes them.
"""
@recipe function plot(s::SampledSignal; t0=0.0, downsample=:auto, pooling=:auto)
  ticks --> :native
  seriestype = get(plotattributes, :seriestype, :line)
  pooling == :auto && seriestype ∉ (:line, :path) && (pooling = nothing)
  legend --> ndims(s) > 1 && size(s, 2) > 1
  s1 = samples(s)
  eltype(s1) <: Complex && (s1 = abs.(s1))
  t = domain(s) .+ t0
  if maximum(t) <= 1.0
    t *= 1000.0
    xguide --> "Time (ms)"
  else
    xguide --> "Time (s)"
  end
  n = size(s1, 1)
  if downsample == :auto
    downsample = n >= 10000 ? ceil(Int, n/5000) : 1
    #downsample > 1 && @warn "Downsampling by $downsample"
  end
  if downsample !== nothing && downsample != 1
    if pooling === nothing || (pooling == :auto && isanalytic(s))
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
        m = fld(n, downsample)
        rep, rem = fldmod(size(s1, 1), m)
        rem == 0 || throw(ArgumentError("Pooling function must sample the signal evenly"))
        t = t[1:downsample:end]
        rep > 1 && (t = repeat(t, inner=rep))
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
        m = fld(n, downsample)
        rep, rem = fldmod(size(s1, 1), m)
        rem == 0 || throw(ArgumentError("Pooling function must sample the signal evenly"))
        t = t[1:downsample:end]
        rep > 1 && (t = repeat(t, inner=rep))
      end
    end
  end
  @series begin
    seriestype --> :line
    seriestype == :sticks && (markershape --> :circle)
    t[1:size(s1,1)], s1
  end
end

"""
    plot(X::SampledSpectrum; kwargs...)

Plots the magnitude (in dB) of a spectrum against frequency.
"""
@recipe function plot(X::SampledSpectrum)
  ticks --> :native
  legend --> ndims(X) > 1 && size(X, 2) > 1
  f = domain(X)
  if maximum(abs, f) >= 5000.0
    f = f ./ 1000.0
    xguide --> "Frequency (kHz)"
  else
    xguide --> "Frequency (Hz)"
  end
  yguide --> "Magnitude (dB)"
  @series begin
    seriestype --> :line
    f, amp2db.(abs.(samples(X)))
  end
end

"""
    plot(data::TFR; kwargs...)

Displays time frequency representation.

# Optional keyword arguments
- `t0=0.0`: start time
- `cscale=:log`: color scale type (`:log`, `:linear`)
- `crange=:auto`: color scale dynamic range (`:auto`, number)
- other `kwargs` are passed on to `plot`
"""
@recipe function plot(tfr::Periodograms.TFR; t0=0.0, cscale=:log, crange=:auto)
  ticks --> :native
  t = time(tfr) .+ t0
  f = freq(tfr)
  z = power(tfr)
  if maximum(t) <= 1.0
    t *= 1000.0
    xguide --> "Time (ms)"
  else
    xguide --> "Time (s)"
  end
  if maximum(f) >= 5000.0
    f /= 1000.0
    yguide --> "Frequency (kHz)"
  else
    yguide --> "Frequency (Hz)"
  end
  z = abs.(z)
  if cscale == :log
    z = pow2db.(z)
    cmax = ceil(maximum(z)/5)*5
    crange == :auto && (crange = 50)
    cmin = max(cmax-crange, floor(minimum(z)/5)*5)
    clims --> (cmin, cmax)
  elseif crange != :auto
    cmax = maximum(z)
    cmin = max(cmax-crange, minimum(z))
    clims --> (cmin, cmax)
  end
  @series begin
    seriestype := :heatmap
    t, f, z
  end
end


"""
    psd(data; kwargs...)

Plots the power spectral density of data.

# Optional keyword arguments
- `fs=1.0`: derived from the data if a `SampledSignal` is provided as input
- `nfft=512`: size of FFT window
- `noverlap=nfft÷2`: window overlap size
- `window=hamming(nfft)`: accepts any window from [`DSP.jl`](https://juliadsp.github.io/DSP.jl/stable/windows/)
- `xscale=:auto`: one of `:auto`, `:identity` or `:log10`
- `yrange=50`: y-scale dB range for automatic scaling
- other `kwargs` are passed on to `plot`
"""
@userplot PSD
@recipe function plot(s::PSD; fs=1.0, nfft=512, noverlap=div(nfft,2),
                      window=hamming(nfft), xscale=:auto, yrange=50)
  if length(s.args) != 1 || !(s.args[1] isa AbstractArray)
    throw(ArgumentError("psd should be provided timeseries data"))
  end
  s.args[1] isa SampledSignal && (fs = framerate(s.args[1]))
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
    ymin = max(ymax-yrange, floor(Int, minimum(y)/5)*5)
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
- `fs=1.0`: derived from the data if a `SampledSignal` is provided as input
- `nfft=256`: size of FFT window
- `noverlap=nfft÷2`: window overlap size
- `window=hamming(nfft)`: accepts any window from [`DSP.jl`](https://juliadsp.github.io/DSP.jl/stable/windows/)
- `t0=0.0`: start time
- `downsample=:auto`: downsampling factor (integer) for time axis
- `pooling=:mean`: pooling mode (`:min`, `:max`, `:mean`, `nothing` or function)
- `crange=50`: color scale dB range for automatic scaling
- other `kwargs` are passed on to `plot`
"""
@userplot Specgram
@recipe function plot(s::Specgram; fs=1.0, nfft=256, noverlap=div(nfft,2), crange=50,
                      window=hamming(nfft), t0=0.0, downsample=:auto, pooling=:mean)
  if length(s.args) != 1 || !(s.args[1] isa AbstractArray)
    throw(ArgumentError("specgram should be provided timeseries data"))
  end
  s.args[1] isa SampledSignal && (fs = framerate(s.args[1]))
  s1 = samples(s.args[1])
  if ndims(s1) > 1 && size(s1, 2) > 1
    throw(ArgumentError("specgram does not support multichannel data"))
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
  cmax = ceil(maximum(z)/5)*5
  cmin = max(cmax-crange, floor(minimum(z)/5)*5)
  clims --> (cmin, cmax)
  @series begin
    seriestype := :heatmap
    t, f, z
  end
end

"""
    freqz(filter; kwargs...)
    freqz(num; kwargs...)
    freqz(num, den; kwargs...)

Plots frequency response of a digital filter. The filter may be a DSP.jl filter
(such as one designed by [`iir`](@ref)), a vector of FIR taps (such as one
designed by [`fir`](@ref)), or numerator and denominator coefficients.

# Optional keyword arguments
- `fs=1.0`: sampling frequency
- `nfreq=256`: number of frequency points to evaluate filter response at
- `xscale=:auto`: one of `:auto`, `:identity` or `:log10`
- `show=:both`: one of `:both`, `:mag` or `:phase`
- `dbrange=30`: dynamic range (dB) of magnitude plot, or `nothing` for automatic
- `unwrap=true`: unwrap the phase response
- other `kwargs` are passed on to `plot`

Plot attributes apply to both magnitude and phase subplots. To set attributes
for each subplot separately, use a row matrix, e.g. `ylims=[(-60, 0) (-π, π)]`.
"""
@userplot FreqZ
@recipe function plot(s::FreqZ; fs=1.0, nfreq=256, xscale=:auto, show=:both, dbrange=30, unwrap=true)
  if length(s.args) ∉ (1, 2)
    throw(ArgumentError("freqz should be provided filter details"))
  end
  show ∈ (:both, :mag, :phase) || throw(ArgumentError("show should be one of :both, :mag or :phase"))
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
  z = Filters.freqresp(filt, 2π * f/fs)
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
  show == :both && (layout := @layout [mag; phase])
  if show ∈ (:both, :mag)
    mag = amp2db.(abs.(z))
    @series begin
      seriestype := :line
      subplot := 1
      xguide --> "Frequency ("*funit*")"
      yguide --> "Response (dB)"
      legend --> false
      if dbrange !== nothing
        maxdb = maximum(filter(isfinite, mag))
        ylims --> (maxdb - dbrange, maxdb + 3)
      end
      f, mag
    end
  end
  if show ∈ (:both, :phase)
    ϕ = angle.(z)
    @series begin
      seriestype := :line
      subplot := show == :both ? 2 : 1
      xguide --> "Frequency ("*funit*")"
      yguide --> "Phase (radians)"
      legend --> false
      f, unwrap ? DSP.unwrap(ϕ) : ϕ
    end
  end
end

# forward the userplot entry points to the SignalAnalysis stubs, since the
# @userplot-generated functions are local to this extension module
SignalAnalysis.psd(args...; kwargs...) = psd(args...; kwargs...)
SignalAnalysis.psd!(args...; kwargs...) = psd!(args...; kwargs...)
SignalAnalysis.specgram(args...; kwargs...) = specgram(args...; kwargs...)
SignalAnalysis.specgram!(args...; kwargs...) = specgram!(args...; kwargs...)
SignalAnalysis.freqz(args...; kwargs...) = freqz(args...; kwargs...)
SignalAnalysis.freqz!(args...; kwargs...) = freqz!(args...; kwargs...)

# helper functions

orderedextrema(x) = argmin(x) < argmax(x) ? (minimum(x), maximum(x)) : (maximum(x), minimum(x))

end # module
