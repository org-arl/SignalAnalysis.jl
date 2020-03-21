using Interact

export iplot, ispecgram

"""
$(SIGNATURES)
Plots interactive timeseries of the signal.
"""
function iplot(s; fs=framerate(s), t0=0.0, legend=false, size=(800,400), ylabel=nothing, kwargs...)
  fs = inHz(fs)
  s = samples(s)
  if isanalytic(s)
    @warn "Plotting only real part of complex signal"
    s = real.(s)
  end
  siglen = Base.size(s, 1)/fs
  panback = button("<")
  panfwd = button(">")
  pan = widget((0:(length(s)-1))/fs; value=0.0, label="Pan:") |> onchange
  zoomout = button("-")
  zoomin = button("+")
  zoom = widget((100:length(s))/fs; value=1000/fs, label="Zoom:") |> onchange
  #pooling = widget(OrderedDict("minmax"=>orderedextrema, "sample"=>nothing, "min"=>minimum, "max"=>maximum, "mean"=>mean); label="Pooling:")
  ymin = minimum(s)
  ymax = maximum(s)
  f = (t, z) -> begin
    t = round(Int, t*fs)
    z = round(Int, z*fs)
    t1 = max(t - z÷2, 1)
    t2 = min(t1 + z, Base.size(s,1))
    t1 = max(t2 - z, 1)
    ds = ceil(Int, (t2-t1)/(10*size[1]))
    #s1 = @view s[t1:t2,:]
    s1 = s[t1:t2,:]
    plot(signal(s1, fs); size=size, t0=t0+float(t1-1)/fs, downsample=ds, leg=legend, kwargs...)
    #timeseries(s[t1:t2,:]; size=size, fs=fs, t0=t0+float(t1-1)/fs, downsample=ds, pooling=pooling[], legend=legend, kwargs...)
    #ds == 1 || annotate!(t0+(t1+t2)/2/fs, ymin+(ymax-ymin)/20, text("Downsampled by $ds", :red, 8))
    ylabel == nothing || ylabel!(ylabel)
    ylims!(ymin, ymax)
  end
  on(n -> pan[] = max(pan[]-zoom[]/2, 0.0), panback)
  on(n -> pan[] = min(pan[]+zoom[]/2, siglen), panfwd)
  on(n -> zoom[] = min(zoom[]*1.5, siglen), zoomout)
  on(n -> zoom[] = max(zoom[]/1.5, 100/fs), zoomin)
  #on(n -> zoom[] = zoom[], pooling)
  vbox(
    hbox(panback, hskip(5px), panfwd, pan, zoomout, hskip(5px), zoomin, zoom),
    #hbox(pooling),
    map(f, pan, zoom)
  )
end

"""
$(SIGNATURES)
Plots interactive spectrogram of the signal.
"""
function ispecgram(s; fs=framerate(s), nfft=min(div(length(s),8),256), noverlap=div(nfft,2), window=nothing, t0=0.0, size=(800,400), kwargs...)
  @assert Base.size(s,2) == 1 "ispecgram only works with vectors"
  fs = inHz(fs)
  siglen = length(s)/fs
  panback = button("<")
  panfwd = button(">")
  pan = widget(0:(length(s)-1)/fs; value=0.0, label="Pan:") |> onchange
  zoomout = button("-")
  zoomin = button("+")
  zoom = widget(0.1:0.1:(length(s)-1)/fs; value=1.0, label="Zoom:") |> onchange
  pooling = widget(OrderedDict("mean"=>mean, "median"=>median, "min"=>minimum, "max"=>maximum, "sample"=>nothing); label="Pooling:")
  f = (t, z) -> begin
    t = round(Int, t*fs)
    z = max(round(Int, z*fs), nfft)
    t1 = max(t - z÷2, 1)
    t2 = min(t1 + z, Base.size(s,1))
    t1 = max(t2 - z, 1)
    ds = ceil(Int, (t2-t1)/(10*size[1]))
    s1 = @view s[t1:t2]
    specgram(s1; size=size, fs=fs, t0=t0+float(t1-1)/fs, downsample=ds, pooling=pooling[], nfft=nfft, noverlap=noverlap, window=window, kwargs...)
    #ds == 1 || annotate!(t0+(t1+t2)/2/fs, ymin+(ymax-ymin)/20, text("Downsampled by $ds", :red, 8))
    #ylims!(ymin, ymax)
  end
  on(n -> pan[] = max(pan[]-zoom[]/2, 0.0), panback)
  on(n -> pan[] = min(pan[]+zoom[]/2, siglen), panfwd)
  on(n -> zoom[] = min(zoom[]*1.5, siglen), zoomout)
  on(n -> zoom[] = max(zoom[]/1.5, 100/fs), zoomin)
  on(n -> zoom[] = zoom[], pooling)
  vbox(
    hbox(panback, hskip(5px), panfwd, pan, zoomout, hskip(5px), zoomin, zoom, pooling),
    map(f, pan, zoom)
  )
end
