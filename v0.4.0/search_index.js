var documenterSearchIndex = {"docs":
[{"location":"generate.html#Generating-signals","page":"Generating signals","title":"Generating signals","text":"","category":"section"},{"location":"generate.html","page":"Generating signals","title":"Generating signals","text":"Modules = [SignalAnalysis]\nPages   = [\"generate.jl\"]","category":"page"},{"location":"generate.html#SignalAnalysis.chirp","page":"Generating signals","title":"SignalAnalysis.chirp","text":"chirp(freq1, freq2, duration)\nchirp(freq1, freq2, duration, fs; shape, phase, window)\n\n\nGenerates a frequency modulated chirp signal from freq1 to freq2 and specified duration at frame rate fs. The type of frequency modulation may be controlled using shape (:linear (default) or :hyperbolic). The starting phase and window type may be optionally specified.\n\nExamples:\n\njulia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz)\nSampledSignal @ 44100.0 Hz, 4410-element Array{Complex{Float64},1}:\n  ⋮\n\njulia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz; phase=45°, window=hamming)\nSampledSignal @ 44100.0 Hz, 4410-element Array{Complex{Float64},1}:\n  ⋮\n\njulia> x = chirp(5kHz, 7kHz, 100ms, 44.1kHz; shape=:hyperbolic, window=(tukey,0.05))\nSampledSignal @ 44100.0 Hz, 4410-element Array{Complex{Float64},1}:\n  ⋮\n\n\n\n\n\n","category":"function"},{"location":"generate.html#SignalAnalysis.cw-Tuple{Any,Any,Any}","page":"Generating signals","title":"SignalAnalysis.cw","text":"cw(freq, duration, fs; phase, window)\n\n\nGenerates a sinusoidal signal with specified freq and duration at frame rate fs. The starting phase and window type may be optionally specified.\n\nExamples:\n\njulia> x = cw(5kHz, 200ms, 44.1kHz)\nSampledSignal @ 44100.0 Hz, 8820-element Array{Complex{Float64},1}:\n  ⋮\n\njulia> x = cw(5kHz, 200ms, 44.1kHz; window=hamming)\nSampledSignal @ 44100.0 Hz, 8820-element Array{Complex{Float64},1}:\n  ⋮\n\njulia> x = cw(-5kHz, 200ms, 44.1kHz; phase=45°, window=(tukey, 0.05))\nSampledSignal @ 44100.0 Hz, 8820-element Array{Complex{Float64},1}:\n  ⋮\n\n\n\n\n\n","category":"method"},{"location":"iplot.html#Interactive-plotting","page":"Interactive plotting","title":"Interactive plotting","text":"","category":"section"},{"location":"iplot.html","page":"Interactive plotting","title":"Interactive plotting","text":"Interactive plots use InteractiveViz.jl to allow viewing of long signals.","category":"page"},{"location":"iplot.html","page":"Interactive plotting","title":"Interactive plotting","text":"Modules = [SignalAnalysis]\nPages   = [\"iplot.jl\"]","category":"page"},{"location":"iplot.html#InteractiveViz.iplot!-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Interactive plotting","title":"InteractiveViz.iplot!","text":"iplot!(s; kwargs...)\n\n\nPlots interactive timeseries of the signal over a previous plot.\n\n\n\n\n\n","category":"method"},{"location":"iplot.html#InteractiveViz.iplot-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T,Vararg{Any,N} where N}","page":"Interactive plotting","title":"InteractiveViz.iplot","text":"iplot(s, args; kwargs...)\n\n\nPlots interactive timeseries of the signal.\n\n\n\n\n\n","category":"method"},{"location":"iplot.html#SignalAnalysis.ispecgram-Tuple{Any}","page":"Interactive plotting","title":"SignalAnalysis.ispecgram","text":"ispecgram(s; fs, nfft, noverlap, window, pooling, kwargs...)\n\n\nPlots interactive spectrogram of the signal.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#Basic-signal-analysis","page":"Basic signal analysis","title":"Basic signal analysis","text":"","category":"section"},{"location":"basic.html","page":"Basic signal analysis","title":"Basic signal analysis","text":"Modules = [SignalAnalysis]\nPages   = [\"basic.jl\"]","category":"page"},{"location":"basic.html#SignalAnalysis.energy-Tuple{AbstractArray{T,1} where T}","page":"Basic signal analysis","title":"SignalAnalysis.energy","text":"energy(s; fs)\n\n\nComputes total signal energy.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#SignalAnalysis.ifrequency-Tuple{Any}","page":"Basic signal analysis","title":"SignalAnalysis.ifrequency","text":"ifrequency(s; fs)\n\n\nComputes instantaneous frequency of the signal.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#SignalAnalysis.meanfrequency-Tuple{AbstractArray{T,2} where T}","page":"Basic signal analysis","title":"SignalAnalysis.meanfrequency","text":"meanfrequency(s; fs, nfft, window)\n\n\nComputes mean frequency of a signal.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#SignalAnalysis.meantime-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Basic signal analysis","title":"SignalAnalysis.meantime","text":"meantime(s)\n\n\nComputes mean time of the signal.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#SignalAnalysis.rmsbandwidth-Union{Tuple{AbstractArray{T,2} where T}, Tuple{T}} where T","page":"Basic signal analysis","title":"SignalAnalysis.rmsbandwidth","text":"rmsbandwidth(s; fs, nfft, window)\n\n\nComputes RMS bandwidth of a signal.\n\n\n\n\n\n","category":"method"},{"location":"basic.html#SignalAnalysis.rmsduration-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Basic signal analysis","title":"SignalAnalysis.rmsduration","text":"rmsduration(s)\n\n\nComputes RMS duration of the signal.\n\n\n\n\n\n","category":"method"},{"location":"array.html#Array-signal-processing","page":"Array signal processing","title":"Array signal processing","text":"","category":"section"},{"location":"array.html","page":"Array signal processing","title":"Array signal processing","text":"Modules = [SignalAnalysis]\nPages   = [\"array.jl\"]","category":"page"},{"location":"array.html#SignalAnalysis.Bartlett","page":"Array signal processing","title":"SignalAnalysis.Bartlett","text":"Frequency-domain Bartlett beamformer.\n\n\n\n\n\n","category":"type"},{"location":"array.html#SignalAnalysis.Capon","page":"Array signal processing","title":"SignalAnalysis.Capon","text":"Frequency-domain Capon beamformer with diagonal loading factor γ.\n\n\n\n\n\n","category":"type"},{"location":"array.html#SignalAnalysis.Music","page":"Array signal processing","title":"SignalAnalysis.Music","text":"Frequency-domain MUSIC beamformer with nsignals signals.\n\n\n\n\n\n","category":"type"},{"location":"array.html#SignalAnalysis.beamform-NTuple{4,Any}","page":"Array signal processing","title":"SignalAnalysis.beamform","text":"beamform(s, f, n, sd; fs=framerate(s), method=Bartlett())\nbeamform(s, f, sd; fs=framerate(s), method=Bartlett())\n\nNarrowband frequency-domain beamformer. Takes in passband signals s and produces beamformer output for all directions specified by the steering delays sd. The beamformer output is an energy estimate (or equivalent) for each steering direction. The beamforming only uses a narrowband signal cenetered at frequency f with a bandwidth of about fs/n.\n\nIf n is not specified, or is 1, then the input signal is assumed to be narrowband, and centered at frequency f.\n\nThe narrowband assumption requires that the bandwidth be no greater than about 5/T, where T is the maximum time taken for a signal to propagate through the array.\n\nSeveral beamforming methods are available:\n\nBartlett()\nCapon(γ)\nMusic(nsignals)\n\nCustom beamformers can be implemented by creating a subtype of SignalAnalysis.Beamformer and implementing the SignalAnalysis.beamformer() method dispatched on that type.\n\nExample:\n\njulia> x = cw(100.0, 1.0, 44100.0);\njulia> sd = steering(0.0:1.0:3.0, 1500.0, range(0.0, π; length=181));\njulia> bfo = beamform([x x x x], 100.0, 4096, sd; method=Capon(0.1))\n181-element Array{Float64,1}:\n 0.12406290296318974\n 0.1240975045516605\n ⋮\n 0.12406290296318974\n\n\n\n\n\n","category":"method"},{"location":"array.html#SignalAnalysis.beamform-Tuple{Any,Any}","page":"Array signal processing","title":"SignalAnalysis.beamform","text":"beamform(s, sd; fs=framerate(s))\n\nBroadband time-domain delay-and-sum beamformer. Takes in passband or baseband signals s and produces beamformer output for all directions specified by the steering delays sd. The beamformer output is a timeseries signal for each steering direction.\n\nExample:\n\njulia> x = cw(100.0, 1.0, 44100.0);\njulia> sd = steering(0.0:1.0:3.0, 1500.0, range(0.0, π; length=181));\njulia> bfo = beamform([x x x x], sd)\nSampledSignal @ 44100.0 Hz, 44100×181 Array{Complex{Float64},2}:\n  ⋮\n\n\n\n\n\n","category":"method"},{"location":"array.html#SignalAnalysis.steering-Tuple{AbstractArray{T,2} where T,Any,AbstractArray{T,1} where T}","page":"Array signal processing","title":"SignalAnalysis.steering","text":"steering(rxpos, c, θ)\n\n\nComputes steering delays for specified receiver positions rxpos, signal propagation speed c, and angles θ. Array θ can be a 1D array of angles or 2D array with (azimuth, elevation) pair in each row. The delays are computed with a far-field assumption, i.e., for plane incoming waves.\n\nExamples:\n\njulia> steering(0.0:1.0:5.0, 1500.0, range(0.0, π; length=181))\n6×181 Array{Float64,2}:\n  0.00166667    0.00166641   …  -0.00166641   -0.00166667\n  0.001         0.000999848     -0.000999848  -0.001\n  0.000333333   0.000333283     -0.000333283  -0.000333333\n -0.000333333  -0.000333283      0.000333283   0.000333333\n -0.001        -0.000999848      0.000999848   0.001\n -0.00166667   -0.00166641   …   0.00166641    0.00166667\n\njulia> rxpos = [  # can be 2D or 3D coordinates\n  0.0  1.0  2.0  3.0  4.0  5.0\n  0.0  0.0  0.0  0.0  0.0  0.0\n];\njulia> steering(rxpos, 1500.0, range(0.0, π; length=181))\n6×181 Array{Float64,2}:\n  0.00166667    0.00166641   …  -0.00166641   -0.00166667\n  0.001         0.000999848     -0.000999848  -0.001\n  0.000333333   0.000333283     -0.000333283  -0.000333333\n -0.000333333  -0.000333283      0.000333283   0.000333333\n -0.001        -0.000999848      0.000999848   0.001\n -0.00166667   -0.00166641   …   0.00166641    0.00166667\n\njulia> rxpos = [\n  0.0  0.0  0.5  0.5\n  0.0  0.5  0.0  0.5\n  0.0  0.0  0.0  0.0\n];\njulia> θ = deg2rad.(reduce(vcat, hcat.(LinRange(-20,20,41)', LinRange(-10,10,21)))) # 2D array with (azimuth, elevation) in each row\njulia> steering(rxpos, 1500.0, θ)\n4×861 Array{Float64,2}:\n  9.80987e-5    9.83857e-5    9.86427e-5  …   0.00021154   0.000210989   0.000210373\n  0.000210373   0.000210989   0.00021154      9.86427e-5   9.83857e-5    9.80987e-5\n -0.000210373  -0.000210989  -0.00021154     -9.86427e-5  -9.83857e-5   -9.80987e-5\n -9.80987e-5   -9.83857e-5   -9.86427e-5     -0.00021154  -0.000210989  -0.000210373\n\n\n\n\n\n","category":"method"},{"location":"plot.html#Plot-recipes","page":"Plot recipes","title":"Plot recipes","text":"","category":"section"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"Plot recipes are enabled by importing Plots:","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"using Plots","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"We provide a plot recipe to display signals:","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"plot(data::SampledSignal; t0=0.0, downsample=:auto, pooling=:auto, kwargs...)","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"Optional arguments:","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"t0=0.0: start time (for labeling only)\ndownsample=:auto: downsampling factor (integer or :auto)\npooling=:auto: pooling mode (:auto, :min, :max, :mean, :minmax, nothing or function)","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"If the signal is too long, it is automatically downsampled in a perceptually meaningful way. The downsampling can be controlled using the downsample and pooling keywords.","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"We also provide convenience plot recipes for common signal processing plots:","category":"page"},{"location":"plot.html","page":"Plot recipes","title":"Plot recipes","text":"psd\nspecgram\nfreqresp","category":"page"},{"location":"plot.html#SignalAnalysis.psd","page":"Plot recipes","title":"SignalAnalysis.psd","text":"psd(data; kwargs...)\n\nPlots the power spectral density of data.\n\nOptional keyword arguments\n\nfs=1.0: derived from the data if a SampledSignal is provided as input\nnfft=512: size of FFT window\nnoverlap=nfft÷2: window overlap size\nwindow=hamming(nfft): accepts any window from DSP.jl\nxscale=:auto: one of :auto, :identity or :log10\nyrange=50: y-scale dB range for automatic scaling\nother kwargs are passed on to plot\n\n\n\n\n\n","category":"function"},{"location":"plot.html#SignalAnalysis.specgram","page":"Plot recipes","title":"SignalAnalysis.specgram","text":"specgram(data; kwargs...)\n\nPlots a spectrogram of the data.\n\nOptional keyword arguments\n\nfs=1.0: derived from the data if a SampledSignal is provided as input\nnfft=256: size of FFT window\nnoverlap=nfft÷2: window overlap size\nwindow=hamming(nfft): accepts any window from DSP.jl\nt0=0.0: start time\ndownsample=:auto: downsampling factor (integer) for time axis\npooling=:mean: pooling mode (:min, :max, :mean, nothing or function)\ncrange=50: color scale dB range for automatic scaling\nother kwargs are passed on to plot\n\n\n\n\n\n","category":"function"},{"location":"plot.html#SignalAnalysis.freqresp","page":"Plot recipes","title":"SignalAnalysis.freqresp","text":"freqresp(filter; kwargs...)\nfreqresp(num; kwargs...)\nfreqresp(num, den; kwargs...)\n\nPlots frequency response of a digital filter.\n\nOptional keyword arguments\n\nfs=1.0: sampling frequency\nnfreq=256: number of frequency points to evaluate filter response at\nxscale=:auto: one of :auto, :identity or :log10\nother kwargs are passed on to plot\n\n\n\n\n\n","category":"function"},{"location":"signals.html#Creating-and-managing-signals","page":"Creating & managing signals","title":"Creating & managing signals","text":"","category":"section"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"Signals are represented by array-like data types with a series (time, frequency or space domain) per column. The SampledSignal data type wraps an array and carries sampling rate as metadata with it. This allows the API to be used without having to specify sampling rate at each call. However, if a user prefers to use other array-like data types, the sampling rate may be provided as a fs keyword argument for API calls that require sampling rate.","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"All code examples in this manual assume you have imported SignalAnalysis and SignalAnalysis.Units:","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"using SignalAnalysis, SignalAnalysis.Units","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"The SignalAnalysis.Units package re-exports commonly used units (s, ms, Hz, kHz, etc) from Unitful.jl. In addition, a variable 𝓈 is always exported and is an alias for the s unit from SignalAnalysis.Units. This allows indexing of signals using time in seconds:","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"# create a 1-second signal sampled at 20 kHz\nx = signal(randn(20000), 20kHz)\n\n# get a signal segment from 0.25 to 0.5 seconds:\ny = x[0.25:0.5𝓈]","category":"page"},{"location":"signals.html#Creating-and-wrapping-signals","page":"Creating & managing signals","title":"Creating and wrapping signals","text":"","category":"section"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"Signals wrapped in a SampledSignal data type may be easily created using signal():","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"x = signal(data, fs)","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"Properties such as frame rate (sampling rate), number of channels, etc may be accessed using the SignalBase API. The signals can be treated as arrays, but carry sampling rate metadata with them. While most operations infer metadata from the input signal, some operations may be unable to automatically infer the frame rate of the output signal. We provide some handy wrappers around common DSP.jl functions to aid with rate inference:","category":"page"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"# create a 1-second signal sampled at 20 kHz\ny = signal(randn(20000), 20kHz)\n\n# use DSP.filt() to filter a signal but retainin sampling rate\ny = sfilt(lpf, y)\n\n# use DSP.filtfilt() to filter a signal but retainin sampling rate\ny = sfiltfilt(lpf, y)\n\n# use DSP.resample() to resample a signal and infer sampling rate\ny = sresample(y, 2//3)","category":"page"},{"location":"signals.html#API-reference","page":"Creating & managing signals","title":"API reference","text":"","category":"section"},{"location":"signals.html","page":"Creating & managing signals","title":"Creating & managing signals","text":"Modules = [SignalAnalysis]\nPages   = [\"signals.jl\"]","category":"page"},{"location":"signals.html#Base.Colon-Tuple{Union{Unitful.Quantity{T,𝐓,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓,U}} where S where L} where U where T,Union{Unitful.Quantity{T,𝐓,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓,U}} where S where L} where U where T}","page":"Creating & managing signals","title":"Base.Colon","text":"(:)(start::Unitful.Time, stop::Unitful.Time)\n\nGenerates a time range index for a signal.\n\nExamples:\n\njulia> x = signal(randn(2000), 8kHz)\njulia> x[0.2𝓈:0.201𝓈]\nSampledSignal @ 8000.0 Hz, 9-element Array{Float64,1}:\n -0.08671384898800058\n -0.665143340284631\n -0.3955367460364236\n  1.2386430598616671\n -0.4882254309443194\n -1.080437097803303\n  0.8209785486953832\n  1.3477512734963886\n -0.27722340584395494\n\n\n\n\n\n","category":"method"},{"location":"signals.html#Base.Iterators.partition-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T,Integer}","page":"Creating & managing signals","title":"Base.Iterators.partition","text":"partition(x, n; step=n, flush=true)\n\nIterates over the signal x, n samples at a time, with a step size of step. If flush is enabled, the last partition may be smaller than n samples.\n\nWhen applied to a multichannel signal x, each partition contains samples from all channels.\n\nExamples:\n\n```julia-repl julia> x = signal(collect(1:10), 1.0); julia> collect(partition(x, 5)) 2-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:  [1, 2, 3, 4, 5]  [6, 7, 8, 9, 10]\n\njulia> collect(partition(x, 5; step=2)) 5-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:  [1, 2, 3, 4, 5]  [3, 4, 5, 6, 7]  [5, 6, 7, 8, 9]  [7, 8, 9, 10]  [9, 10]\n\njulia> collect(partition(x, 5; step=2, flush=false)) 3-element Array{SubArray{Int64,1,Array{Int64,1},Tuple{UnitRange{Int64}},true},1}:  [1, 2, 3, 4, 5]  [3, 4, 5, 6, 7]  [5, 6, 7, 8, 9]\n\njulia> x = signal(hcat(collect(1:10), collect(11:20)), 1.0); julia> collect(partition(x, 5))[1] 5×2 view(::Array{Int64,2}, 1:5, :) with eltype Int64:  1  11  2  12  3  13  4  14  5  15  ```\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.analytic-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Creating & managing signals","title":"SignalAnalysis.analytic","text":"analytic(s)\n\n\nConverts a signal to analytic representation.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.domain-Tuple{Any}","page":"Creating & managing signals","title":"SignalAnalysis.domain","text":"domain(x)\n\n\nReturns the domain of the signal.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.isanalytic-Tuple{Any}","page":"Creating & managing signals","title":"SignalAnalysis.isanalytic","text":"isanalytic(s)\n\n\nChecks if a signal is analytic.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.padded-Union{Tuple{T}, Tuple{AbstractArray{T,1},Any}} where T","page":"Creating & managing signals","title":"SignalAnalysis.padded","text":"padded(s, padding; delay, fill)\n\n\nGenerates a padded view of a signal with optional delay/advance.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.samples-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Creating & managing signals","title":"SignalAnalysis.samples","text":"samples(s)\n\n\nGets the underlying samples in the signal.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{AbstractArray,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(x::AbstractArray, fs::Any) -> MetaArrays.MetaArray{_A,_B,_C,_D} where _D where _C where _B where _A\n\n\nCreates a signal with frame rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{AbstractString}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(filename::AbstractString; start, nsamples) -> MetaArrays.MetaArray{_A,_B,_C,_D} where _D where _C where _B where _A\n\n\nLoads a signal from a WAV file.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(fs::Any) -> MetaArrays.MetaArray{_A,_B,_C,_D} where _D where _C where _B where _A\n\n\nCreates a curried function that takes in an array and creates a signal with sampling rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{Int64,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(n::Int64, fs::Any) -> MetaArrays.MetaArray{Array{Float64,1},SignalAnalysis.SamplingInfo,Float64,1}\n\n\nCreates an empty signal of length n samples, and frame rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{Int64,Int64,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(n::Int64, ch::Int64, fs::Any) -> MetaArrays.MetaArray{Array{Float64,2},SignalAnalysis.SamplingInfo,Float64,2}\n\n\nCreates an empty signal of length n samples, ch channels, and frame rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(x::MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T, fs::Any) -> MetaArrays.MetaArray{_A,_B,_C,_D} where _D where _C where _B where _A\n\n\nCreates a signal with frame rate fs. If the original signal's frame rate is the same as fs, this method simply returns the original signal. Otherwise, it creates a new signal with the specified frame rate and data from the original signal. Do note that this method does not resample the signal.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{Type,Int64,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(T::Type, n::Int64, fs::Any) -> MetaArrays.MetaArray{_A,SignalAnalysis.SamplingInfo,_B,1} where _B where _A\n\n\nCreates an empty signal of type T, length n samples, and frame rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.signal-Tuple{Type,Int64,Int64,Any}","page":"Creating & managing signals","title":"SignalAnalysis.signal","text":"signal(T::Type, n::Int64, ch::Int64, fs::Any) -> MetaArrays.MetaArray{_A,SignalAnalysis.SamplingInfo,_B,2} where _B where _A\n\n\nCreates an empty signal of type T, length n samples, ch channels, and frame rate fs.\n\n\n\n\n\n","category":"method"},{"location":"signals.html#SignalAnalysis.toframe-Tuple{Any,MetaArrays.MetaArray{#s12,SignalAnalysis.SamplingInfo,T,N} where N where #s12 where T}","page":"Creating & managing signals","title":"SignalAnalysis.toframe","text":"toframe(t, s)\n\n\nConverts time to signal frame number.\n\nExamples:\n\njulia> x = signal(randn(2000), 8kHz);\njulia> toframe(0.2s, x)\n1601\n\njulia> toframe([0.2s, 201ms], x)\n2-element Array{Int64,1}:\n 1601\n 1609\n\njulia> toframe(0.2:0.01:0.3, x)\n 11-element Array{Int64,1}:\n  1601\n  1681\n  1761\n   ⋮\n\n\n\n\n\n","category":"method"},{"location":"tfa.html#Time-frequency-analysis","page":"Time-frequency analysis","title":"Time-frequency analysis","text":"","category":"section"},{"location":"tfa.html","page":"Time-frequency analysis","title":"Time-frequency analysis","text":"Modules = [SignalAnalysis]\nPages   = [\"tfa.jl\"]","category":"page"},{"location":"tfa.html#SignalAnalysis.tfd-Tuple{Any,Spectrogram}","page":"Time-frequency analysis","title":"SignalAnalysis.tfd","text":"tfd(s::Any, kernel::Spectrogram; onesided, fs) -> SignalAnalysis.TFD{_A} where _A\n\n\nComputes a spectrogram time-frequency distribution of signal s sampled as sampling rate fs.\n\nExamples:\n\njulia> x = real.(chirp(1kHz, 10kHz, 1s, 44.1kHz));\njulia> y = tfd(x, Spectrogram());\njulia> plot(y)\njulia> y = tfd(x, Spectrogram(nfft=512, noverlap=256, window=hamming));\njulia> plot(y)\n\n\n\n\n\n","category":"method"},{"location":"tfa.html#SignalAnalysis.tfd-Tuple{Any,Wigner}","page":"Time-frequency analysis","title":"SignalAnalysis.tfd","text":"tfd(s::Any, kernel::Wigner; onesided, fs) -> SignalAnalysis.TFD{_A} where _A\n\n\nComputes a Wigner-Ville time-frequency distribution of signal s sampled as sampling rate fs.\n\nExamples:\n\njulia> x = real.(chirp(1kHz, 10kHz, 0.01s, 44.1kHz));\njulia> y = tfd(x, Wigner());\njulia> plot(y; clim=(0,20))\njulia> y = tfd(x, Wigner(nfft=512, smooth=10, method=:CM1980, window=hamming));\njulia> plot(y; clim=(0,20))\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#Signal-processing","page":"Signal processing","title":"Signal processing","text":"","category":"section"},{"location":"dsp.html","page":"Signal processing","title":"Signal processing","text":"Modules = [SignalAnalysis]\nPages   = [\"dsp.jl\"]","category":"page"},{"location":"dsp.html#SignalAnalysis.circconv","page":"Signal processing","title":"SignalAnalysis.circconv","text":"circconv(x)\ncircconv(x, y)\n\n\nComputes the circular convolution of x and y. Both vectors must be the same length.\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.demon-Tuple{Any}","page":"Signal processing","title":"SignalAnalysis.demon","text":"demon(x; fs, downsample, method, cutoff)\n\n\nEstimates DEMON spectrum. The output is highpass filtered with a cutoff frequency and downsampled. Supported downsampling methods are :rms (default), :mean and :fir.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.downconvert","page":"Signal processing","title":"SignalAnalysis.downconvert","text":"downconvert(s, sps, fc)\ndownconvert(s, sps, fc, pulseshape; fs)\n\n\nConverts passband signal centered around carrier frequency fc to baseband, and downsamples it by a factor of sps. If the pulseshape is specified to be nothing, downsampling is performed without filtering.\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.fir","page":"Signal processing","title":"SignalAnalysis.fir","text":"fir(n, f1)\nfir(n, f1, f2; fs, method)\n\n\nDesigns a n-tap FIR filter with a passband from f1 to f2 using the specified method. If frame rate fs is not specified, f1 and f2 are given in normalized units (1.0 being Nyquist). If f1 is 0, the designed filter is a lowpass filter, and if f2 is nothing then it is a highpass filter.\n\nThis method is a convenience wrapper around DSP.digitalfilter.\n\nExamples:\n\njulia> lpf = fir(127, 0, 10kHz; fs=44.1kHz)   # design a lowpass filter\n127-element Array{Float64,1}:\n  ⋮\n\njulia> hpf = fir(127, 10kHz; fs=44.1kHz)      # design a highpass filter\n127-element Array{Float64,1}:\n  ⋮\n\njulia> bpf = fir(127, 1kHz, 5kHz; fs=44.1kHz) # design a bandpass filter\n127-element Array{Float64,1}:\n  ⋮\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.gmseq","page":"Signal processing","title":"SignalAnalysis.gmseq","text":"gmseq(m)\ngmseq(m, θ)\n\n\nGenerates an generalized m-sequence of length 2^m-1 or tap specification m.\n\nGeneralized m-sequences are related to m-sequences but have an additional parameter θ. When θ = π/2, generalized m-sequences become normal m-sequences. When θ < π/2, generalized m-sequences contain a DC-component that leads to an exalted carrier after modulation. When θ is atan(√(2^m-1)), the m-sequence is considered to be period matched. Period matched m-sequences are complex sequences with perfect discrete periodic auto-correlation properties, i.e., all non-zero lag periodic auto-correlations are zero. The zero-lag autocorrelation is 2^m-1, where m is the shift register length.\n\nThis function currently supports shift register lengths between 2 and 30.\n\nExamples:\n\njulia> x = gmseq(3)         # generate period matched m-sequence\n7-element Array{Complex{Float64},1}:\n 0.3535533905932738 + 0.9354143466934853im\n 0.3535533905932738 + 0.9354143466934853im\n 0.3535533905932738 + 0.9354143466934853im\n 0.3535533905932738 - 0.9354143466934853im\n 0.3535533905932738 + 0.9354143466934853im\n 0.3535533905932738 - 0.9354143466934853im\n 0.3535533905932738 - 0.9354143466934853im\n\njulia> x = gmseq(3, π/4)    # generate m-sequence with exalted carrier\n7-element Array{Complex{Float64},1}:\n 0.7071067811865476 + 0.7071067811865475im\n 0.7071067811865476 + 0.7071067811865475im\n 0.7071067811865476 + 0.7071067811865475im\n 0.7071067811865476 - 0.7071067811865475im\n 0.7071067811865476 + 0.7071067811865475im\n 0.7071067811865476 - 0.7071067811865475im\n 0.7071067811865476 - 0.7071067811865475im\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.goertzel-Tuple{AbstractArray{T,1} where T,Any,Any}","page":"Signal processing","title":"SignalAnalysis.goertzel","text":"goertzel(x, f, n; fs)\n\n\nDetects frequency f in input signal using the Goertzel algorithm.\n\nThe detection metric returned by this function is the complex output of the Goertzel filter at the end of the input block. Typically, you would want to compare the magnitude of this output with a threshold to detect a frequency.\n\nWhen a block size n is specified, the Goertzel algorithm in applied to blocks of data from the original time series.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.mfilter-Tuple{Any,Any}","page":"Signal processing","title":"SignalAnalysis.mfilter","text":"mfilter(r, s)\n\n\nMatched filter looking for reference signal r in signal s.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.mseq","page":"Signal processing","title":"SignalAnalysis.mseq","text":"mseq(m)\nmseq(m, θ)\n\n\nGenerates an m-sequence of length 2^m-1 or tap specification m.\n\nm-sequences are sequences of +1/-1 values with near-perfect discrete periodic auto-correlation properties. All non-zero lag periodic auto-correlations are -1. The zero-lag autocorrelation is 2^m-1, where m is the shift register length.\n\nThis function currently supports shift register lengths between 2 and 30.\n\nExamples:\n\njulia> x = mseq(3)                  # generate regular m-sequence\n7-element Array{Float64,1}:\n  1.0\n  1.0\n  1.0\n -1.0\n  1.0\n -1.0\n -1.0\n\njulia> x = mseq((1,3))              # generate m-sequence with specification (1,3)\n7-element Array{Float64,1}:\n  1.0\n  1.0\n  1.0\n -1.0\n  1.0\n -1.0\n -1.0\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.pll","page":"Signal processing","title":"SignalAnalysis.pll","text":"pll(x)\npll(x, bandwidth; fs)\n\n\nPhased-lock loop to track dominant carrier frequency in input signal.\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.rcosfir","page":"Signal processing","title":"SignalAnalysis.rcosfir","text":"rcosfir(β, sps)\nrcosfir(β, sps, span)\n\n\nRaised cosine filter.\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.removedc!-Tuple{Any}","page":"Signal processing","title":"SignalAnalysis.removedc!","text":"removedc!(s; α)\n\n\nDC removal filter. Parameter α controls the cutoff frequency. Implementation based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).\n\nSee also: removedc\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.removedc-Tuple{Any}","page":"Signal processing","title":"SignalAnalysis.removedc","text":"removedc(s; α)\n\n\nDC removal filter. Parameter α controls the cutoff frequency. Implementation based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d).\n\nSee also: removedc!\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.rrcosfir","page":"Signal processing","title":"SignalAnalysis.rrcosfir","text":"rrcosfir(β, sps)\nrrcosfir(β, sps, span)\n\n\nRoot-raised cosine filter.\n\n\n\n\n\n","category":"function"},{"location":"dsp.html#SignalAnalysis.sfilt-Tuple{Any,Any,Vararg{Any,N} where N}","page":"Signal processing","title":"SignalAnalysis.sfilt","text":"sfilt(f, x, args)\n\n\nSame as filt, but retains sampling rate information.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.sfiltfilt-Tuple{Any,Any}","page":"Signal processing","title":"SignalAnalysis.sfiltfilt","text":"sfiltfilt(coef, x)\n\n\nSame as filtfilt, but retains sampling rate information.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.sresample-Tuple{Any,Any,Vararg{Any,N} where N}","page":"Signal processing","title":"SignalAnalysis.sresample","text":"sresample(x, rate, args)\n\n\nSame as resample, but correctly handles sampling rate conversion.\n\n\n\n\n\n","category":"method"},{"location":"dsp.html#SignalAnalysis.upconvert","page":"Signal processing","title":"SignalAnalysis.upconvert","text":"upconvert(s, sps, fc)\nupconvert(s, sps, fc, pulseshape; fs)\n\n\nConverts baseband signal with sps symbols per passband sample to a real passband signal centered around carrier frequency fc.\n\n\n\n\n\n","category":"function"},{"location":"index.html#SignalAnalysis.jl","page":"Home","title":"SignalAnalysis.jl","text":"","category":"section"},{"location":"index.html#Signal-analysis-toolbox-for-Julia","page":"Home","title":"Signal analysis toolbox for Julia","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule = SignalAnalysis","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"While a few great signal processing packages (e.g. DSP.jl, SignalOperators.jl) are available, they have limited time-frequency analysis, sonar analysis, and baseband analysis capabilities. This SignalAnalysis.jl package aims to fill that gap. The package has grown out of my research needs, but I hope to expand it over time to provide a wide variety of time-frequency analysis, and baseband signal processing tools.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"While this package works with most array-like data types, it uses the SignalBase.jl API to represent multichannel 1D signals (time, frequency or spatial domain). While the package adopts a SampledSignal data type to carry sampling rate information with the sampled signal, the API design allows sampling rate to be provided as a keyword argument in most cases, enabling the user to pass in any array-like data.","category":"page"},{"location":"index.html#APIs","page":"Home","title":"APIs","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Creating & managing signals\nGenerating signals\nBasic signal analysis\nSignal processing\nTime-frequency analysis\nArray signal processing\nRandom noise generation\nPlot recipes\nInteractive plotting","category":"page"},{"location":"index.html#Quick-links","page":"Home","title":"Quick links","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"SignalAnalysis.jl\nSignalBase.jl\nSampledSignals.jl\nDSP.jl\nPlots.jl\nUnitful.jl","category":"page"},{"location":"rand.html#Random-noise-generation","page":"Random noise generation","title":"Random noise generation","text":"","category":"section"},{"location":"rand.html","page":"Random noise generation","title":"Random noise generation","text":"PinkGaussian\nRedGaussian","category":"page"},{"location":"rand.html#SignalAnalysis.PinkGaussian","page":"Random noise generation","title":"SignalAnalysis.PinkGaussian","text":"struct PinkGaussian{T} <: Distributions.Distribution{Distributions.Multivariate,Distributions.Continuous}\n\nPink Gaussian noise distribution for random variate generation.\n\nExample:\n\njulia> rand(PinkGaussian(1000))\n1000-element Array{Float64,1}:\n   ⋮\n\njulia> rand(PinkGaussian(1000, 2.0))\n1000-element Array{Float64,1}:\n   ⋮\n\n\n\n\n\n","category":"type"},{"location":"rand.html#SignalAnalysis.RedGaussian","page":"Random noise generation","title":"SignalAnalysis.RedGaussian","text":"struct RedGaussian{T} <: Distributions.Distribution{Distributions.Multivariate,Distributions.Continuous}\n\nRed Gaussian noise distribution for random variate generation.\n\nExample:\n\njulia> rand(RedGaussian(1000))\n1000-element Array{Float64,1}:\n   ⋮\n\njulia> rand(RedGaussian(1000, 2.0))\n1000-element Array{Float64,1}:\n   ⋮\n\n\n\n\n\n","category":"type"}]
}