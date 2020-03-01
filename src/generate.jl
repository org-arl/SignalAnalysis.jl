export cw, chirp
export hanning, hamming, tukey, cosine, lanczos, triang, bartlett, gaussian, bartlett_hann, blackman, kaiser, dpss

"Generate a sinusoidal signal."
function cw(freq, duration; fs=deffs, phase=0.0, window=nothing, analytic=true)
    t = 0:1/freqQ(fs):timeQ(duration)
    x = exp.(2im*π*freqQ(freq).*t .+ phase)
    if window != nothing
        x .*= getwindow(window, length(t))
    end
    analytic ? x : real(x)
end

"Generate a frequency modulated chirp signal."
function chirp(freq1, freq2, duration; fs=deffs, shape=:linear, phase=0.0, window=nothing, analytic=true)
    freq1 = freqQ(freq1)
    freq2 = freqQ(freq2)
    duration = timeQ(duration)
    t = 0:1/freqQ(fs):duration
    if shape == :linear
        cby2 = (freq2-freq1)/duration/2.0
        x = exp.(2im*π .* (cby2.*t.^2 .+ freq1.*t) .+ phase)
    elseif shape == :hyperbolic
        s = freq2*duration/(freq2-freq1)
        x = exp.(-2im*π*s*freq1 .* log.(abs.(1.0 .- t./s)))
    else
        @assert false "Unknown chirp shape"
    end
    if window != nothing
        x .*= getwindow(window, length(t))
    end
    analytic ? x : real(x)
end

### utility functions

function getwindow(window, n)
    window == nothing && return ones(n)
    window isa Tuple && return window[1](n, window[2])
    window(n)
end
