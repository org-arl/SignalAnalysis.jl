export fir, removedc, removedc!, demon

"Design FIR filter."
function fir(n, f1, f2=nothing; fs=deffs[], method=FIRWindow(hanning(n)))
  fs = freqQ(fs)
  if f1 == 0
    f = Lowpass(freqQ(f2); fs=fs)
  elseif f2 == nothing || freqQ(f2) == fs/2
    f = Highpass(freqQ(f1); fs=fs)
  else
    f = Bandpass(freqQ(f1), freqQ(f2); fs=fs)
  end
  return digitalfilter(f, method)
end

"DC removal filter."
function removedc!(s; α=0.95)
  # based on Lyons 2011 (3rd ed) real-time DC removal filter in Fig. 13-62(d)
  for k = 1:size(s,2)
    for j = 2:size(s,1)
      s[j,k] += α*s[j-1,k]
    end
    s[2:end,k] .+= -s[1:end-1,k]
  end
  s *= sqrt(α)
  return s
end

"DC removal filter."
removedc(s; α=0.95) = removedc!(copy(s); α=α)

"Estimate DEMON spectrum."
function demon(x; fs=deffs[], downsample=250, method=:rms, cutoff=1.0)
  local y
  fs /= downsample
  for k = 1:size(x,2)
    if downsample == 1
      y1 = abs.(hilbert(x[:,k]))
    elseif method == :rms
      y1 = sqrt.(mean.(Periodograms.arraysplit(abs2.(hilbert(x[:,k])), downsample, 0)))
    elseif method == :mean
      y1 = mean.(Periodograms.arraysplit(abs.(hilbert(x[:,k])), downsample, 0))
    elseif method == :fir
      aaf = fir(127, 0, 0.48fs; fs=fs)
      y1 = filtfilt(aaf, abs.(hilbert(x[:,k])))[1:downsample:end]
    else
      throw(ArgumentError("Unknown method"))
    end
    if k == 1
      y = zeros(length(y1), size(x,2))
    end
    y[:,k] .= y1
  end
  hpf = fir(127, cutoff; fs=fs)
  filtfilt(hpf, y)
end
