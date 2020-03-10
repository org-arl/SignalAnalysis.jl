export fir, removedc, removedc!

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
  # based on https://www.embedded.com/dsp-tricks-dc-removal/ real-time DC removal filter (d)
  s[2:end] .+= (α-1)*s[1:end-1]
  return s
end

"DC removal filter."
removedc(s; α=0.95) = removedc!(copy(s); α=α)

