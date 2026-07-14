module WAVExt

using SignalAnalysis
using WAV: wavread

function loadsignal(filename::AbstractString; start=1, nsamples=missing)
  if start == 1 && nsamples === missing
    data, fs = wavread(filename)
  elseif start == 1
    data, fs = wavread(filename; subrange=nsamples)
  else
    if nsamples === missing
      n, ch = wavread(filename; format="size")
      nsamples = n - start + 1
    end
    data, fs = wavread(filename; subrange=start:(start+nsamples-1))
  end
  signal(data, fs)
end

end # module
