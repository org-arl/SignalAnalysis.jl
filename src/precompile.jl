# toy load to force precompilation of typical use cases

function _precompile(x)
  length(x), size(x), nframes(x), framerate(x), nchannels(x), isanalytic(x)
  domain(x), time([1,2,3], x), time(1:3, x), samples(x), analytic(x)
  toframe(0.2𝓈, x), toframe([0.2𝓈, 201ms], x), toframe((0.2𝓈, 201ms), x), toframe(0.2:0.201s, x)
  padded(x, (10, 5); delay=1), samerateas(x)
  collect(partition(x, 5))
  ifrequency(x), energy(x), meantime(x), rmsduration(x), meanfrequency(x), rmsbandwidth(x)
  removedc(x), removedc!(x)
  hb = rand(eltype(x), 127)
  goertzel(x, 1000.0), filt(hb, x), filtfilt(hb, x)
  upconvert(downconvert(x, 4, 2000.0), 4, 2000.0)
  if x isa AbstractVector
    mfilter(x[1:10], x), findsignal(x[1:10], x)
    delay(x, 0.5), delay(x, 1), delay!(x, 1), delay!(x, 0.5), delay!(x, 0.5𝓈), delay(x, 0.5𝓈)
  end
end

for T ∈ (Float64, Float32, ComplexF64, ComplexF32)
  _precompile(signal(randn(T, 200), 1000))
  _precompile(signal(randn(T, (200, 2)), 1000))
end

cw(1000, 0.1, 44100)
cw(1000.0, 0.1, 44100.0)
chirp(100.0, 5000.0, 0.1, 44100.0)
chirp(100, 5000, 0.1, 44100)
mseq(3), gmseq(3)
fir(127, 100.0, 5000.0; fs=44100.0)
