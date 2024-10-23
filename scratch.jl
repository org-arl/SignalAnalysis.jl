using SignalAnalysis, BenchmarkTools
c = 340
rxpos = [-1.0 0.0 1.0 0.0; 0.0 1.0 0.0 -1.0]
fc = 500.0
fs = 44100.0
x = cw(fc, 30, fs)
x̄ = padded(x, 0; delay=round(Int, √2 / c * fs))
x4 = [x x x̄ x̄]
n = 4096;
θ = deg2rad.(1:360)
sd1 = steering(rxpos, c, θ)

music() = beamform(x4, fc, 4096, sd1; method=Music())
music();
@btime music();

capon() = beamform(x4, fc, 4096, sd1; method=Capon())
capon();
@btime capon();

go() = goertzel(x4, fc, n; fs=fs);
go();
@btime go();

bartlett() = beamform(x4, fc, 4096, sd1; method=Bartlett())
bartlett();
@btime bartlett()
