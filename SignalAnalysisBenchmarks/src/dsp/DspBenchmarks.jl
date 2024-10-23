module DspBenchmarks
    using SignalAnalysis: goertzel
    using BenchmarkTools
    
    const SUITE = BenchmarkGroup()
    SUITE["goertzel"] = BenchmarkGroup(["math","dsp"])
    # function goertzel(x::AbstractVector, f, n; fs=framerate(x))
    #   function goertzel(x::AbstractVector, f; fs=framerate(x))
    #   function goertzel(x::AbstractMatrix, f, n; fs=framerate(x))k
    goertzel_containers = [Vector,Matrix];
    datatypes = [Float64];
    #can feed a vector with a block size or not
    Fs = 34100;
    dims = [(Fs* duration,n_sensors) for n_sensors in 2:2:10 for duration in 1:2:10];
    vector_sizes = (1:2:10) .* Fs
    ns = [13 14 15].^2;
    Ts =[Float64,Float32,Float64];
    for n in ns
        for vector_size in vector_sizes
            for T ∈ Ts 
                V = rand(T,vector_size)
                SUITE["goertzel_vector"][string(vector_size),string(eltype(T))] = @benchmarkable goertzel($V,2200,$n;fs = $Fs)
            end
        end 
    end
end