### interface functions

export isanalytic, samplingrate, sampledsignal

"Create a signal from samples, converting to/from analytic representation, if needed."
function sampledsignal(s::AbstractArray, samplingrate=1.0; starttime=0.0, analytic=nothing)
    t = starttime .+ (0.0:1/samplingrate:(size(s,1)-1)/samplingrate)
    sampledsignal(AxisArray(s, Axis{:time}(t)); analytic=analytic)
end

function sampledsignal(s::AxisArray; analytic=nothing)
    analytic == true && !isanalytic(s) && return AxisArray(hilbert(s), s.axes...)
    analytic == false && isanalytic(s) && return AxisArray(real(s), s.axes...)
    return s
end

"Check if signal is analytic."
isanalytic(s) = eltype(s) <: Complex

"Get sampling rate of signal."
samplingrate(s) = 1.0
samplingrate(s::AxisArray) = 1.0/Float64(s.axes[1].val.step)

"Get time vector corresponding to each sample in signal."
Base.time(s::AbstractVector) = float(0:size(s,1)-1)
Base.time(s::AxisArray) = s.axes[1].val
