### plotting recipes

@recipe function plot(s::AxisArray; apply=nothing)
    ticks --> :native
    s1 = s
    slims = nothing
    if apply != nothing
        s1 = apply.(s.data)
        if apply == pow2db || apply == amp2db
            slims = (maximum(s1)-30, maximum(s1)+5)
        end
    end
    if ndims(s) == 2 && axisname(s.axes[1]) == "Frequency" && axisname(s.axes[2]) == "Time"
        xlabel --> "Time"
        ylabel --> "Frequency"
        if slims != nothing
            clims --> slims
        end
        @series begin
            seriestype := :heatmap
            s.axes[2].val, s.axes[1].val, real(s1)
        end
    else
        xlabel --> axisname(s.axes[1])
        legend --> ndims(s) > 1
        if slims != nothing
            ylims --> slims
        end
        @series begin
            seriestype := :line
            time(s), real(s1)
        end
        if s1 isa AbstractArray{<:Complex}
            @series begin
                seriestype := :line
                time(s), imag(s1)
            end
        end
    end
end

### utility functions

axisname(a::Axis{N,<:Any}) where N = uppercase(string(N)[1:1]) * string(N)[2:end]
