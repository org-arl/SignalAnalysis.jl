import Unitful: Quantity, ustrip, upreferred, 𝐓
import Base: (:)

const TimeQuantity = Quantity{T,𝐓,S} where {T,S}
const FreqQuantity = Quantity{T,𝐓^-1,S} where {T,S}

timeQ(q) = float(q)
timeQ(q::TimeQuantity) = float(ustrip(upreferred(q)))

freqQ(q) = float(q)
freqQ(q::FreqQuantity) = float(ustrip(upreferred(q)))

Base.float(q::TimeQuantity) = float(ustrip(upreferred(q)))
Base.float(q::FreqQuantity) = float(ustrip(upreferred(q)))

Base.getindex(s::AbstractArray, t::Quantity{T,𝐓,U}, ndx...) where {T,U} = s[toindex(t),ndx...]

(:)(start::TimeQuantity, stop::TimeQuantity) = toindex(start):toindex(stop)
(:)(start::Int, stop::TimeQuantity) = start:toindex(stop)
(:)(start::TimeQuantity, stop::Int) = toindex(start):stop
(:)(start::TimeQuantity, step::Int, stop::TimeQuantity) = toindex(start):step:toindex(stop)
(:)(start::Int, step::Int, stop::TimeQuantity) = start:step:toindex(stop)
(:)(start::TimeQuantity, step::Int, stop::Int) = toindex(start):step:stop
