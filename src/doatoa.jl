using FindPeaks1D: findpeaks1d
using ImageFiltering: findlocalmaxima
using Optim

export snapdetect, snapdoatoa

"""
$(SIGNATURES)
Envelope of signal `x`.
"""
envelope(x) = abs.(hilbert(x))

"""
$(SIGNATURES)
Detect time of arrivals (ToAs) of snaps at each sensor. The column vector of `data` is sensor data
with sampling rate `fs`. Percentile `p` as a threshold, minimum time distance between two snaps
`tdist` (in seconds), minimum time length of a snap `twidth` (in seconds) are defined for the
peak detection of each sensor recording.
"""
function snapdetect(data::AbstractMatrix, 
                    fs::Real; 
                    p::Real = 99.5, 
                    tdist::Real = 2e-3, 
                    twidth::Union{Nothing,Real} = nothing)
  # M. Chitre, S. Kuselan, and V. Pallayil, The Journal of the Acoustical Society of America 838–847, 2012. 
  snaps = Vector{Int}[]
  for i ∈ axes(data, 2)
    @views x = data[:,i]
    env = envelope(x)
    height = quantile(env, p / 100)
    width = !isnothing(twidth) ? round.(Int, twidth .* fs) : nothing
    distance = round(Int, tdist .* fs)
    crds, _ = findpeaks1d(env; height=height, distance=distance, width=width)
    push!(snaps, crds)
  end
  snaps
end

snapdetect(data::SignalAnalysis.SampledSignal; kwargs...) = snapdetect(samples(data), framerate(data); kwargs...)

"""
$(SIGNATURES)
Return DoA-ToA hough space of the detected snaps.
"""
function houghtransform(snaps::AbstractVector{Vector{Int}}, 
                        fs::Real, 
                        θ::AbstractArray, 
                        rxpos::AbstractArray, 
                        c::Real)
  numsensors = size(rxpos, 2)
  numbeams = size(θ, 1)
  samplesd = round.(Int, steering(rxpos, c, θ) .* fs)
  maxtimesamples = maximum(maximum.(snaps))
  vote_numtime = maxtimesamples + maximum(samplesd)
  votetype = numsensors < typemax(Int8) ? Int8 : Int16 # use Int8 or Int16 to reduce the memory
  votes = zeros(votetype, numbeams, vote_numtime)
  Γs = (0:vote_numtime-1)
  for i ∈ 1:numbeams
      for j ∈ 1:numsensors
          δs = snaps[j] .- samplesd[j,i]
          indices = δs .+ 1 
          filter!(x -> (x ≥ 1) && (x ≤ vote_numtime), indices)
          votes[i,indices] .+= one(votetype)
      end
  end
  votes, Γs
end

"""
$(SIGNATURES)
Estimate coarse DoA-ToA based on the detected snaps using Hough transform.
"""
function coarse_doatoas(snaps::AbstractVector{Vector{Int}}, 
                        fs::Real, 
                        θ::AbstractArray{T}, 
                        rxpos::AbstractArray, 
                        c::Real; 
                        minvotes::Union{Nothing,Int} = nothing) where {T}
  numsensors = size(rxpos, 2)
  isnothing(minvotes) && (minvotes = numsensors * 2 ÷ 3)
  votes, Γs = houghtransform(snaps, fs, θ, rxpos, c)  # hough space
  pkindices = findlocalmaxima(votes)
  doatoas = if size(θ, 2) == 1
    Tuple{T,Int}[]
  else
    Tuple{T,T,Int}[]
  end
  for pkindex ∈ pkindices
    if votes[pkindex] ≥ minvotes
      push!(doatoas, (θ[pkindex.I[1],:]..., Γs[pkindex.I[2]]))
    end
  end
  doatoas
end

"""
$(SIGNATURES)
Return snap ToA of each sensor, which is associated with the coarse DoA-Toa `doatoa`. If the smallest 
discrepancy between the coarse ToAs and steering delay is larger than 5 samples, ToA = -1 which will be 
omitted from DoA-ToA refinement.
"""
function gather_snap(doatoa::Union{Tuple{T,Int},Tuple{T,T,Int}}, 
                     snaps::AbstractVector{Vector{Int}}, 
                     fs::Real, 
                     rxpos::AbstractArray, 
                     c::Real) where {T}
  numsensors = size(rxpos, 2)
  m = length(doatoa)
  θ′, Γ′ = [doatoa[1:m-1]...], last(doatoa)
  m == 3 && (θ′ = transpose(θ′))
  samplesd = steering(rxpos, c, θ′) .* fs 
  associated_snaps = Vector{Int}()
  for i ∈ 1:numsensors
    @views snap = snaps[i]
    @views minval, index = findmin(abs.(snap .- (Γ′ + samplesd[i,1])))
    if minval ≤ 5 # error in samples
      push!(associated_snaps, snap[index])
    else
      push!(associated_snaps, -1)
    end
  end
  associated_snaps
end

"""
$(SIGNATURES)
Refine DoA-ToA of the coarse estimate `doatoa`.
"""
function refine_doatoa(doatoa::Union{Tuple{T,Int},Tuple{T,T,Int}}, 
                       snaps::AbstractVector{Vector{Int}}, 
                       fs::Real, 
                       rxpos::AbstractArray, 
                       c::Real) where {T}
  anglebnd = deg2rad(T(2))
  samplebnd = T(5)
  associated_snaps = gather_snap(doatoa, snaps, fs, rxpos, c)
  b = associated_snaps .> 0
  # loss function
  function mse1d(ζ)
    mean(abs2, associated_snaps[b] .- (last(ζ) .+ steering(rxpos, c, ζ[1:1])[b,1] .* fs))
  end
  function mse2d(ζ)
    mean(abs2, associated_snaps[b] .- (last(ζ) .+ steering(rxpos, c, transpose(ζ[1:2]))[b,1] .* fs))
  end
  # lower and upper bounds
  m = length(doatoa)
  is1d = m == 2 ? true : false
  if is1d
    lower = [doatoa[1] - anglebnd, T(doatoa[2] - 2)]
    upper = [doatoa[1] + anglebnd, T(doatoa[2] + 2)]
    mse = mse1d
  else
    lower = [doatoa[1] - anglebnd, doatoa[2] - anglebnd, T(doatoa[3] - samplebnd)]
    upper = [doatoa[1] + anglebnd, doatoa[2] + anglebnd, T(doatoa[3] + samplebnd)]
    mse = mse2d
  end
  # refine DoA-ToA
  doatoa1 = [T(doatoa1) for doatoa1 ∈ doatoa] # vector of 3 elements
  result = optimize(mse, lower, upper, doatoa1, Fminbox(LBFGS()), Optim.Options(g_tol = 1e-9))
  if (Optim.minimum(result) < 1) # mse is smaller than one sample.  
    return NTuple{m,T}(Optim.minimizer(result))
  else
    return nothing
  end
end

"""
$(SIGNATURES)
Refine DoA-ToAs of the coarse estimates `doatoas`.
"""
function refine_doatoas(doatoas::Union{Vector{Tuple{T,Int}},Vector{Tuple{T,T,Int}}}, 
                        snaps::AbstractVector{Vector{Int}}, 
                        fs::Real, 
                        rxpos::AbstractArray, 
                        c::Real) where {T}
  m = eltype(doatoas) |> fieldcount
  refinedoatoas = NTuple{m,T}[]
  for doatoa ∈ doatoas
    refinedoatoa = refine_doatoa(doatoa, snaps, fs, rxpos, c)
    !isnothing(refinedoatoa) && push!(refinedoatoas, refinedoatoa)
  end
  refinedoatoas
end

"""
$(SIGNATURES)
Detect and estimate direction of arrivals and time of arrivals (DoAs-ToAs) of snaps from an array
of acoustic recordngs `data` with sampling rate `fs`. The DoA-ToAs output is a Vector of Tuple
(azimuth, elevation, ToA). Sensor positions `rxpos` and all the possible coarse direction of
arrivals `θ` are required for DoA-Toa estimation. By default, the sound speed `c` is 1540.

Percentile `p` as a threshold, minimum time distance (sec) between two snaps `tdist`, minimum time
length (sec) of a snap `twidth` are defined for the peak detection of each sensor recording.

On the hough space, only DoA-ToAs with number of votes larger than or equal to `minvotes` are
selected for the further processing.

# Example
```julia-repl
julia> data = randn(10000, 4)
10000×4 Matrix{Float64}:
 -1.35577   -0.133958   1.83179    0.779208
 -0.674991   0.085191   1.45824   -0.245136
  0.366752   0.398075   0.144273   0.46261
  1.38094    0.986693  -0.466958  -0.306699
  ⋮                               
  0.682338   0.741276  -1.3476    -0.355784
  0.11414    0.686721   0.256951  -0.936501
 -1.07116   -0.570336   1.26326    1.23492
  0.381281  -1.62439    0.348399   0.176433

julia> fs = 200000
200000

julia> rxpos = [0.0 1.0 0.0 0.0; 
                0.0 0.0 1.0 0.0;
                0.0 0.0 0.0 1.0]
3×4 Matrix{Float64}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> c = 1540.0
1540.0
  
julia> p = 99.999
99.999
  
julia> tdist = 2e-3
0.002
  
julia> minvotes = 4
4

julia> Γ₀ = 1000
1000

julia> data[Γ₀ - round(Int, 0.75 / c * fs),2] += 100 
99.05055855900952

julia> data[Γ₀ + round(Int, 0.25 / c * fs), [1,3,4]] .+= 100
3-element view(::Matrix{Float64}, 1130, [1, 3, 4]) with eltype Float64:
 100.53221300670761
  99.77840685117116
  99.8155047017162

julia> snaps = snapdetect(data, fs; p = p, tdist = tdist)
4-element Vector{Vector{Int64}}:
 [1032]
 [903]
 [1032]
 [1032]

julia> doatoas = snapdoatoa(data, fs, θ, rxpos, c; p = p, tdist = tdist, minvotes = minvotes) # (azimuth, elevation, ToA)
1-element Vector{Tuple{Float64, Float64, Float64}}:
 (0.00336680173802057, 0.0033667826560101035, 999.7500000000016) 
```
"""
function snapdoatoa(data::AbstractMatrix, 
                    fs::Real, 
                    θ::AbstractArray{T}, 
                    rxpos::AbstractArray, 
                    c::Real = 1540.0; 
                    p::Real = 99.5, 
                    tdist::Real = 2e-3, 
                    twidth::Union{Nothing,Real} = nothing,
                    minvotes::Union{Nothing,Int} = nothing) where {T}
  # M. Chitre, S. Kuselan, and V. Pallayil, The Journal of the Acoustical Society of America 838–847, 2012. 
  snaps = snapdetect(data, fs; p = p, tdist = tdist, twidth = twidth)
  doatoas = coarse_doatoas(snaps, fs, θ, rxpos, c; minvotes=minvotes)
  refine_doatoas(doatoas, snaps, fs, rxpos, c)
end

snapdoatoa(data::SignalAnalysis.SampledSignal, args...; kwargs...) = snapdoatoa(samples(data), framerate(data), args...; kwargs...)