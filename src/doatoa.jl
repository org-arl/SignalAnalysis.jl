using FindPeaks1D: findpeaks1d
using ImageFiltering: findlocalmaxima
using Optim

export snapdetect, snapdoatoa

"""
Envelope of signal `x`.
"""
envelope(x) = abs.(hilbert(x))

"""
$(SIGNATURES)
Detect snaps from an array of acoustic recordings `data` with sampling rate `fs`. 
Percentile `p` as a threshold, minimum time distance between two snaps `tdist` (in seconds), 
minimum time length of a snap `twidth` (in seconds) are defined for the peak detection of each 
sensor recording.
"""
function snapdetect(data::AbstractMatrix{T}, 
                    fs::Real; 
                    p::Real = 99.5, 
                    tdist::Real = 2e-3, 
                    twidth::Union{Nothing,Real} = nothing) where {T}
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

function houghtransform(snaps, rxpos, θ, fs, c)
  numsensors = size(rxpos, 2)
  numbeams = size(θ, 1)
  samplesd = round.(Int, steering(rxpos, c, θ) .* fs)
  maxtimesamples = maximum(maximum.(snaps))
  vote_numtime = maxtimesamples + maximum(samplesd)
  votes = zeros(Int16, numbeams, vote_numtime) # use Int16 to reduce the memory
  Γs = (0:vote_numtime-1)
  for i ∈ 1:numbeams
      for j ∈ 1:numsensors
          δs = snaps[j] .- samplesd[j,i]
          indices = [δ + 1 for δ ∈ δs]
          filter!(x -> (x ≥ 0) && (x ≤ vote_numtime), indices)
          @views votes[i,indices] .+= one(Int16)
      end
  end
  votes, Γs
end

"""
$(SIGNATURES)
Estimate coarse DoA-ToA based on the detected snaps using Hough transform.
"""
function coarseDoAToAs(snaps, rxpos, θ, fs, c; minvotes::Union{Nothing,Int} = nothing)
  votes, Γs = houghtransform(snaps, rxpos, θ, fs, c)
  numsensors = size(rxpos, 2)
  isnothing(minvotes) && (minvotes = numsensors * 2 ÷ 3)

  pkindices = findlocalmaxima(votes)
  doatoas = if size(θ, 2) == 1
    Tuple{Float64,Int}[]
  else
    Tuple{Float64,Float64,Int}[]
  end
  for pkindex ∈ pkindices
    if votes[pkindex] ≥ minvotes
      θ1 = if size(θ, 2) == 1
        (θ[pkindex.I[1]], Γs[pkindex.I[2]])
      else
        (θ[pkindex.I[1],:]..., Γs[pkindex.I[2]])
      end
      push!(doatoas, θ1)
    end
  end
  doatoas
end

"""
Return snap ToA of each sensor, which is associated with the coarse DoA-Toa `doatoa`. 
If the smallest discrepancy between the coarse ToAs and steering delay is
larger than 5 samples, ToA = -1 which will be omitted for DoA-ToA refinement.
"""
function get_associatedsnap(doatoa::Tuple{T,T,Int}, snaps, rxpos, fs, c) where {T}
  numsensors = size(rxpos, 2)
  numdoatoa = length(doatoa)
  θ′, Γ′ = doatoa[1:numdoatoa-1], doatoa[numdoatoa]
  samplesd = round.(Int, steering(rxpos, c, [θ′...]') .* fs)
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
Refine DoA-ToA of the coarse estimate `doatoa`.
"""
function refineDoAToA(doatoa::Tuple{T,T,Int}, snaps, rxpos, fs, c) where {T}
  anglebnd = deg2rad(T(2))
  samplebnd = T(5)
  
  associated_snaps = get_associatedsnap(doatoa, snaps, rxpos, fs, c)
  b = associated_snaps .> 0
  function mse(ζ)
    numζ = length(ζ)
    mean(abs2, associated_snaps[b] .- (ζ[numζ] .+ steering(rxpos, c, ζ[1:numζ-1]')[b,1] .* fs))
  end
  if length(doatoa) == 2
    lower = [doatoa[1] - anglebnd, Float64(doatoa[2] - 2)]
    upper = [doatoa[1] + anglebnd, Float64(doatoa[2] + 2)]
  else
    lower = [doatoa[1] - anglebnd, doatoa[2] - anglebnd, T(doatoa[3] - samplebnd)]
    upper = [doatoa[1] + anglebnd, doatoa[2] + anglebnd, T(doatoa[3] + samplebnd)]
  end
  result = optimize(mse, lower, upper, [T(doatoa1) for doatoa1 ∈ doatoa], Fminbox(LBFGS()), Optim.Options(g_tol = 1e-9))
  if (result.minimum ≤ 1) # mse is smaller or equal to one sample.  
    return Tuple(result.minimizer)
  else
    return nothing
  end
end

"""
Refine DoA-ToAs of the coarse estimates `doatoas`.
"""
function refineDoAToAs(doatoas::Vector{Tuple{T,T,Int}}, snaps, rxpos, fs, c) where {T}
  refinedoatoas = NTuple{3,T}[]
  for doatoa ∈ doatoas
    refinedoatoa = refineDoAToA(doatoa, snaps, rxpos, fs, c)
    !isnothing(refinedoatoa) && push!(refinedoatoas, refinedoatoa)
  end
  refinedoatoas
end

"""
$(SIGNATURES)
Detect and estimate direction of arrivals and time of arrivals (DoA-ToA) of snaps from an array of acoustic 
recordngs `data` with sampling rate `fs`. The DoA-ToAs output is a Vector of Tuple (azimuth, elevation, ToA).
Sensor positions `rxpos` and all the possible coarse direction of arrivals `θ` are required for DoA-Toa estimation. 
By default, the sound speed `c` is 1540.

Percentile `p` as a threshold, minimum time distance (sec) between two snaps `tdist`, minimum time length (sec)
of a snap `twidth` are defined for the peak detection of each sensor recording.

On the hough space, only DoA-ToAs with number of votes larger than or equal to `minvotes` are selected for the 
further processing.

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

julia> doatoas = snapdoatoa(data, fs, rxpos, θ, c; p = p, tdist = tdist, minvotes = minvotes) # (azimuth, elevation, ToA)
1-element Vector{Tuple{Float64, Float64, Float64}}:
 (0.00336680173802057, 0.0033667826560101035, 999.7500000000016) 
```
"""
function snapdoatoa(data, 
                    fs, 
                    rxpos, 
                    θ, 
                    c::Real = 1540.0; 
                    p::Real = 99.5, 
                    tdist::Real = 2e-3, 
                    twidth::Union{Nothing,Real} = nothing,
                    minvotes::Union{Nothing,Int} = nothing)
  # M. Chitre, S. Kuselan, and V. Pallayil, The Journal of the Acoustical Society of America 838–847, 2012. 
  snaps = snapdetect(data, fs; p = p, tdist = tdist, twidth = twidth)
  doatoas = coarseDoAToAs(snaps, rxpos, θ, fs, c; minvotes=minvotes)
  refineDoAToAs(doatoas, snaps, rxpos, fs, c)
end
detect_doatoa(data::SignalAnalysis.SampledSignal, args...; kwargs...) = detect_doatoa(samples(data), framerate(data), args...; kwargs...)