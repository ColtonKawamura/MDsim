module GranMA

using StaticArrays
using LinearAlgebra: norm
using Distributions
using Statistics
using Plots
using MAT

include("types.jl")
include("makeParticles.jl")
include("forceLaws.jl")
include("forceCalulators.jl")
include("velocities.jl")
include("simulate.jl")
include("visualize.jl")
include("boundaries.jl")
include("matlabConvert.jl")

end