using StaticArrays
using LinearAlgebra: norm
using Distributions
using Statistics
using Plots

struct Pos2D{T} <: FieldVector{2,T}
    x::T
    y::T
end
 
struct Pos3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

mutable struct Particle{VecType}
    position::VecType
    diameter::Float64
end

const k = 100

function main()
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,10], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, 1 ,1)
    trajectory = md_verlet(particle_list, VelInitial, mass, .05, 100, 1, forces!, ForceHooke)
    plot_trajectory(trajectory)
end

# Function to generate a random position within a given range
function PosRandom(::Type{VecType}, range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    return VecType(range[1] + rand(T)*(range[2]-range[1]) for _ in 1:dimensions)
end

# Function to generate a random particle with random position and diameter
function ParticleRandom(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end

function EnergyHooke(p_i::Particle{VecType}, p_j::Particle{VecType}) where VecType
    r_vector = p_j.position - p_i.position
    r = norm(r_vector)
    cutoff = .5 .* (p_i.diameter + p_i.diameter)
    if r > cutoff
        PotEnergy  = zero(T)
    else
        PotEnergy = .5 * k * (cutoff - r)^2 * (r_vector / r)
    end
    return PotEnergy
end


function ForceHooke(p_i::Particle{VecType}, p_j::Particle{VecType}) where VecType
    r_vector = p_j.position - p_i.position
    r = norm(r_vector)
    cutoff = 0.5 * (p_i.diameter + p_j.diameter)
    if r > cutoff
        force_i = zero(VecType)  # Ensure we're returning a vector zero
    else
        force_i = -k * (cutoff - r) * (r_vector / r)
    end
    return force_i
end

force_list = similar(map( p -> p.position, particle_list))
# position_list = map( p -> p.position, particle_list)

function forces!(force_list::Vector{VecType}, particle_list::Vector{Particle{VecType}}, ForceLaw::F) where {VecType, F}
    fill!(force_list, zero(VecType)) 
    n = length(particle_list)
    for i in 1:n-1
        for j in i+1:n
            force_i = ForceLaw(particle_list[i], particle_list[j])
            force_list[i] += force_i
            force_list[j] -= force_i
        end
    end
    return force_list
end

function VelocitiesInit(particle_list::Vector{Particle{VecType}}, temperature::Real, mass::Real) where VecType
    # Determine the number of dimensions based on the first particle's position
    dims = length(fieldnames(VecType))
    N = length(particle_list)  # Number of particles

    # Initialize velocities based on equipartition theorem
    initial_velocities = sqrt(temperature / mass) .* rand(N, dims)
    
    # Ensure net velocity is zero
    initial_velocities .-= mean(initial_velocities, dims=1)
    
    # Return the velocities as a vector of position-like objects
    if dims == 2
        VelInitial = map(i -> VecType(initial_velocities[i, 1], initial_velocities[i, 2]), 1:N)
    elseif dims == 3
        VelInitial =  map(i -> VecType(initial_velocities[i, 1], initial_velocities[i, 2], initial_velocities[i, 3]), 1:N)
    else
        error("Unsupported dimensionality: $dims. Only 2D and 3D are supported.")
    end
    return VelInitial
end


function md_verlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, ForceLaw) where {VecType}
    p = getfield.(particle_list, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    trajectory = [(map(p -> copy(p.position), particle_list), map(p -> p.diameter, particle_list))] 
 
    for step in 1:nsteps
        forces!(f, particle_list, ForceHooke)
        
        a .= f ./ mass

        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        # Compute new forces at updated positions
        setfield!.(particle_list, :position, p)
        forces!(f, particle_list, ForceHooke)

        a_new = f ./ mass

        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(p -> copy(p.position), particle_list), map(p -> p.diameter, particle_list)))
        end
    end

    return trajectory
end


function plot_trajectory(trajectory)
    initial_positions, initial_diameters = trajectory[1]
    xlims = (minimum([pos.x for pos in initial_positions]), maximum([pos.x for pos in initial_positions]))
    ylims = (minimum([pos.y for pos in initial_positions]), maximum([pos.y for pos in initial_positions]))

    marker_size_scale = 22  # Adjust this if the sizes don't match visually as expected.

    @gif for i in 1:length(trajectory)
        positions, diameters = trajectory[i]
        x = [pos.x for pos in positions]
        y = [pos.y for pos in positions]
        
        scatter(x, y, ms=diameters .* marker_size_scale, legend=false, xlims=xlims, ylims=ylims, markershape=:circle)
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end

function periodic(x,side)
    x = rem(x,side)
    if x >= side/2
        x -= side
    elseif x < -side/2
        x += side
    end
    return x
end