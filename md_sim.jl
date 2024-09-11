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
    mass = 1
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    # trajectory = md_verlet(particle_list, VelInitial, mass, .05, 100, 1, forces!, ForceHooke)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.01, 1000, 10, 
           (force_list, particle_list) -> forces!(force_list, particle_list, (p_i, p_j) -> ForceHooke(p_i, p_j)), side)
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

particle_list = [ ParticleRandom(Pos2D{Float64}, [0,10], [.2,1]) for i in 1:100]
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

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
 
    for step in 1:nsteps
        forces!(f, particle_list, ForceLaw)
        
        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        setfield!.(particle_list, :position, p) # update positions
        forces!(f, particle_list, ForceLaw) # update forces

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end


function plot_trajectory(trajectory)
    initial_positions, initial_diameters = trajectory[1]
    xlims = (-5, 15) 
    ylims = (-5, 15)

    @gif for i in 1:length(trajectory)
        positions, diameters = trajectory[i]
        plot(; xlims=xlims, ylims=ylims, legend=false, aspect_ratio=:equal)

        for (pos, diameter) in zip(positions, diameters)
            circle_shape = Shape([(pos.x + diameter/2 * cos(θ), pos.y + diameter/2 * sin(θ)) for θ in range(0, 2π, length=50)])
            plot!(circle_shape, lw=0, c=:blue)
        end
        
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end

# periodic on all sides, square
function periodic(p::VecType, side::T) where {VecType<:FieldVector, T}
    return VecType(mod.(p, side))
end

#  periodic boundary only to the y-component (2nd component)
function periodic(p::VecType, side::T) where {VecType<:FieldVector, T}
    return VecType(p[1], mod(p[2], side), p[3:end]...)  # Apply mod only to the 2nd component (y)
end

# reflect x
function reflect(p::VecType, v::VecType, side::T) where {VecType<:FieldVector, T}
    # Reflect the x-component of the position and velocity if the particle hits the boundaries
    x_reflected = p[1]
    v_reflected = v[1]
    
    if x_reflected < 0
        x_reflected = -x_reflected  # Reflect at 0 boundary
        v_reflected = -v_reflected  # Reverse velocity in the x-direction
    elseif x_reflected > side
        x_reflected = 2*side - x_reflected  # Reflect at side boundary
        v_reflected = -v_reflected  # Reverse velocity in the x-direction
    end
    
    # Return the updated position and velocity vectors
    return VecType(x_reflected, p[2:end]...), VecType(v_reflected, v[2:end]...)
end

const side = 10

# Regular method ... need to be able to pass the periodic into just the force law?
function md_verlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side) where {VecType}
    p = getfield.(particle_list, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
 
    for step in 1:nsteps
        forces!(f, particle_list)
        
        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        setfield!.(particle_list, :position, p) # update positions
        forces!(f, particle_list) # update forces

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end


# periodic method
function md_verlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side) where {VecType}
    p = getfield.(particle_list, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
 
    for step in 1:nsteps
        forces!(f, particle_list)
        
        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        p .= periodic.(p, side)
        setfield!.(particle_list, :position, p) # update positions
        forces!(f, particle_list) # update forces

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end

# reflect method
function md_verlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side) where {VecType}
    p = getfield.(particle_list, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    force_list = similar(p)

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
 
    for step in 1:nsteps
        # Calculate forces on each particle using the provided ForceLaw
        forces!(force_list, particle_list)
        
        # Compute accelerations from forces
        a .= force_list ./ mass
        
        # Update positions using Verlet integration
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        # Apply reflection in accordance with boundary conditions
        p, v = reflect.(p, v, side)
        
        # Update the positions of the particles
        setfield!.(particle_list, :position, p)
        
        # Update velocities using Verlet integration
        v .= v .+ 0.5 .* a .* dt

        # Save the trajectory at specified intervals
        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end

