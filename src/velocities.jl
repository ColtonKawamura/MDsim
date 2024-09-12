export
    VelocitiesInit
    
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