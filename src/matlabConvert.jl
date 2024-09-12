export
    convertInput,
    convertSplit

function convertInput(file_path::String)
    # Open the .mat file
    mat_data = matread(file_path)
    
    # Extract variables 'x', 'y', and 'Dn'
    x = mat_data["x"]
    y = mat_data["y"]
    Dn = mat_data["Dn"]
    
    # Check if x, y, and Dn have the same length
    num_particles = length(x)
    if length(y) != num_particles || length(Dn) != num_particles
        error("x, y, and Dn must have the same length")
    end
    
    # Create a list of Particle{Pos2D} instances
    particle_list = [Particle(Pos2D(x[i], y[i]), Dn[i]) for i in 1:num_particles]
    
    return particle_list
end

function convertSplit(file_path::String)
    # Open the .mat file and read the necessary variables
    mat_data = matread(file_path)
    
    # Extract variables 'x', 'y', 'Dn', and 'Lx'
    x = vec(mat_data["x"])  # Treat 'x' as a 1D array
    y = vec(mat_data["y"])  # Treat 'y' as a 1D array
    Dn = vec(mat_data["Dn"])  # Treat 'Dn' as a 1D array
    Lx = mat_data["Lx"]  # Read 'Lx' from the mat file
    
    # Check if x, y, and Dn have the same length
    num_particles = length(x)
    if length(y) != num_particles || length(Dn) != num_particles
        error("x, y, and Dn must have the same length")
    end
    
    # Create a list of Particle{Pos2D} instances
    particle_list = [Particle(Pos2D(x[i], y[i]), Dn[i]) for i in 1:num_particles]
    
    # Split particles into left wall, right wall, and flow groups based on x positions
    left_wall_list = findall(x .< Dn ./ 2)  # Particles near the left wall
    right_wall_list = findall(x .> (Lx .- Dn ./ 2))  # Particles near the right wall
    bulk_list = findall(.! (x .< Dn ./ 2 .|| x .> (Lx .- Dn ./ 2)))  # Particles in the bulk
    
    # Use the linear indices to extract particles
    particlesLeftWall = particle_list[left_wall_list]
    particlesRightWall = particle_list[right_wall_list]
    particlesFlow = particle_list[bulk_list]
    
    return particlesLeftWall, particlesRightWall, particlesFlow
end