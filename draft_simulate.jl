using Random
using Plots
using LaTeXStrings
using Statistics
using BenchmarkTools
plotlyjs()

function main()
    ###### Stuff needed for packing
    N = 100
    W = 5
    diameter_average = 1
    diameter_spread = 0.4
    bi_disperse = true
    particle_diameters, initial_positions = create_packings_3D(N, W, diameter_average, diameter_spread, bi_disperse)
    
    ###### Plotting
    # plot_packing(initial_positions, particle_diameters)

    ###### Stuff needed for simulation
    K = 100
    temperature = 1
    mass = 1
    initial_velocities = sqrt(temperature / mass) .* rand(N, 3)  # Based on equipartition theorem
    initial_velocities .-= mean(initial_velocities, dims=1)  # Ensure net velocity is zero
    
    # Initialize the acceleration as zero
    initial_accelerations = zeros(N, 3)
    
    # Steps
    number_of_steps = Int(1E5)  # Convert to integer
    dt = pi * sqrt(mass / K) * 0.05
    
    # Initialize positions array
    positions = Array{Float64}(undef, N, 3, number_of_steps)
    
    # Run the Verlet integration
    positions = verlet_integration(N, K, particle_diameters, initial_positions, initial_velocities, initial_accelerations, number_of_steps, dt)
    
    # Return or save the final positions as needed
    return positions
end

function verlet_integration(N, K, particle_diameters, initial_positions, initial_velocities, initial_accelerations, number_of_steps, dt)
    # Initialize arrays
    positions = Array{Float64}(undef, N, 3, number_of_steps)
    current_positions = initial_positions
    velocities = initial_velocities
    accelerations = initial_accelerations
    
    for step = 1:number_of_steps
        # Update positions
        current_positions .+= velocities .* dt .+ 0.5 .* accelerations .* dt^2
        
        # Store the current positions
        positions[:, :, step] = current_positions
        
        # Update velocities and accelerations here (not shown in the code snippet)
        # Example:
        # velocities .+= 0.5 .* accelerations .* dt
        # Update accelerations based on forces here (e.g., based on interactions between particles)
    end
    
    return positions
end


function create_packings_3D(N,W, diameter_average, diameter_spread, bi_disperse)

    if bi_disperse
        diameter_large = diameter_average + diameter_spread
        diameter_small = diameter_average
    end

    # Calulate the heigh to of the simulation box and scale to the size of the largest particle
    x_limit_upper = W*diameter_large
    y_limit_upper = W*diameter_large
    z_limit_upper = N/(W^2)*diameter_large

    ##  Create a rank 3 tensor that is x_limit_upper by y_limit_upper by z_limit_upper where each point is the x,y,z coordinate that is evenly spaced by diameter_large 
    # Generate ranges for x, y, and z coordinates so they don't extend past the boundary
    x_range = collect(range(diameter_large/2, x_limit_upper-diameter_large/2, step=diameter_large)) 
    y_range = collect(range(diameter_large/2, y_limit_upper-diameter_large/2, step=diameter_large))
    z_range = collect(range(diameter_large/2, z_limit_upper-diameter_large/2, step=diameter_large))

    # Create 3D grid of coordinates
    x_coords = [x for x in x_range, y in y_range, z in z_range] # creates length(z_range) matricies that have dimensions length(x) by length(y) where each element is the x-position
    y_coords = [y for x in x_range, y in y_range, z in z_range] # puts in all the values for x (in x_range) into each element element of y (in y_range) and each one of those into each element of z (in z_range)
    z_coords = [z for x in x_range, y in y_range, z in z_range] 

    # concatanate into a matrix where each row is the x,y,z coordinate of each particle has 
    # reshape x_coods into rows (:) with one column (1)
    initial_positions = hcat(reshape(x_coords, :, 1), reshape(y_coords, :, 1), reshape(z_coords, :, 1))

    # Readjust the number of particles to those that actually fit in the box
    N_actual = size(initial_positions,1)
    
    ## Take the new N and then randomly assign a partcile diamter to reach between diameter_large and diameter_small
    particle_diameters = diameter_small .+ (diameter_large - diameter_small) .* rand(N_actual)

    return particle_diameters, initial_positions
end

function plot_packing(initial_positions, particle_diameters)
    # Create the 3D plot
    x = initial_positions[:, 1]
    y = initial_positions[:, 2]
    z = initial_positions[:, 3]
    marker_sizes = particle_diameters * 10  # Adjust the scaling factor if needed for better visualization

    # Readjust the number of particles to those that actually fit in the box
    N_actual = size(initial_positions,1)

    # Generate hover texts including particle number and coordinates
    hover_texts = ["Particle $i: (x, y, z) = ($(x[i]), $(y[i]), $(z[i]))" for i in 1:N_actual]
    plot_before_compression = scatter3d(x, y, z, markersize=marker_sizes, label="", xlabel="X", ylabel="Y", zlabel="Z", title="3D Scatter Plot of Particles", hover=hover_texts)
    display(plot_before_compression)
end



function othermain()
    ###### Stuff needed for packing
    N = 100
    W = 5
    diameter_average = 1
    diameter_spread = 0.4
    bi_disperse = true
    particle_diameters, initial_positions = create_packings_3D(N, W, diameter_average, diameter_spread, bi_disperse)
    
    ###### Plotting
    # plot_packing(initial_positions, particle_diameters)

    ###### Stuff needed for simulation
    K = 100
    temperature = 1
    mass = 1
    initial_velocities = sqrt(temperature / mass) .* rand(N, 3)
    initial_velocities .-= mean(initial_velocities, dims=1)
    
    # Convert initial positions and velocities to SVectors
    initial_positions_s = [SVector{3, Float64}(initial_positions[i, :]) for i in 1:N]
    initial_velocities_s = [SVector{3, Float64}(initial_velocities[i, :]) for i in 1:N]
    initial_accelerations_s = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:N]

    # Steps
    number_of_steps = Int(1E5)
    dt = pi * sqrt(mass / K) * 0.05
    
    # Run the Verlet integration, storing results in memory
    positions = otherverlet_integration(N, K, particle_diameters, initial_positions_s, initial_velocities_s, initial_accelerations_s, number_of_steps, dt)
    
    # Optional: Save the final positions to disk
    # jldsave("final_positions.jld2", "positions" => positions)
    
    return positions
end

function otherverlet_integration(N, K, particle_diameters, initial_positions_s, initial_velocities_s, initial_accelerations_s, number_of_steps, dt)
    # Initialize storage for positions across all time steps
    all_positions = Vector{Vector{SVector{3, Float64}}}(undef, number_of_steps)
    current_positions = initial_positions_s
    velocities = initial_velocities_s
    accelerations = initial_accelerations_s
    for step in 1:number_of_steps
        # Update positions
        new_positions = [current_positions[i] .+ velocities[i] .* dt .+ 0.5 .* accelerations[i] .* dt^2 for i in 1:N]
        
        # Store the positions for this time step
        all_positions[step] = copy(new_positions)
        
        # Update velocities and accelerations here (example, pseudocode):
        # forces = compute_forces(new_positions, particle_diameters, N, K)
        # new_accelerations = forces ./ mass
        # velocities .= velocities .+ 0.5 .* (accelerations .+ new_accelerations) .* dt
        
        # Update current positions, velocities, and accelerations for the next step
        current_positions = new_positions
        # accelerations = new_accelerations
    end
    
    return all_positions
end

M = [
SVector(1, 2),
SVector(1.1, 2.2),
SVector(1.11, 2.22)
]

n_steps = 1E5

function bench_svector(M, n_steps)
    for step in 1:n_steps
        delta_M = step / n_steps
        out =  M .+ (SVector(delta_M, delta_M),)
    end
end

M = [
    1.0  2.0;
    1.1  2.2;
    1.11 2.22
]
n_steps = 1E5

function bench_matrix(M, n_steps)
    for step in 1:n_steps
        delta_M = step / n_steps
        M = M .+ delta_M
    end
end
