export
    md_verlet,
    md_verlet_Acoustic


# method 1 without boundaries 
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

# periodic method
function md_verlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, ForceLaw, side) where {VecType}
    p = getfield.(particle_list, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
 
    for step in 1:nsteps
        forces!(f, particle_list, ForceLaw)
        
        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        p .= periodic.(p, side)
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

function md_verlet_Acoustic(flow_particles::Vector{Particle{VecType}}, 
    left_wall_particles::Vector{Particle{VecType}}, 
    right_wall_particles::Vector{Particle{VecType}}, 
    VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, ForceLaw, side) where {VecType}

    # Combine all particles (flow + wall particles) into a single list for force calculation
    all_particles = vcat(flow_particles, left_wall_particles, right_wall_particles)

    # Extract initial positions and velocities of flow particles
    p_flow = getfield.(flow_particles, :position)
    v_flow = copy(VelInitial)
    a_flow = similar(p_flow)
    f_flow = similar(getfield.(all_particles, :position))  # Create a force vector for all particles

    # Create masks for identifying which particles belong to the walls
    left_wall_list = [false for _ in 1:length(all_particles)]
    right_wall_list = [false for _ in 1:length(all_particles)]

    # Mark wall particles in the masks
    for i in 1:length(left_wall_particles)
        left_wall_list[i + length(flow_particles)] = true  # Mark left wall particles
    end
    for i in 1:length(right_wall_particles)
        right_wall_list[i + length(flow_particles) + length(left_wall_particles)] = true  # Mark right wall particles
    end

    # Initialize trajectory storage for the flow particles
    trajectory = [(map(p -> copy(p.position), flow_particles), map(p -> p.diameter, flow_particles))]

    for step in 1:nsteps
        # Calculate forces on all particles (including wall particles)
        forces!(f_flow, all_particles, ForceLaw)

        # Zero out the forces on the left and right wall particles so they don't move
        for i in 1:length(all_particles)
            if left_wall_list[i] || right_wall_list[i]
            f_flow[i] = zero(VecType)  # Set force to zero for wall particles
            end
        end

        # Update accelerations of flow particles only (no update for wall particles)
        a_flow .= f_flow[1:length(flow_particles)] ./ mass

        # Verlet integration: Update positions of flow particles only
        p_flow .= p_flow .+ v_flow .* dt .+ 0.5 .* a_flow .* dt^2

        # Set the new positions of the flow particles
        setfield!.(flow_particles, :position, p_flow)
        p_flow .= periodic.(p_flow, side)
        # Recalculate forces after position update
        forces!(f_flow, all_particles, ForceLaw)

        # Zero out forces on wall particles again to keep them stationary   
        for i in 1:length(all_particles)
            if left_wall_list[i] || right_wall_list[i]
                    f_flow[i] = zero(VecType)  # Set force to zero for wall particles
            end
        end

        # Update velocities of flow particles
        a_flow_new = f_flow[1:length(flow_particles)] ./ mass
        v_flow .= v_flow .+ 0.5 .* (a_flow .+ a_flow_new) .* dt

        # Save trajectory at intervals
        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(p -> copy(p.position), flow_particles), map(p -> p.diameter, flow_particles)))
        end 
    end

    return trajectory
end
