export
    md_verlet,
    md_verlet_Acoustic,
    md_verlet_AcousticCL,
    md,
    md_verletCL,
    md_verlet_AcousticCL_2out,
    md_verletCLosc,
    md_verletAc


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

# Acoustic simulation with oscillation. Works. No cell list so too slow.
function md_verletAc(particleList::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, ForceLaw, side) where {VecType}
    p = getfield.(particleList, :position)  # extract just the positions as initial positions
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)
    p0 = copy(p)

    trajectory = [(map(element -> copy(element.position), particleList), map(element -> element.diameter, particleList))]
    A = .1 # added to to feed into sin
    omega = 100 # added to feed into sin
    leftIndex = [particle.position.x < particle.diameter / 2 for particle in particleList] # added this to ID leftwall
 
 
    for step in 1:nsteps
        forces!(f, particleList, ForceLaw)
        
        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2
        p[leftIndex] .= [Pos2D(p0[i][1] + A * sin(omega * step * dt), p0[i][2]) for i in findall(leftIndex)] # added update positions left wall


        p .= periodic.(p, side)
        setfield!.(particleList, :position, p) # update positions
        forces!(f, particleList, ForceLaw) # update forces
        f[leftIndex] .= [Pos2D(0.0, 0.0) for _ in findall(leftIndex)] # added to zero out forces on left wall

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particleList), map(element -> element.diameter, particleList)))
        end
    end

    return trajectory
end



# periodic method
function md_verletCL(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side) where {VecType}
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



# periodic method
function md_verletCLosc(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side) where {VecType}
    p = getfield.(particle_list, :position) 
    p0 = copy(p) # added to get initial positions, used for the oscilation
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    A = .1 # added to to feed into sin
    omega = 10 # added to feed into sin

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 
    
    leftIndex = [particle.position.x < particle.diameter / 2 for particle in particle_list] # added this to ID leftwall
    rightIndex = [particle.position.x > 1000 - particle.diameter / 2 for particle in particle_list] # added this to ID leftwall 
    
    for step in 1:nsteps
        forces!(f, particle_list)

        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2
        p[leftIndex] .= [Pos2D(p0[i][1] + A * sin(omega * step * dt), p0[i][2]) for i in findall(leftIndex)] # added update positions left wall
        p[rightIndex] .= [Pos2D(p0[i][1], p0[i][2]) for i in findall(rightIndex)] # added update positions left wall

        p .= periodic.(p, side)
        setfield!.(particle_list, :position, p) 
        forces!(f, particle_list) 

        f[leftIndex] .= [Pos2D(0.0, 0.0) for _ in findall(leftIndex)] # added to zero out forces on left wall
        f[rightIndex] .= [Pos2D(0.0, 0.0) for _ in findall(rightIndex)] # added to zero out forces on left wall

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end