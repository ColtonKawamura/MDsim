export
    mdVerlet

function mdVerlet(particle_list::Vector{Particle{VecType}}, VelInitial::Vector{VecType}, mass, dt, nsteps, save_interval, forces!, side; kwargs...) where {VecType}
    p = getfield.(particle_list, :position) 
    p0 = copy(p)  # Added to get initial positions, used for the oscillation
    v = copy(VelInitial)
    a = similar(p)
    f = similar(p)

    trajectory = [(map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list))] 

    # kwargs for groups
    group1 = get(kwargs, :group1, nothing)
    group2 = get(kwargs, :group2, nothing)

    # Retrieve the move function for group1
    moveFunc = get(kwargs, :moveFunc, (step, dt) -> 0.0)  # Default to no movement if moveFunc not passed

    for step in 1:nsteps
        forces!(f, particle_list, v)

        a .= f ./ mass
        p .= p .+ v .* dt .+ 0.5 .* a .* dt^2

        # Move group1 based on the moveFunc
        movement = moveFunc(step, dt)
        p[group1] .= [Pos2D(p0[i][1] + movement, p0[i][2]) for i in findall(group1)]  # Move group1 according to moveFunc

        p[group2] .= [Pos2D(p0[i][1], p0[i][2]) for i in findall(group2)]  # Keep group2 fixed

        p .= periodic.(p, side)
        setfield!.(particle_list, :position, p)
        forces!(f, particle_list, v)

        f[group1] .= [Pos2D(0.0, 0.0) for _ in findall(group1)]  # Zero out forces on group1
        f[group2] .= [Pos2D(0.0, 0.0) for _ in findall(group2)]  # Zero out forces on group2

        a_new = f ./ mass
        v .= v .+ 0.5 .* (a .+ a_new) .* dt

        if mod(step, save_interval) == 0
            println("Saved trajectory at step: ", step)
            push!(trajectory, (map(element -> copy(element.position), particle_list), map(element -> element.diameter, particle_list)))
        end
    end

    return trajectory
end
