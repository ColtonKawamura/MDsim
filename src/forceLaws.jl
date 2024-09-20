export
    EnergyHooke,
    ForceHooke,
    ForceHookeCL,
    ForceHookeDamped,
    forceHookeDamp


function EnergyHooke(p_i::Particle{VecType}, p_j::Particle{VecType}, k) where VecType
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


function ForceHooke(p_i::Particle{VecType}, p_j::Particle{VecType}, k) where VecType
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

function ForceHookeCL(particle_list::Vector{Particle{VecType}}, p_i, p_j, k, i, j, d2, f::Vector{VecType}, box::Box) where {VecType}
    r_vector = p_j - p_i
    r = sqrt(d2)
    cutoff = 0.5 * (particle_list[i].diameter + particle_list[j].diameter)
    
    if r > cutoff
        force_i = zero(VecType)  # Ensure we're returning a zero vector of the correct dimension (2D or 3D)
    else
        force_i = -k * (cutoff - r) * (r_vector / r)
    end
    
    f[i] += force_i
    f[j] -= force_i
    return f
end
function ForceHookeDamped(particle_list::Vector{Particle{VecType}}, p_i, p_j, v_i::VecType, v_j::VecType, k, gamma, i, j, d2, f::Vector{VecType}, box::Box) where {VecType}
    r_vector = p_j - p_i
    r = sqrt(d2)
    cutoff = 0.5 * (particle_list[i].diameter + particle_list[j].diameter)
    
    # Hookean spring force
    if r > cutoff
        force_spring = zero(VecType)  # No force if beyond cutoff
    else
        force_spring = -k * (cutoff - r) * (r_vector / r)  # Spring force
    end
    
    # Relative velocity between particles
    v_rel = v_j - v_i
    
    # Damping force, proportional to relative velocity
    if r > cutoff
        force_damping = zero(VecType)  
    else
        # force_damping = gamma * StaticArrays.dot(v_rel, r_vector) / r * (r_vector / r)
        force_damping = gamma * v_rel # not negative because negative is captured by v_rel
    end
    
    # Total force
    force_i = force_spring + force_damping

    # Apply forces to both particles
    f[i] += force_i
    f[j] -= force_i
    
    return f
end

# Create a function that returns the Hooke damped force function with k and gamma set
function forceHookeDamp(k,gamma, box, cl)
    return (forceList, particleList, v) -> forcesDamped!(k, gamma, forceList, particleList, v, ForceHookeDamped, box, cl) 
end