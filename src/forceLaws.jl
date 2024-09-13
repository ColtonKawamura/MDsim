export
    EnergyHooke,
    ForceHooke


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