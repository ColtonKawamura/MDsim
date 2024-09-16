export
    EnergyHooke,
    ForceHooke,
    fₓ,
    fpair_cl,
    ForceHookeCL


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

function ForceHookeCL(particle_list, p_i, p_j, k, i , j , d2, f, box::Box) # using CellLists
    r_vector = p_j - p_i
    r = sqrt(d2)
    cutoff = 0.5 * (particle_list[i].diameter + particle_list[j].diameter)
    force_i = -k * (cutoff - r) * (r_vector / r)
    f[i] += force_i
    f[j] -= force_i
    return f
end

# Examples

function fₓ(x::T,y::T,cutoff,side) where T
    Δv = wrap.(y - x, side)
    d = norm(Δv)
    if d > cutoff
        fₓ = zero(T)
    else
        fₓ = 2*(d - cutoff)*(Δv/d)
    end
    return fₓ
end

function fpair_cl(x,y,i,j,d2,f,box::Box)
    Δv = y - x
    d = sqrt(d2)
    fₓ = 2*(d - box.cutoff)*(Δv/d)
    f[i] += fₓ
    f[j] -= fₓ
    return f
end

function ForceHookeCLDis(particle_list, p_i, p_j, k, i , j , d2, f, box::Box) # using CellLists
    r_vector = p_j - p_i
    r = sqrt(d2)
    cutoff = 0.5 * (particle_list[i].diameter + particle_list[j].diameter)
    force_i = -k * (cutoff - r) * (r_vector / r)
    f[i] += force_i
    f[j] -= force_i
    return f
end