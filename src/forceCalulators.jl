#forceCalulators

export
    forces!,
    exforces!,
    forces_cl!,
    forces_CL!,
    forcesDamped!

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


function forces_CL!(k, force_list::Vector{VecType}, particle_list::Vector{Particle{VecType}}, ForceLaw::F, box::Box, cl::CellList ) where {VecType, F}
    fill!(force_list, zero(VecType))
    cl = UpdateCellList!([p.position for p in particle_list], box, cl, parallel=false)
    map_pairwise!(
        (p_i, p_j, i , j, d2, force_list) -> ForceLaw(particle_list, p_i, p_j, k, i, j, d2, force_list, box),
        force_list, box, cl, parallel=false
    )
    return force_list
end
function forcesDamped!(k, gamma, force_list::Vector{VecType}, particle_list::Vector{Particle{VecType}}, velocity_list::Vector{VecType}, ForceLaw::F, box::Box, cl::CellList) where {VecType, F}
    fill!(force_list, zero(VecType))
    
    # Update cell list to ensure particles are mapped correctly
    cl = UpdateCellList!([p.position for p in particle_list], box, cl, parallel=false)
    
    # Apply force calculation for every pair of particles
    map_pairwise!(
        (p_i, p_j, i, j, d2, force_list) -> ForceLaw(particle_list, p_i, p_j, velocity_list[i], velocity_list[j], k, gamma, i, j, d2, force_list, box),
        force_list, box, cl, parallel=false
    )
    
    return force_list
end
