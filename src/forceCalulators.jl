#forceCalulators

export
    forces!,
    fPair_cl,
    exforces!,
    forces_cl!

# function forces!(force_list::Vector{VecType}, particle_list::Vector{Particle{VecType}}, ForceLaw::F) where {VecType, F}
#     fill!(force_list, zero(VecType)) 
#     n = length(particle_list)
#     for i in 1:n-1
#         for j in i+1:n
#             force_i = ForceLaw(particle_list[i], particle_list[j])
#             force_list[i] += force_i
#             force_list[j] -= force_i
#         end
#     end
#     return force_list
# end



function fPair_cl(p_i,p_j,i,j,d2,f,box::Box)
    Δv = p_j - p_i
    d = sqrt(d2)
    fₓ = 2*(d - box.cutoff)*(Δv/d)
    f[i] += fₓ
    f[j] -= fₓ
    return f
end


# Examples

function forces!(f::Vector{T},x::Vector{T},fₓ::F) where {T,F}
    fill!(f,zero(T))
    n = length(x)
    for i in 1:n-1
        for j in i+1:n
            fᵢ = fₓ(i,j,x[i],x[j])
            f[i] += fᵢ 
            f[j] -= fᵢ
        end
    end
    return f
end

function forces_cl!(f::Vector{T},x,box::Box,cl::CellList,fpair::F) where {T,F}
    fill!(f,zero(T))
    cl = UpdateCellList!(x,box,cl,parallel=false)
    map_pairwise!(
        (x,y,i,j,d2,f) -> fpair(x,y,i,j,d2,f,box),
        f, box, cl, parallel=false
    )
    return f
end