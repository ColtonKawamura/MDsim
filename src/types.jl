export
    Pos2D,
    Pos3D,
    Particle

struct Pos2D{T} <: FieldVector{2,T}
    x::T
    y::T
end
 
struct Pos3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

mutable struct Particle{VecType}
    position::VecType
    diameter::Float64
end