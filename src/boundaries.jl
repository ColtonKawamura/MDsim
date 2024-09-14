export
    periodic,
    reflect,
    wrap


# periodic on all sides, square
function periodic(p::VecType, side::T) where {VecType<:FieldVector, T}
    return VecType(mod.(p, side))
end

#  method 2 with just y boundary...need to add way to flag
function periodic(p::VecType, side::T) where {VecType<:FieldVector, T}
    return VecType(p[1], mod(p[2], side), p[3:end]...)  # Apply mod only to the 2nd component (y)
end

# reflect x
function reflect(p::VecType, v::VecType, side::T) where {VecType<:FieldVector, T}
    # Reflect the x-component of the position and velocity if the particle hits the boundaries
    x_reflected = p[1]
    v_reflected = v[1]
    
    if x_reflected < 0
        x_reflected = -x_reflected  # Reflect at 0 boundary
        v_reflected = -v_reflected  # Reverse velocity in the x-direction
    elseif x_reflected > side
        x_reflected = 2*side - x_reflected  # Reflect at side boundary
        v_reflected = -v_reflected  # Reverse velocity in the x-direction
    end
    
    # Return the updated position and velocity vectors
    return VecType(x_reflected, p[2:end]...), VecType(v_reflected, v[2:end]...)
end

#Examples
function wrap(x,side)
    x = rem(x,side)
    if x >= side/2
        x -= side
    elseif x < -side/2
        x += side
    end
    return x
end