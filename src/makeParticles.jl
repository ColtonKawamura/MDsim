export 
    PosRandom,
    ParticleRandom,
    PosLine,
    MakeParticlesAlongLine,
    MakeParticle



# Function to generate a random position within a given range
function PosRandom(::Type{VecType}, range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    return VecType(range[1] + rand(T)*(range[2]-range[1]) for _ in 1:dimensions)
end

# Function to generate a random particle with random position and diameter
function ParticleRandom(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end

# Function to generate evenly spaced positions along a line at x_fixed
function PosLine(::Type{VecType}, x_fixed::Float64, y_range::Vector{Float64}, num_particles::Int) where VecType
    # Generate evenly spaced y positions within the y_range
    y_positions = range(y_range[1], y_range[2], length=num_particles)
    
    # Create a list of positions with fixed x and evenly spaced y
    positions = [VecType(x_fixed, y) for y in y_positions]
    
    return positions
end

# Function to create particles with evenly spaced positions along the line at x_fixed
function MakeParticlesAlongLine(::Type{VecType}, x_fixed::Float64, y_range::Vector{Float64}, DiamRange::Vector{Float64}, num_particles::Int) where VecType
    # Get evenly spaced positions along the line
    positions = PosLine(VecType, x_fixed, y_range, num_particles)
    
    # Create particles with these positions and random diameters
    particles = [Particle(position, DiamRange[1] + rand() * (DiamRange[2] - DiamRange[1])) for position in positions]
    
    return particles
end

function MakeParticle(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end