include("src/GranMA.jl")

using .GranMA

const side = 10

# Works, but should probably change ForceHooke so that it takes k = 100 (current inside GranMA)
function periodicY()
    mass = 1
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.001, 1000, 10, forces!, ForceHooke, side)
    plotTrajectory(trajectory)
end

function noPeriodic()
    mass = 1
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.001, 1000, 10, forces!, ForceHooke)
    plotTrajectory(trajectory)
end

function acoustics()
    mass = 1
    particlesLeftWall, particlesRightWall, particlesFlow = convertSplit("2D_N5000_P0.1_Width5_Seed1.mat")
    VelInitial = VelocitiesInit(particlesFlow, 1, 1)
    trajectory = md_verlet_Acoustic(particlesFlow, particlesLeftWall, particlesRightWall, VelInitial, 1, .1, 200, 10, forces!, ForceHooke, 5)
    plotTrajectory(trajectory)
end