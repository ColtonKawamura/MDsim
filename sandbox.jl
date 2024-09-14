include("src/GranMA.jl")

using .GranMA
using CellListMap


# Probably need to change this so it takes k
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

function periodicY()
    mass = 1
    k = 100
    side = 5
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.001, 1000, 10, forces!, (p_i, p_j) -> ForceHooke(p_i, p_j, k), side)
    plotTrajectory(trajectory)
end

function celllist()
    mass = 1
    boxX = 1000
    boxY = 5
    side = boxY
    cutoff = 1.4 # max particle diameter
    box = Box([boxX,boxY],cutoff)
    
    particlesLeftWall, particlesRightWall, particlesFlow = convertSplit("2D_N5000_P0.1_Width5_Seed1.mat")

    cl = CellList([p.position for p in particlesFlow], box)
    
    VelInitial = VelocitiesInit(particlesFlow, 1, 1)
    trajectory = md_verlet_AcousticCL(particlesFlow, particlesLeftWall, particlesRightWall, VelInitial, 1, .1, 200, 10, (f_flow, all_particles) -> forces_CL!(f_flow, all_particles, ForceHookeCL, box, cl), side)
    plotTrajectory(trajectory)
end


# Commented out m y forces! for theirs right now in forceCalculators.jl
function example()
    mass = 1
    boxX = 1000
    boxY = 5
    side = boxY
    cutoff = 1.4 # max particle diameter
    box = Box([boxX,boxY],cutoff)
    n_large = 1000
    
    particlesLeftWall, particlesRightWall, particlesFlow = convertSplit("2D_N5000_P0.1_Width5_Seed1.mat")
    particlesFlowPos = [p.position for p in particlesFlow]

    cl = CellList(particlesFlowPos,box)
    
    VelInitial = VelocitiesInit(particlesFlow, 1, 1)
    trajectory  = md((
        x0 = particlesFlowPos, 
        v0 = VelInitial, 
        mass = 1,
        dt = 0.1,
        nsteps = 1000,
        isave = 10,
        forces! = (f,x) -> forces_cl!(f,x,box,cl,fpair_cl)
    )...)
    plotTrajectoryEx(trajectory)
end
