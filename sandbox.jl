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

# Oscillaitng and plotting. Doens't look like left ins interacting with flow
function cellList()
    mass = 1
    boxX = 1000
    boxY = 5
    side = boxY
    cutoff = 1.4 # max particle diameter
    box = Box([boxX,boxY],cutoff)
    k = 100
    
    particlesLeftWall, particlesRightWall, particlesFlow = convertSplit("2D_N5000_P0.1_Width5_Seed1.mat")
    all_particles = vcat(particlesFlow, particlesLeftWall, particlesRightWall)

    # cl = CellList([p.position for p in particlesFlow], box)
    # Combine flow and left wall particles in the cell list
    cl = CellList([p.position for p in vcat(particlesFlow, particlesLeftWall)], box)

    
    VelInitial = VelocitiesInit(all_particles, 1, 1)
    trajectory_flow, trajectory_left_wall = md_verlet_AcousticCL_2out(particlesFlow, particlesLeftWall, particlesRightWall, VelInitial, 1, .0001, 2000, 10, (f_flow, all_particles) -> forces_CL!(k, f_flow, all_particles, ForceHookeCL, box, cl), side)
    plotTrajectoryAcoustic(trajectory_flow, trajectory_left_wall)
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

# Workls try to incorprated boundaries
function periodicYCL()
    mass = 1
    k = 100
    side = 5
    boxX = 10
    boxY = 5
    side = boxY
    cutoff = 1 # max particle diameter
    box = Box([boxX,boxY],cutoff)

    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    cl = CellList([p.position for p in particle_list], box)

    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    trajectory = md_verletCL(particle_list, VelInitial, 1, 0.0001, 1000, 10, (f_flow, all_particles) -> forces_CL!(k, f_flow, all_particles, ForceHookeCL, box, cl), side)
    plotTrajectory(trajectory)
end



# Oscillaitng and plotting. Doens't look like left ins interacting with flow
function cellAcoust()
    mass = 1
    boxX = 1000
    boxY = 5
    side = boxY
    cutoff = 1.4 # max particle diameter
    box = Box([boxX,boxY],cutoff)
    k = 100
    
    particleList = convertInput("2D_N5000_P0.1_Width5_Seed1.mat")

    # cl = CellList([p.position for p in particlesFlow], box)
    # Combine flow and left wall particles in the cell list
    cl = CellList([p.position for p in particleList], box)

    
    VelInitial = VelocitiesInit(particleList, 1, 1)
    trajectory = md_verletCLosc(particleList, VelInitial, 1, 0.0001, 500, 10, (f_flow, particleList) -> forces_CL!(k, f_flow, particleList, ForceHookeCL, box, cl), side)
    plotTrajectoryAcoustic(trajectory)
end