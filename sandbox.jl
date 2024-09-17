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


function periodicY()
    mass = 1
    k = 100
    side = 5
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.5,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.001, 1000, 10, forces!, (p_i, p_j) -> ForceHooke(p_i, p_j, k), side)
    plotTrajectory(trajectory)
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
    trajectory = md_verletCL(particle_list, VelInitial, 1, 0.001, 1000, 10, (f_flow, all_particles) -> forces_CL!(k, f_flow, all_particles, ForceHookeCL, box, cl), side)
    plotTrajectory(trajectory)
end

# Acoustic simulation with oscillation. Works. No cell list so too slow.
function Acou()
    mass = 1
    k = 100
    side = 5
    particleList = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:50]
    # particleList = convertInput("2D_N5000_P0.1_Width5_Seed1.mat")
    VelInitial = VelocitiesInit(particleList, 1 ,1)
    trajectory = md_verletAc(particleList, VelInitial, 1, 0.001, 1000, 10, forces!, (p_i, p_j) -> ForceHooke(p_i, p_j, k), side)
    plotTrajectoryAcoustic(trajectory)
end


# Oscillaitng and plotting. Doens't look like left ins interacting with flow
function cellAcoust()
    mass = 1
    boxX = 1000
    # boxY = 5.6896
    boxY = 6.9205
    side = boxY
    cutoff = 1.4 # max particle diameter
    box = Box([boxX,boxY],cutoff)
    k = 100
    
    particleList = convertInput("2D_N5000_P0.001_Width5_Seed1.mat")

    cl = CellList([p.position for p in particleList], box)

    
    VelInitial = VelocitiesInit(particleList, 0, 1)
    trajectory = md_verletCLosc(particleList, VelInitial, 1, 0.001, 500, 10, (forceList, particleList) -> forces_CL!(k, forceList, particleList, ForceHookeCL, box, cl), side)
    plotTrajectoryAcoustic(trajectory)
end

function general()
    mass = 1
    boxX = 10
    boxY = 10
    side = boxY
    cutoff = .75# max particle diameter
    box = Box([boxX,boxY],cutoff)
    k = 100
    forceLaw = (forceList, particleList) -> forces_CL!(k, forceList, particleList, ForceHookeCL, box, cl)
    
    # particleList = convertInput("2D_N5000_P0.1_Width5_Seed1.mat")
    particleList = [ ParticleRandom(Pos2D{Float64}, [0,9], [.5,1]) for i in 1:50]

    cl = CellList([p.position for p in particleList], box)

    
    VelInitial = VelocitiesInit(particleList, 0, 1)
    trajectory = sim(particleList, VelInitial, 1, 0.001, 2000, 10, forceLaw)
    plotTrajectoryAcoustic(trajectory)
end

