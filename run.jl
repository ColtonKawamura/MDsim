include("src/GranMA.jl")

using .GranMA
using CellListMap

mass = 1
boxX = 1000
boxY = 5.6896
side = boxY
cutoff = 1.4 # max particle diameter
# box = Box([boxX,boxY], cutoff)
k = 100
gamma = 5

particleList = convertInput("2D_N5000_P0.1_Width5_Seed1.mat")

# cl = CellList([p.position for p in particleList], box)
box, cl = makeBox(boxX, boxY, cutoff, particleList) 
forceLaw = foo(k, gamma, box, cl)
VelInitial = VelocitiesInit(particleList, 0, 1)
trajectory = simAcou(particleList, VelInitial, 1, 0.001, 500, 10, forceLaw, side)
plotTrajectoryAcoustic(trajectory)

