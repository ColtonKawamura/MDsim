include("src/GranMA.jl")

using .GranMA
using CellListMap

mass = 1
boxX = 1000
boxY = 5.6896
side = boxY
cutoff = 1.4 # max particle diameter
k = 100
gamma = 5

particleList = convertInput("2D_N5000_P0.1_Width5_Seed1.mat")

leftWall = [particle.position.x <= particle.diameter / 1.9 for particle in particleList] # added this to ID leftwall
rightWall = [particle.position.x >= 1000 - particle.diameter / 2 for particle in particleList] # added this to ID leftwall 

A = .1
omega = 5
moveLeftWall = (step, dt) -> A * sin(omega * step * dt)

box, cl = makeBox(boxX, boxY, cutoff, particleList) 
forceLaw = foo(k, gamma, box, cl)
VelInitial = VelocitiesInit(particleList, 0, 1)

trajectory = mdVerlet(particleList, VelInitial, 1, 0.001, 500, 10, forceLaw, side, group1=leftWall, group2=rightWall, moveFunc = moveLeftWall)

plotTrajectoryAcoustic(trajectory)

