include("src/GranMA.jl")

using .GranMA


function old_method()
    mass = 1
    k = 100
    side = 10
    particle_list = [ ParticleRandom(Pos2D{Float64}, [0,9], [.2,1]) for i in 1:100]
    VelInitial = VelocitiesInit(particle_list, .5 ,1)
    # trajectory = md_verlet(particle_list, VelInitial, mass, .05, 100, 1, forces!, ForceHooke)
    trajectory = md_verlet(particle_list, VelInitial, 1, 0.001, 1000, 10, 
           (force_list, particle_list) -> forces!(force_list, particle_list, (p_i, p_j) -> ForceHooke(p_i, p_j)), side)
    plot_trajectory(trajectory)
end