include("src/GranMA.jl")

using .GranMA
using CellListMap
using MATLAB
using MAT

mat_data = matread("2D_N5000_P0.001_Width5_Seed1.mat")
boxY = mat_data["Ly"]
mass = 1
boxX = 1000
# boxY = 5.6896
side = boxY
cutoff = 1.4 # max particle diameter
k = 100
gamma = 0

particleList = convertInput("2D_N5000_P0.001_Width5_Seed1.mat")


leftWall = [particle.position.x <= particle.diameter / 1.9 for particle in particleList] # added this to ID leftwall
rightWall = [particle.position.x >= 1000 - particle.diameter / 2 for particle in particleList] # added this to ID leftwall 

A = .1
omega = 5
# moveLeftWall = (step, dt) -> A * sin(omega * step * dt)
tau = 1 / (omega / (2 * pi)) * 2   # Pulse duration (e.g., 10 cycles)
sigma = tau / sqrt(2 * log(2))   # Spread of the pulse
t_max = 0 * tau  # Start at maximum
gaussian_envelope = (t, t_max, sigma) -> exp(-(t - t_max)^2 / (sigma^2))
moveLeftWall = (step, dt) -> A * cos(omega * (step * dt - t_max)) * gaussian_envelope(step * dt, t_max, sigma)



box, cl = makeBox(boxX, boxY, cutoff, particleList) 
forceLaw = foo(k, gamma, box, cl)
VelInitial = VelocitiesInit(particleList, 0, 1)

trajectory = mdVerlet(particleList, VelInitial, 1, 0.01, 1000, 10, forceLaw, side, group1=leftWall, group2=rightWall, moveFunc = moveLeftWall)

plotTrajectoryAcoustic(trajectory)
A = 2
closestIndex = argmin(abs.([trajectory[1][1][i][1] for i in 1:length(trajectory[1][1])] .- A))
particleIndex = closestIndex
posX = map(step -> step[1][particleIndex][1], trajectory)
t = collect(1:length(posX))
plot(t,posX)

# trajector[2][1] = xy positions at timestep 2
# trajectory[3][2] = diameter at timestep 3
# trajectory[3][1][2] = xy position atiem timestep 3 of particle 2
# trajectory[3][1][2][1] = x position atiem timestep 3 of particle 2
# trajectory[step][xy,diam][particleIndex][x,y]

# Step 1: Define the new target x-positions
target_positions = [15.0, 30, 50, 70]  # Updated target x-positions

# Step 2: Find the indices of the particles closest to the target positions at step = 1
initial_x_positions = [trajectory[1][1][i][1] for i in 1:length(trajectory[1][1])]
closest_indices = [argmin(abs.(initial_x_positions .- target_position)) for target_position in target_positions]

# Step 3: Extract the x-positions of the closest particles over all time steps
x_all = [map(step -> step[1][closest_index][1], trajectory) for closest_index in closest_indices]

# Generate a time vector (assuming 1 unit per step)
time_vector = collect(1:length(x_all[1]))

# Step 4: Call MATLAB to plot the data, passing the target_positions variable
mat"""
x_all = $(x_all);
time_vector = $(time_vector);
target_positions = $(target_positions);  % Pass the target positions to MATLAB

% Calculate the maximum y-axis value for consistent scaling
max_y_value = 0; % Initialize
for i = 1:length(x_all)
    index_particle_to_plot = i;
    
    % Get max values for x positions
    max_y_value = max(max_y_value, max(abs(x_all{index_particle_to_plot} - mean(x_all{index_particle_to_plot}))));
end

% Set up tiled layout for the plots
tiledlayout(length(x_all), 1); % Create tiled layout with rows equal to the number of particles

% Loop over each selected particle to plot its position over time in separate tiles
for i = 1:length(x_all)
    nexttile;
    
    % Plot x position
    position_particles_x = x_all{i};
    plot(time_vector, position_particles_x - mean(position_particles_x), ...
        'DisplayName', sprintf('Particle closest to %.3f', target_positions(i)));
    
    % Set consistent y-axis limits
    ylim([-max_y_value, max_y_value]);
    
    % Box, labels, and legends for each tile
    box on;
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('X Position', 'Interpreter', 'latex');
    legend show;
    hold off;
end
"""
