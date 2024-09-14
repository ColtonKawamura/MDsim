export
    plotTrajectory,
    plotTrajectoryAcoustic

function plotTrajectory(trajectory)
    # Extract initial positions and diameters
    initial_positions, initial_diameters = trajectory[1]
    
    # Compute the min and max x and y values from the initial positions
    x_min = minimum(p.x for p in initial_positions)
    x_max = maximum(p.x for p in initial_positions)
    y_min = minimum(p.y for p in initial_positions)
    y_max = maximum(p.y for p in initial_positions)
    
    # Set the x and y limits to be 1 unit beyond the initial min/max values
    xlims = (x_min - 1, x_max + 1)
    ylims = (y_min - 1, y_max + 1)
    
    # Generate the plot animation
    @gif for i in 1:length(trajectory)
        positions, diameters = trajectory[i]
        plot(; xlims=xlims, ylims=ylims, legend=false, aspect_ratio=:equal)

        # Plot each particle as a circle
        for (pos, diameter) in zip(positions, diameters)
            circle_shape = Shape([(pos.x + diameter/2 * cos(θ), pos.y + diameter/2 * sin(θ)) for θ in range(0, 2π, length=50)])
            plot!(circle_shape, lw=2, c=:blue, fillalpha=0)
        end
        
        # Annotate with the current step number
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end
    

# Examples
function plotTrajectoryAcoustic(trajectory)
    # Extract initial positions and diameters
    initial_positions, initial_diameters = trajectory[1]
    
    # Compute the min and max x and y values from the initial positions
    x_min = minimum(p.x for p in initial_positions)
    x_max = 10
    y_min = minimum(p.y for p in initial_positions)
    y_max = maximum(p.y for p in initial_positions)
    
    # Set the x and y limits to be 1 unit beyond the initial min/max values
    xlims = (x_min - 1, x_max + 1)
    ylims = (y_min - 1, y_max + 1)
    
    # Generate the plot animation
    @gif for i in 1:length(trajectory)
        positions, diameters = trajectory[i]
        plot(; xlims=xlims, ylims=ylims, legend=false, aspect_ratio=:equal)

        # Plot each particle as a circle
        for (pos, diameter) in zip(positions, diameters)
            circle_shape = Shape([(pos.x + diameter/2 * cos(θ), pos.y + diameter/2 * sin(θ)) for θ in range(0, 2π, length=50)])
            plot!(circle_shape, lw=2, c=:blue, fillalpha=0)
        end
        
        # Annotate with the current step number
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end

function plotTrajectoryAcoustic(trajectory_flow, trajectory_left_wall=nothing)
    # Extract initial positions and diameters of flow particles
    initial_positions_flow, initial_diameters_flow = trajectory_flow[1]
    
    # Compute the min and max x and y values from the initial positions
    x_min = minimum(p.x for p in initial_positions_flow)
    x_max = 10
    y_min = minimum(p.y for p in initial_positions_flow)
    y_max = maximum(p.y for p in initial_positions_flow)
    
    # Set the x and y limits to be 1 unit beyond the initial min/max values
    xlims = (x_min - 1, x_max + 1)
    ylims = (y_min - 1, y_max + 1)
    
    # Generate the plot animation
    @gif for i in 1:length(trajectory_flow)
        positions_flow, diameters_flow = trajectory_flow[i]
        plot(; xlims=xlims, ylims=ylims, legend=false, aspect_ratio=:equal)

        # Plot flow particles
        for (pos, diameter) in zip(positions_flow, diameters_flow)
            circle_shape = Shape([(pos.x + diameter/2 * cos(θ), pos.y + diameter/2 * sin(θ)) for θ in range(0, 2π, length=50)])
            plot!(circle_shape, lw=2, c=:blue, fillalpha=0)
        end

        # If left wall particle trajectory is provided, plot them in a different color
        if trajectory_left_wall != nothing
            positions_left_wall, diameters_left_wall = trajectory_left_wall[i]
            for (pos, diameter) in zip(positions_left_wall, diameters_left_wall)
                circle_shape = Shape([(pos.x + diameter/2 * cos(θ), pos.y + diameter/2 * sin(θ)) for θ in range(0, 2π, length=50)])
                plot!(circle_shape, lw=2, c=:red, fillalpha=0)  # Use red color for left wall particles
            end
        end
        
        # Annotate with the current step number
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end
