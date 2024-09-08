
# Granular Material Analysis (GranMA)

In an effort to learn more about particle simulations, I've taken notes go through [this](https://m3g.github.io/2021_FortranCon/) and taken notes. I have modified them to fit my specific needs for molecular dynamics.

## Creating Position Vector Structures

We're only dealing with with vectors that are either 2D or 3D for position. So, we can save significant computational overhead by using `FieldVector`'s in to define the spatial component of particles.

```julia
using StaticArrays

struct Pos2D{T} <: FieldVector{2,T}
    x::T
    y::T
end
```

We define a specific structure for our position vectors called `Pos2D` that is parameterized be a type `T`.  The type T could be anything, like Float64, Int, etc. The struct will contain two fields, `x` and `y`, which will both have the type `T`. The `<:` symbol indicates that `Vec2D` is a subtype of `FieldVector{2,T}`, which is an abstract type from `StaticArrays.jl`. We do this so that our `Pos2D` inherits  all the functionalities of a FieldVector, such as element access, iteration, and mathematical operations.

In practice, we simply feed the structure `Pos2D` a type `T` and an `x` and `y` that are type `T` and a vector is defined.

```julia
julia> Pos2D{Float64}(1,2)
2-element Pos2D{Float64} with indices SOneTo(2):
 1.0
 2.0
```
 
## Creating Random Particles

Next, for the sake of making examples, let's create a function that creates random vectors using the `Vec2D` structure.

```julia
function PosRandom(::Type{VecType}, range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    return VecType(range[1] + rand(T)*(range[2]-range[1]) for _ in 1:dimensions)
end
```


The arguments for our function `PosRandom` are:
* `::Type{VecType}` - the structure we want the output vector `VecType` to be
* `range`- the range of the points in that vector. 

`where VecType` introduces a type  `VecType` that can be used throughout the function to refer to the specific type passed as the first argument (`::Type{VecType}`). This allows us to use different vector types later on, for example 3D. 

***Important:*** `VecType` will be coming up alot in this tutorial, so remember it represents our choice of the `Pos2D{T}` structure.

Here's an example of creating a random vector with this function.

```julia
julia> PosRandom(Pos2D{Float64}, [0,1])
2-element Pos2D{Float64} with indices SOneTo(2):
 0.28591930991040737
 0.7042845413535875
```
We passed `Pos2D{Float64}` as the type of strucutre we want the output to be and `[0,1]`is the range of the random value for each element of `Pos2D{Float64}`.


Now that we have a position for the particle, we'll want to keep track of other information of the particle as well, such as the diameter. We can do this by creating a `stuct` that encompases all that informaiton.

```julia
struct Particle{VecType}
    position::VecType
    diameter::Float64
end
```
The structure `Particle` will take on the type `VecType` which is based on the position vector and has fields of `position` and `diameter`. Later, we can include more information if we want each particle do have different masses or spring contants, etc.

Again, to make things easier, let's make a function that will generate particles with random positions and random diamters given a range for each.

```julia
function ParticleRandom(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end
```

This has a similar flow as `PosRandom` function, except it has another argument for the range of the diameter. Let's make a random particle using this function.

```julia
julia> ExampleParticle = ParticleRandom(Pos2D{Float64}, [0,1], [0,1])
Particle{Pos2D{Float64}}([0.7197480175465489, 0.36301736551172714], 0.0846518956414718)
```

We can access the position of the `ExampleParticle` with

```julia
julia> ExampleParticle.position
2-element Pos2D{Float64} with indices SOneTo(2):
 0.7197480175465489
 0.36301736551172714
 ```
 and the diameter,

 ```julia
 julia> ExampleParticle.diameter
0.0846518956414718
 ```
## Forces

Now that we have a computationally efficient way to define and create particles, we can move on to defining the forces between these points. Let's start with a Hookean spring force law where $\delta$ is the overlap. For particles that have no overlap ($\delta = 0$), there is no force. If the distance between particles $i$ and $j$ is given by:

$$
r_{i,j} = \sqrt{(x_j - x_i)^2 + (y_j - y_i)^2}
$$

where $x_i, y_i$ are the positions of particle $i$, and $x_j, y_j$ are the positions of particle $j$, then the potential energy function can be defined as:

$$
U(r) = \frac{k}{2}\left[\left(\frac{d_i}{2} + \frac{d_j}{2}\right) - r_{i,j} \right]^2 \quad \text{for} \quad r_{i,j} \leq \left(\frac{d_i}{2} + \frac{d_j}{2}\right)
$$

and

$$
U(r) = 0.0 \quad \text{for} \quad r_{i,j} > \left(\frac{d_i}{2} + \frac{d_j}{2}\right)
$$

where $d_i$ and $d_j$ are the diameter of particles $i$ and $j$ respectively. Computationally, this is


```julia
function EnergyHooke(p_i::Particle{VecType}, p_j::Particle{VecType}) where VecType
    r_vector = p_j.position - p_i.position
    r = norm(r_vector)
    cutoff = .5 .* (p_i.diameter + p_i.diameter)
    if r > cutoff
        PotEnergy  = zero(T)
    else
        PotEnergy = .5 * k * (cutoff - r)^2 * (r_vector / r)
    end
    return PotEnergy
end
```

using the `norm` function from `LinearAlgebra` helps us handle the 3D positions without having to change our function. Here `p_i` and `p_j` are the structures `Particle{VecType}` for particles $i$ and $j$ respectively. The output of `EnergyHooke` is a vector that gives you the energy of each dimensional component of `VecType`.

The force exerted on the particles due to this potential is given by:

$$
F_i(r) = k\left[\left(\frac{d_i}{2} + \frac{d_j}{2}\right) - r_{i,j} \right]\quad \text{for} \quad r_{i,j} \leq \left(\frac{d_i}{2} + \frac{d_j}{2}\right)
$$

and $0$ otherwise. For which the corresponding forces are:

$$
\vec{F}_j = -\vec{F}_i.
$$

Computationally,

```julia
function ForceHooke(p_i::Particle{VecType}, p_j::Particle{VecType}) where VecType
    r_vector = p_j.position - p_i.position
    r = norm(r_vector)
    cutoff = 0.5 * (p_i.diameter + p_j.diameter)
    if r > cutoff
        force_i = zero(VecType) 
    else
        force_i = k * (cutoff - r) * (r_vector / r)
    end
    return force_i
end
```
Which gives you the force on particle $$i$$ due to it's interaction with particle $j$.　`ForceHooke` has the same type logic as `EnergyHooke` and outputs forces in the dimensions of `VecType`.

## Calculating Forces for All Particles in Simulation

Let's generate a number of particles and create a function that will calculate the forces between them all. 

### Creating Multiple Random Particles.

For the example, let's create 100 random particles and store all their information into a vector called `particle_list`

```julia
particle_list = [ ParticleRandom(Pos2D{Float64}, [0,50], [1,1.5]) for i in 1:100]
```

This creates a 100-element vector where element of the vector is the `Particle` structure that contains the position and the diameter of the particle. If we want the position of the particle in the 2nd element of the vector,

```julia
julia> particle_list[2].position
2-element Pos2D{Float64} with indices SOneTo(2):
 24.696201789919336
 41.01577656742174
```
and do the same same for diameter.

## Designing the Total Force Function

Now we have our list, let's create a function that will go through each non-repeating pair of particles in `particle_list` and calculate the forces in each direction. We know the inputs of the function should be `particle_list` and some function `ForceLaw` that we choose as the force interaction between the particles. We should also have an output list of forces `force_list` which should have the same dimensions as `VecType`. We can do that with 

```julia
force_list = similar(map( p -> p.position, particle_list))
```
Here, `map( element -> element.position, particle_list)` transforms `particle_list` by applying the anomalous function `element -> element.position` elementwise. `similar` creates an uninitialized version of it's argument. `force_list` will be have a type `Vector{VecType}` with the same length as `particle_list` and each element contains a vector of type `VecType`. For this example, it'll be a vector of position vectors.

Ok, now that we have our inputs and our output, let's create a function that goes initializes `force_list`, and calculates the $x$ and $y$ force components for each particle pair, and finally adds those components to `force_list`

```julia
function forces!(force_list::Vector{VecType}, particle_list::Vector{Particle{VecType}}, ForceLaw::F) where {VecType, F}
    fill!(force_list, zero(VecType)) 
    n = length(particle_list)
    for i in 1:n-1
        for j in i+1:n
            force_i = ForceLaw(particle_list[i], particle_list[j])
            force_list[i] += force_i
            force_list[j] -= force_i
        end
    end
    return force_list
end

```
To call this function,

```julia
forces!(force_list, particle_list, ForceHooke)
```
Later, we can create new force functions and simply use those as our `ForceLaw` argument. Just make sure it outputs the same type as `VecType`!

## Simulation Algorithm

Now it's time to get to the actual simulation. We use a Velocity Verlet integration algorithm where positions of particles are updated. Here's the flow:

1. Given the inital positions and velocities of each particle, calculate the intital forces $$F_0$$ on them.

2. Using those initial forces, calculate the initial accelerations $$a_0$$ on each particle using $$F_0 = ma_0$$.

3. Calculate what the position $$p$$ of each particle will be after $$dt$$ time has passed *assuming no interactions happen in* $$dt$$ using

$$
p = 
$$


## Plotting
I'm just using the simple plotting package `using Plots` to see our trajectories.
```julia

function plot_trajectory(trajectory)
    # Use the first set of positions to determine the limits
    initial_positions, initial_diameters = trajectory[1]
    xlims = (minimum([pos.x for pos in initial_positions]), maximum([pos.x for pos in initial_positions]))
    ylims = (minimum([pos.y for pos in initial_positions]), maximum([pos.y for pos in initial_positions]))

    # This scale factor aligns the ms with the data units.
    marker_size_scale = 22  # Adjust this if the sizes don't match visually as expected.

    @gif for i in 1:length(trajectory)
        positions, diameters = trajectory[i]
        x = [pos.x for pos in positions]
        y = [pos.y for pos in positions]
        
        # ms is the marker size in points, we calculate it based on data units
        scatter(x, y, ms=diameters .* marker_size_scale, legend=false, xlims=xlims, ylims=ylims, markershape=:circle)
        annotate!([(xlims[1], ylims[1], text("step: $i", :bottom))])
    end every 1
end
```

## Adding Periodic Boundaries

As you can see, our particles will all eventually find their own free paths and just continue to go in those directions. Let's set up a boundary condition that allows for more collisions.

We could do hard boundaries, but first let's do periodic boundaries where they go out on one side and come back in on the other side.