# SimpleTectonics

C++-17 implementation of a clustered convection based plate tectonics simulation for procedural terrain generation

Implemented using [TinyEngine](https://github.com/weigert/TinyEngine)

[Link to Blog Post](https://nickmcd.me/2020/12/03/clustered-convection-for-simulating-plate-tectonics/)

![Tectonics Banner with Plates and Terrain](https://github.com/weigert/SimpleTectonics/blob/master/screenshots/banner1.png)

## Compilation

Use the makefile to compile the program

## Dependencies

    - TinyEngine
    - libnoise

## Usage

    ./tectonics [SEED]

If no seed is specified, it will pick a random one.

## Controls

    Zoom and Rotate: Scroll X / Scroll Y
    Toggle Pause: P (Paused by Default)
    Tilt Camera Vertically: Up / Down Key
    View Control Panel: ESC
    Toggle Map View: M
    Move Camera Anchor: WASD / Space / C

## Screenshots

Here are two example generations:

![Plates Example 0](https://github.com/weigert/SimpleTectonics/blob/master/screenshots/plates0.png)
![Terrain Example 0](https://github.com/weigert/SimpleTectonics/blob/master/screenshots/terrain0.png)
![Plates Example 1](https://github.com/weigert/SimpleTectonics/blob/master/screenshots/plates1.png)
![Terrain Example 1](https://github.com/weigert/SimpleTectonics/blob/master/screenshots/terrain1.png)

## Reading

The main file contains all the relevant code for initializing the system and rendering the data to screen.

A number of header files execute the code logic:

    source/poisson.h - Naive Poisson Disc Sampler for generating the centroids
    source/scene.h - Rendering Data (Camera / Lighting / Colors) and Heightmap Mesher (i.e. just visualization stuff)
    source/cluster.h - Clustered Convection system in generalized form
    soruce/tectonics.h - System that utilizes cluster.h to simulate the plate tectonics physical process

Additionally, there are a number of important shaders:

    source/shader/cascading.fs - Simple Shader-Based Sand-Pile Cascading Implementation
    source/shader/diffusion.fs - Simple Shader-Based Diffusion Algorithm (for Heatmap)
    source/shader/subduction.fs - Shader for turning colliding centroids into heatmap

    source/shader/voronoi.* - GPU Accelerated Voronoise Texture - https://github.com/weigert/TinyEngine/tree/master/examples/11_Voronoi

    source/shader/depth.* - Depth Map Generation for Shadow Mapping
    source/shader/flat.* - Shader for Rendering Billboards to FBO


## License

MIT License
