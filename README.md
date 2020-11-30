#SimpleTectonics

This is my second / third attempt at good plate tectonics.

I have learned a lot from my last attempt using large voronoi based plates.

The simple approach will not work, needs more complex dynamics.

## Core Concepts in Plate Tectonics

- Lithosphere (Upper Crust Layer) is broken into tectonic plates
- The surface on which the crust moves (asthenosphere moves elastically like viscous flow), deforms plastically
- Key concept: Brittle / strong lithosphere on a weak / hot / fluid asthenosphere
- Thickness of the lithosphere depends on the temperature of the asthenosphere beneath, because it will cool and crystallize accordingly
- Oceanic lithosphere thickens as it ages, and is more dense, due to cooling
- Oceaninc lithosphere moves under land because of higher density, not just height
- Oceaninc lithosphere is younger than land lithosphere because it subduces
- Continental lithosphere is the thickest and most stable, less dense

## Ideas
- GPU Accelerated Voronoise could be used for clustering
- Multiple centroids could constitute a plate and allow for upwards / downwards shifting
- Centroids have to be associated with a thickness and a density
- Denser plates will subduce, less dense plates will move over.

I need a method to model the density and thickness of plates in some kind of equilibrium
system, so when they are forced up and down by buoyancy somehow this is all coupled.

Plates move as a whole based on a sampling of forces of the individual guys.

When two plate segments collide, the one with lower density is forced downwards.

Somehow this process creates sediment which then cascades.



- The movement of the centroids up and down could lead to the raising / lowering of sediment which then cascades very roughly to give the underlying rock formations
- Centroids have to be created and destroyed in an appropriate fashion

- At the boundaries plates go down while they stay in the middle, so having a representation where a single plate can have variable height appears to make sense.

- Individual sections of a plate are therefore also associated with a height or thickness
- Thickness can decrease or increase depending on

## Concept

Plate Tectonics Idea:

Use GPU Accelerated Voronoise to simulate the movement of centroids

The individual centroids belong to the same plate but allow for more granular simluation control.

They are swimming on a surface and move around.

Multiple centroids belong to the same "clusters" basically, which means same soil.

Points move around randomly and squish and squash.

The question is: what happens when they collide?
Can centroids have a height value which means they have stronger influence?

/*
    Plates track their Centroids

    //Plate Movement
    Based on a heat map, we can compute centroid convection
    Convection is proportional to the gradient obvisously.
    The heat map is initialized with noise.
    The centroid convections are used to compute the force and torque of the plate
    we can compute the plate center, direction and rotation
    we can then move the centroids

    //Centroid Height
    Centroids have a thickness and a density.
    Age of a segment dictates its thickness, which grows at a certain rate,
    while the heat dictates an equilibrium density which is the density the plate grows at.
    The centroid thereby stores a certain amount of mass over its area.

    Additionally, the centroid's 3D position is determined by its density and volume.
    It is basically floating, and there is a certain amount above and below the surface.
    Density is initially 0.5, so that it floats half above.

    //Plate Collision
    When plates collide, the less dense plate is pushed below and the more dense plate is pushed above.
    This creates sediment and heat. How much sediment and how much heat? Basically we take the combined
    height and density, compute the overshoot on the top and bottom to compute sedimentation and heating rate.

    When plates separate, new plate is created which has the density given by the heat below.
    That is all.

    //Heat Generation
    Heat can initially be seeded with random noise and we can see how plates move.
    Later heat can be generated randomly or can be affected by the centroid growth.
    When centroids collide, they generate heat by destruction of mass.

    The area of the less dense plate is destroyed proportionally to its submerged density

    //Sediment Transport and Heat Generation
    Computing where and how plates move relative to each other,
    we can compute where to generate heat.

    NOTE:
    As a plate with a high density, i.e. sitting at a ratio other than 1:1 in the ground,
    then as the thickness decreases due to collision, the plate will sink further down faster.
    The thickness of the plate needs to decrease as it subduces.

*/
