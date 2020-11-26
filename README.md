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
