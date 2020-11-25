# SimpleTectonics

This is my second / third attempt at good plate tectonics.

I have learned a lot from my last attempt using large voronoi based plates.

The simple approach will not work, needs more complex dynamics.

## Ideas
- GPU Accelerated Voronoise could be used for clustering
- Multiple centroids could constitute a plate and allow for upwards / downwards shifting
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
