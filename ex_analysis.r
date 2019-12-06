## PACKAGES AND DATA

library(spatstat) # For bei dataset
library(INLA) # For model fitting.
#library(maptools) # Provides as.SpatialPolygons method for owin.

# Plot the observed point pattern and its window.
plot(bei)


## MESH CONSTRUCTION

# Get the boundary polygon for the observation window
# and define mesh edge segments for INLA.
bei_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(Window(bei))))

# Create a Delaunay triangulation with maximum edge length of 50 meters to use
# as the mesh.
bei_mesh <- inla.mesh.create(
  boundary = bei_boundary,
  refine = list(max.edge = 50)
)

# Define a SPDE representation of a GP with Matern covariance.
bei_spde <- inla.spde2.matern(mesh = bei_mesh)

# Set up a projection from the SPDE representation to a 400x200 grid.
bei_proj <- inla.mesh.projector(bei_mesh, dims = c(400, 200))


# Plot the mesh and point pattern.
plot(bei_mesh, asp = 1, main = '')
points(bei, pch = 19, cex = 0.25, col = 'red')
title('Mesh Over Bei Data')
