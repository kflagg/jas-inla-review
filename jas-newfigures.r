library(geoR)
library(INLA)

# Create an example mesh over the unit square.
ex_mesh <- inla.mesh.2d(
  loc.domain = cbind(c(0, 0, 1, 1), c(0, 1, 0, 1)),
  max.edge = 0.1,
  offset = 0
)

# Prepare a lattice on which to realize a GRF.
sim_x <- seq(0, 1, 0.001)
sim_y <- seq(0, 1, 0.001)
sim_pts <- as.matrix(expand.grid(x = sim_x, y = sim_y))
sim_positions <- rbind(
  ex_mesh$loc[,1:2],
  sim_pts
)

# Simulate the GRF.
set.seed(3525)
sim_grf <- grf(nrow(sim_positions), sim_positions, cov.pars = c(1, 0.6))

# Project the GRF onto the mesh and interpolate.
ex_spde <- inla.spde2.matern(mesh = ex_mesh)
ex_A <- inla.mesh.project(ex_mesh, sim_positions)$A
ex_interp <- ex_A %*% head(sim_grf$data, n = ex_mesh$n)

# Create vectors of mesh edge start and endpoints.
vert_idx <- cbind(ex_mesh$graph$tv[,c(1:3, 1), drop = FALSE], NA)
edge_x <- ex_mesh$loc[t(vert_idx), 1]
edge_y <- ex_mesh$loc[t(vert_idx), 2]
edge_z <- sim_grf$data[t(vert_idx)]

# Create a plot of the mesh.
plot(ex_mesh, main = '', asp = 1)
points(ex_mesh$loc[,1:2], pch = 16)
points(0.5, 0.5, pch = 19, cex = 1.5)
title('Example Mesh')
# TODO: Highlight node 9 (at (0.5, 0.5)) and shade its dual cell.
