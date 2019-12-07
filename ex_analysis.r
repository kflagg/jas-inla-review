#######################
## PACKAGES AND DATA ##
#######################

library(spatstat) # For bei dataset
library(INLA) # For model fitting.
#library(maptools) # Provides as.SpatialPolygons method for owin.

# Plot the observed point pattern and its window.
# TODO: add an informative plot title.
pdf('figures/bei.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei)
dev.off()


#######################
## MESH CONSTRUCTION ##
#######################

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
pdf('figures/beimesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_mesh, asp = 1, main = '')
points(bei, pch = 19, cex = 0.25, col = 'red')
title('Mesh Over Bei Data')
dev.off()


######################################
## SET UP SIMPSON METHOD PSEUDODATA ##
######################################

# Observed event locations.
bei_pts <- cbind(bei$x, bei$y)

# Get the numbers of mesh nodes and real events.
# The sum of these will be the number of pseudodata points.
bei_mesh_size <- bei_mesh$n
bei_n_events <- nrow(bei_pts)

# Compute the barycentric coordinates of the observed events
# (i.e. project into the space spanned by the basis of mesh nodes).
bei_bary <- inla.mesh.project(bei_mesh, bei_pts)$A

# Compute the barycentric coordinates of the nodes. Because the
# node coordinatess are the basis vectors, this is an identity matrix.
bei_int_matrix <- sparseMatrix(
  i = seq_len(bei_mesh_size),
  j = seq_len(bei_mesh_size),
  x = rep(1, bei_mesh_size)
)

# Bind the node and event coordinates into a single matrix of pseudodata locations
# in barycentric coordinates.
bei_pseudopoints <- rbind(bei_int_matrix, bei_bary)

# Get the numerical integration weights for the SPDE approach.
# Because the whole region was fully observed (constant sampling effort),
# these are left as is. If parts of the region were unobserved, these would be
# multiplied by the observed proportion of the area represented by the node.
bei_int_weights <- diag(inla.mesh.fem(bei_mesh)$c0)

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
bei_pseudodata_exp <- c(bei_int_weights, rep(0, bei_n_events))

# Create the psuedodata. This is a vector giving the count of events at each
# pseudopoint location, that is 0 at each mesh node and 1 at each event.
bei_pseudodata <- c(rep(0, bei_mesh_size), rep(1, bei_n_events))


########################################
## FIT THE MODEL AND PLOT THE RESULTS ##
########################################

# Define the model formula.
# inla() requires a vector of 1s in the data argument for the intercept.
# The f() term specifies the GP with observations indexed by a variable called
# idx. The indices correspond to the indixes of the mesh nodes.
bei_formula <- y ~ -1 + intercept + f(idx, model = bei_spde)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
bei_inla_data <- list(
  y = bei_pseudodata, # The whole pseudodata vector.
  idx = seq_len(bei_mesh_size), # Indices of the nodes.
  intercept = rep(1, bei_mesh_size) # Intercept column.
)

# Fit the model as a Poisson GLM with exposures specified.
bei_result <- inla(
  formula = bei_formula,
  data = bei_inla_data,
  family = 'poisson',
  control.predictor = list(A = bei_pseudopoints),
  E = bei_pseudodata_exp
)

# Plot the posterior mean of the latent surface.
pdf('figures/beimean.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Latent GP')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = 'white')
dev.off()

# Plot the posterior standard deviation of the latent surface.
pdf('figures/beisd.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$sd)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Latent GP')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = 'white')
dev.off()
