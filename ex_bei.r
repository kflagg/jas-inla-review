#######################
## PACKAGES AND DATA ##
#######################

library(spatstat) # For bei dataset.
library(INLA) # For model fitting.

# Plot the observed point pattern and its window.
pdf('figures/bei.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei, main = expression(paste(italic('Beilschmiedia pendula Lauraceae'), ' Locations')))
dev.off()

# Plot the elevation surface.
pdf('figures/beielev.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$elev, main = 'Elevation')
points(bei, pch = '.', col = 'white')
dev.off()

# Plot the gradient surface.
pdf('figures/beigrad.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$grad, main = 'Gradient')
points(bei, pch = '.', col = 'white')
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
  refine = list(max.edge = 25)
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


################################
## PROJECT COVARIATES TO MESH ##
################################

# The covariates are measured on a lattice which does not match up to the mesh
# we'll use for the intensity function. There isn't a well-developed way to
# handle the change-of-support within the model fit in INLA (e.g. by adding a
# level to the model hierarchy to predict the covariates at the nodes).
# However, these covariate surfaces are smooth enough that little will be lost
# by using piecewise linear interpolation to find values at the nodes.
# (Using an overly coarse mesh would be a bigger problem than interpolation).

# Convert the covariate images to data frames.
bei_elev <- as.data.frame(bei.extra$elev)
bei_grad <- as.data.frame(bei.extra$grad)

# Create a temporary mesh with nodes at the observation locations.
bei_covariate_locs <- inla.mesh.create(
  boundary = bei_boundary,
  loc = cbind(bei_elev$x, bei_elev$y)
)

# Compute the barycentric coordinates of the nodes with respect to the
# temporary mesh. This will be the projection matrix that does the
# change-of-support.
bei_change_of_support <- inla.mesh.project(
  bei_covariate_locs,
  bei_mesh$loc[,1:2]
)$A

# Project elevation and gradient onto the nodes.
bei_mesh_elev <- as.vector(bei_change_of_support %*% bei_elev$value)
bei_mesh_grad <- as.vector(bei_change_of_support %*% bei_grad$value)

# Plot the piecewise linear approximation of the elevation surface.
pdf('figures/beielevmesh.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_elev)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        main = 'Piecewise Linear Approximation of Elevation')
plot(Window(bei), border = '#80808080', add = TRUE)
points(bei, pch = '.', col = 'white')
dev.off()

# Plot the piecewise linear approximation of the gradient surface.
pdf('figures/beigradmesh.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_grad)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        main = 'Piecewise Linear Approximation of Gradient')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#80808080')
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

# Create the psuedodata. This is a vector giving the count of events at each
# pseudopoint location, that is 0 at each mesh node and 1 at each event.
bei_pseudodata <- c(rep(0, bei_mesh_size), rep(1, bei_n_events))

# Get the numerical integration weights for the SPDE approach.
# Because the whole region was fully observed (constant sampling effort),
# these are left as is. If parts of the region were unobserved, these would be
# multiplied by the observed proportion of the area represented by the node.
bei_int_weights <- diag(inla.mesh.fem(bei_mesh)$c0)

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
bei_pseudodata_exp <- c(bei_int_weights, rep(0, bei_n_events))

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

# Bind the node and event coordinates into a single matrix of pseudodata
# locations in barycentric coordinates.
bei_pseudopoints <- rbind(bei_int_matrix, bei_bary)


########################################
## FIT THE MODEL AND PLOT THE RESULTS ##
########################################

# Define the model formula.
# inla() requires a vector of 1s in the data argument for the intercept.
# The f() term specifies the GP with observations indexed by a variable called
# idx. The indices correspond to the indixes of the mesh nodes.
bei_formula <- y ~ -1 + intercept + elev + grad + f(idx, model = bei_spde)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
bei_inla_data <- list(
  y = bei_pseudodata, # The whole pseudodata vector.
  elev = bei_mesh_elev, # Elevation at nodes.
  grad = bei_mesh_grad, # Gradient at nodes.
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
     riblab = expression(E(bold(e)(u)*'|'*bold(x))),
     main = 'Posterior Predicted Mean of Latent GP')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
pdf('figures/beisd.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$sd)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))),
     main = 'Posterior Prediction SD of Latent GP')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/beibetas.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad)
        ),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(beta[0]*'|'*bold(x)) +
                         E(beta[1]*'|'*bold(x)) * z[1](u) +
                         E(beta[2]*'|'*bold(x)) * z[2](u)),
     main = 'Posterior Mean of Fixed Effects')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#ffffff80')
dev.off()


# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/beiintensity.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(exp(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad) +
          inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)
        )),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
     main = 'Posterior Intensity Function')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#80808080')
dev.off()


# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/beipost.pdf', width = 18, height = 6)
par(mfrow = c(1, 3), bty = 'n')
plot(inla.smarginal(bei_result$marginals.hyperpar$Theta1), type = 'l',
     yaxt = 'n', xlab = expression(tau),
     main = expression(paste(bold('Posterior Distribution of '), tau)))
plot(inla.smarginal(bei_result$marginals.hyperpar$Theta2), type = 'l',
     yaxt = 'n', xlab = expression(kappa),
     main = expression(paste(bold('Posterior Distribution of '), kappa)))
plot(inla.smarginal(bei_result$marginals.fixed$intercept), type = 'l',
     yaxt = 'n', xlab = expression(beta[0]),
     main = 'Posterior Distribution of the Intercept')
dev.off()

# Plot posterior marginals of the coefficients.
pdf('figures/beipostcoefs.pdf', width = 12, height = 6)
par(mfrow = c(1, 2), bty = 'n')
plot(inla.smarginal(bei_result$marginals.fixed$elev), type = 'l',
     yaxt = 'n', xlab = expression(beta[1]),
     main = 'Posterior Distribution of Elevation Coefficient')
plot(inla.smarginal(bei_result$marginals.fixed$grad), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Gradient Coefficient')
dev.off()


####################
## MODEL CHECKING ##
####################

# Create grid counts for residual calculation.

# Number of grid cells.
NGRID_X <- 40
NGRID_Y <- 20

# Have spatstat find centers for the grid cells.
centers <- gridcenters(
  owin(range(bei_boundary$loc[,1]), range(bei_boundary$loc[,2])),
  NGRID_X, NGRID_Y)

# Compute the grid cell size.
dx <- sum(unique(centers$x)[1:2] * c(-1, 1)) / 2
dy <- sum(unique(centers$y)[1:2] * c(-1, 1)) / 2

# Initialize a data frame to store the counts.
bei_resid_df <- data.frame(x = centers$x, y = centers$y,
                           count = NA_integer_, area = NA_real_)

# Loop through the cells, finding the event count and cell area.
for(r in seq_len(nrow(bei_resid_df))){
  bei_resid_df$count[r] <- sum(bei$x >= bei_resid_df$x[r] - dx &
                               bei$x < bei_resid_df$x[r] + dx &
                               bei$y >= bei_resid_df$y[r] - dy &
                               bei$y < bei_resid_df$y[r] + dy)
  bei_resid_df$area[r] <- area(Window(bei)[owin(c(bei_resid_df$x[r] - dx,
                                                  bei_resid_df$x[r] + dx),
                                                c(bei_resid_df$y[r] - dy,
                                                  bei_resid_df$y[r] + dy))])
}

# Set up a projection from the SPDE representation to the residual grid.
bei_resid_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(as.data.frame(centers)))

# Project the intensity surface to the residual grid.
bei_resid_df$intensity <- as.vector(exp(
  bei_result$summary.fixed['intercept', 'mean'] +
    bei_result$summary.fixed['elev', 'mean'] *
  inla.mesh.project(bei_resid_proj, bei_mesh_elev) +
    bei_result$summary.fixed['grad', 'mean'] *
  inla.mesh.project(bei_resid_proj, bei_mesh_grad) +
  inla.mesh.project(bei_resid_proj, bei_result$summary.random$idx$mean)
))

# Calculate the (approximate) expected counts and the Pearson residuals.
bei_resid_df$expected <- bei_resid_df$intensity * bei_resid_df$area
bei_resid_df$pearson <- (bei_resid_df$count - bei_resid_df$expected) /
  sqrt(bei_resid_df$expected)

# Gridded Pearson residual plot.
pdf('figures/beiresiduals.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(matrix(bei_resid_df$pearson, nrow = length(unique(bei_resid_df$x)))),
        unique(bei_resid_df$x), unique(bei_resid_df$y),
        unitname = c('meter', 'meters')),
     main = 'Gridded Pearson Residuals')
plot(Window(bei), border = 'white', add = TRUE)
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Set up a projection from the SPDE representation to the event locations.
bei_event_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(as.data.frame(bei)))

# Create a copy of bei and mark with 1/sqrt(lambda).
bei_marked <- bei
marks(bei_marked) <- as.vector(1/sqrt(exp(
  bei_result$summary.fixed['intercept', 'mean'] +
    bei_result$summary.fixed['elev', 'mean'] *
  inla.mesh.project(bei_event_proj, bei_mesh_elev) +
    bei_result$summary.fixed['grad', 'mean'] *
  inla.mesh.project(bei_event_proj, bei_mesh_grad) +
  inla.mesh.project(bei_event_proj, bei_result$summary.random$idx$mean)
)))

# Mark plot.
pdf('figures/beimarkplot.pdf', width = 12, height = 6)
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(exp(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad) +
          inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)
        ))),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        main = 'Mark Plot')
plot(Window(bei), border = 'white', add = TRUE)
plot(bei_marked, col = '#ffffff80', add = TRUE)
dev.off()
