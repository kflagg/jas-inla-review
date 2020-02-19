#######################
## PACKAGES AND DATA ##
#######################

library(spatstat)
library(INLA)
library(DSpat) # For dataset.

data(DSpat.obs) # Incompletely observed (distance sampling) dataset.
data(DSpat.lines) # Transect lines used for distance sampling.
data(DSpat.covariates) # Lattice of predictor variables.

# Create spatstat objects for plotting.
site_win <- owin(c(0, 100), c(0, 100))
obs_psp <- psp(DSpat.lines$x0, DSpat.lines$y0, DSpat.lines$x1, DSpat.lines$y1, site_win)
obs_win <- dilation(obs_psp, 0.2)[site_win] # Transects have a half-width of 0.2.
dspat_ppp <- ppp(DSpat.obs$x, DSpat.obs$y, window = obs_win)
dspat_covs_ppp <- as.ppp(DSpat.covariates, site_win)
dspat_river_im <- im(
  matrix(DSpat.covariates$river, ncol = length(unique(DSpat.covariates$x))),
  xcol = unique(DSpat.covariates$x), yrow = unique(DSpat.covariates$y)
)
dspat_habitat_im <- im(
  matrix(DSpat.covariates$habitat, ncol = length(unique(DSpat.covariates$x))),
  xcol = unique(DSpat.covariates$x), yrow = unique(DSpat.covariates$y)
)
dspat_dist_im <- im(t(matrix(
  apply(sapply(unique(DSpat.obs$label), function(l){
      return(dist2line(
        as.ppp(expand.grid(x = 0:100, y = 0:100), site_win),
        DSpat.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance
      )
    }), 1, min), nrow = 101)),
  xrange = site_win$x, yrange = site_win$y)


# Plot the observed point pattern, its window, and the transects.
pdf('figures/dspat.pdf', width = 6, height = 6)
par(mar = c(3, 3, 2, 0))
plot(dspat_ppp, main = 'Event Locations', border = 'grey')
plot(site_win, add = TRUE)
axis(1, at = seq(0, 100, 20))
mtext('Horizontal Coordinate', 1, 2)
axis(2, at = seq(0, 100, 20))
mtext('Vertical Coordinate', 2, 2)
dev.off()

# Plot the river covariate and the transects.
pdf('figures/dspatriver.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(dspat_river_im, main = 'River', border = 'grey',
     riblab = 'River Covariate', ribsep = 0.05)
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the river covariate and the transects.
pdf('figures/dspathabitat.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(dspat_habitat_im, main = 'Habitat', border = 'grey',
     riblab = 'Habitat Covariate', ribsep = 0.05)
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the distance-to-transect surface.
pdf('figures/dspatdist.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(dspat_dist_im, riblab = 'Units', ribsep = 0.05,
     main = 'Distance to Transect')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


#######################
## MESH CONSTRUCTION ##
#######################

#TODO: refine mesh near transects.

# Get the boundary polygon for the study region window
# and define mesh edge segments for INLA.
dspat_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(site_win)))

# Create a Delaunay triangulation with maximum edge length of 2.5 units to use
# as the mesh.
dspat_mesh <- inla.mesh.create(
  boundary = dspat_boundary,
  refine = list(max.edge = 5)
)

# Convert the mesh edges to a psp.
meshloc <- dspat_mesh$loc[,-3]
meshadj <- dspat_mesh$graph$vv
meshadj[lower.tri(meshadj)] <- 0
meshsegidx0 <- do.call(c,
  apply(cbind(seq_len(ncol(meshadj)), apply(meshadj, 2, sum)), 1,
    function(x){return(rep(x[1], x[2]))})
  )
meshsegidx1 <- do.call(c, apply(meshadj == 1, 2, which))
meshseg0 <- meshloc[meshsegidx0,]
meshseg1 <- meshloc[meshsegidx1,]
mesh_psp <- psp(
  x0 = meshseg0[,1],
  y0 = meshseg0[,2],
  x1 = meshseg1[,1],
  y1 = meshseg1[,2],
  window = owin(range(meshloc[,1]), range(meshloc[,2])))
mesh_win <- convexhull(mesh_psp)
Window(mesh_psp) <- mesh_win

# Define a SPDE representation of a GP with Matern covariance.
dspat_spde <- inla.spde2.matern(mesh = dspat_mesh)

# Set up a projection from the SPDE representation to a 400x400 grid.
dspat_proj <- inla.mesh.projector(dspat_mesh, dims = c(400, 400))

# Plot the mesh and point pattern.
pdf('figures/dspatmesh.pdf', width = 6, height = 6)
par(mar = c(0, 0, 2, 0))
plot(dspat_mesh, asp = 1, main = '')
plot(obs_win, add = TRUE, col = '#ff000080', border = NA)
#plot(obs_psp, add = TRUE, col = 'red')
points(dspat_ppp, pch = 20, col = 'red')
title('Mesh Over Event Data')
dev.off()


################################
## PROJECT COVARIATES TO MESH ##
################################

# Create a ppp version of the mesh nodes to help with calculations.
dspat_mesh_ppp <- as.ppp(dspat_mesh$loc[,1:2], site_win)

# Create a temporary mesh with nodes at the observation locations.
dspat_covariate_locs <- inla.mesh.create(
  boundary = dspat_boundary,
  loc = cbind(DSpat.covariates$x, DSpat.covariates$y)
)

dspat_change_of_support <- inla.mesh.project(
  dspat_covariate_locs,
  dspat_mesh$loc[,1:2]
)$A[,-(1:4)]

# Project the river and habitat covariates onto the nodes.
dspat_mesh_river <- as.vector(dspat_change_of_support %*% DSpat.covariates$river)
dspat_mesh_habitat <- as.vector(dspat_change_of_support %*% DSpat.covariates$habitat)

# Get the distance from each node to the nearest transect.
dspat_mesh_dist <- apply(sapply(unique(DSpat.obs$label), function(l){
  return(dist2line(dspat_mesh_ppp, DSpat.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance)
}), 1, min)

# Plot the piecewise linear approximation of the river surface.
pdf('figures/dspatrivermesh.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(dspat_proj, dspat_mesh_river)),
        xrange = site_win$x,
        yrange = site_win$y),
        riblab = 'River Covariate', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of River')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the habitat surface.
pdf('figures/dspathabitatmesh.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(dspat_proj, dspat_mesh_habitat)),
        xrange = site_win$x,
        yrange = site_win$y),
        riblab = 'Habitat Covariate', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Habitat')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the distance surface.
pdf('figures/dspatdistmesh.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(dspat_proj, dspat_mesh_dist)),
        xrange = site_win$x,
        yrange = site_win$y),
        riblab = 'Units', ribsep = 0.05,# zlim = c(0, 75),
        main = 'Piecewise Linear Approximation of Distance to Transect')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
#plot(obs_psp, add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


######################################
## SET UP SIMPSON METHOD PSEUDODATA ##
######################################

# Observed event locations.
dspat_pts <- cbind(DSpat.obs$x, DSpat.obs$y)

# Get the numbers of mesh nodes and real events.
# The sum of these will be the number of pseudodata points.
dspat_mesh_size <- dspat_mesh$n
dspat_n_events <- nrow(dspat_pts)

# Create the psuedodata. This is a vector giving the count of events at each
# pseudopoint location, that is 0 at each mesh node and 1 at each event.
dspat_pseudodata <- c(rep(0, dspat_mesh_size), rep(1, dspat_n_events))

# Take the intersection of obs_psp with each triangle.
mesh_tris <- apply(dspat_mesh$graph$tv, 1, function(x){
  return(owin(poly = meshloc[x,]))
})
tri_areas <- sapply(mesh_tris, area)
lines_subsegs <- lapply(mesh_tris, function(x){return(obs_psp[x])})

# Calculate the area represented by each node.
mesh_area <- rep(0, dspat_mesh$n)
for(i in seq_along(tri_areas)){
  mesh_area[dspat_mesh$graph$tv[i,]] <- mesh_area[dspat_mesh$graph$tv[i,]] +
    tri_areas[i] / 3
}

# Track which triangle each segment came from.
seg_tri_idx <- unlist(lapply(seq_along(lines_subsegs), function(i){return(rep(i, lines_subsegs[[i]]$n))}))

# Get the midpoints
lines_midpoints <- as.matrix(do.call(rbind, lapply(lapply(lines_subsegs, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
lines_bary <- inla.mesh.projector(dspat_mesh, loc = lines_midpoints)$proj$A

# Get the lengths within each triangle.
lines_lengths <- unlist(sapply(lines_subsegs, lengths.psp))

# Calculate the numerical integration weights for the SPDE approach.
# For each row of the barycentric coordinates matrix (which is number of
# segments by number of nodes), divide each entry by the portion of the
# triangle's area represented by that node, which is 1/3rd of the area
# of the triangle the segment is in.
lines_bary_prop <- as(t(sapply(seq_along(seg_tri_idx),
    function(i){return(lines_bary[i,] / tri_areas[seg_tri_idx[i]] * 3)})),
  'sparseMatrix')

# Multiply by the observed area represented in each triangle.
mesh_weights <- as.vector(0.4 * lines_lengths %*% lines_bary_prop)

# Get the unadjusted integration weights
dspat_int_weights <- diag(inla.mesh.fem(dspat_mesh)$c0)

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
dspat_pseudodata_exp <- c(dspat_int_weights * mesh_weights, rep(0, dspat_n_events))

# Compute the barycentric coordinates of the observed events
# (i.e. project into the space spanned by the basis of mesh nodes).
dspat_bary <- inla.mesh.project(dspat_mesh, dspat_pts)$A

# Compute the barycentric coordinates of the nodes. Because the
# node coordinatess are the basis vectors, this is an identity matrix.
dspat_int_matrix <- sparseMatrix(
  i = seq_len(dspat_mesh_size),
  j = seq_len(dspat_mesh_size),
  x = rep(1, dspat_mesh_size)
)

# Bind the node and event coordinates into a single matrix of pseudodata
# locations in barycentric coordinates.
dspat_pseudopoints <- rbind(dspat_int_matrix, dspat_bary)


########################################
## FIT THE MODEL AND PLOT THE RESULTS ##
########################################

# Define the model formula.
# inla() requires a vector of 1s in the data argument for the intercept.
# The f() term specifies the GP with observations indexed by a variable called
# idx. The indices correspond to the indixes of the mesh nodes.
dspat_formula <- y ~ -1 + intercept + river + habitat + I(habitat^2) + dist + f(idx, model = dspat_spde)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
dspat_inla_data <- list(
  y = dspat_pseudodata, # The whole pseudodata vector.
  river = dspat_mesh_river, # River covariate at nodes.
  habitat = dspat_mesh_habitat, # Habitat covariate at nodes.
  dist = dspat_mesh_dist, # Distance at nodes.
  idx = seq_len(dspat_mesh_size), # Indices of the nodes.
  intercept = rep(1, dspat_mesh_size) # Intercept column.
)

# Fit the model as a Poisson GLM with exposures specified.
dspat_result <- inla(
  formula = dspat_formula,
  data = dspat_inla_data,
  family = 'poisson',
  control.predictor = list(A = dspat_pseudopoints),
  E = dspat_pseudodata_exp
)

# Plot the posterior mean of the latent surface.
pdf('figures/dspatmean.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$mean)),
        xrange = site_win$x,
        yrange = site_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Predicted Mean of Latent GP')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
pdf('figures/dspatsd.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$sd)),
        xrange = site_win$x,
        yrange = site_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Prediction SD of Latent GP')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/dspatbetas.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(
          dspat_result$summary.fixed['intercept', 'mean'] +
          dspat_result$summary.fixed['river', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_river) +
          dspat_result$summary.fixed['habitat', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat) +
          dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat^2)
        ),
        xrange = site_win$x,
        yrange = site_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(beta[0]*'|'*bold(x)) +
                         E(beta[1]*'|'*bold(x)) * z[1](u) +
                         E(beta[2]*'|'*bold(x)) * z[2](u) +
                         E(beta[2]*'|'*bold(x)) * z[2](u)^2), ribsep = 0.05,
     main = 'Posterior Mean of Fixed Effects')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the detection surface.
pdf('figures/dspatdetection.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          dspat_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_dist)
        )),
        xrange = site_win$x,
        yrange = site_win$y),
     zlim = 0:1, riblab = 'Detection Probability', ribsep = 0.05,
     main = 'Posterior Detection Surface')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/dspatintensity.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          dspat_result$summary.fixed['intercept', 'mean'] +
          dspat_result$summary.fixed['river', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_river) +
          dspat_result$summary.fixed['habitat', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat) +
          dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat^2) +
          inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$mean)
        )),
        xrange = site_win$x,
        yrange = site_win$y),
     riblab = 'Events per Square Unit', ribsep = 0.05,
     main = 'Posterior Intensity Function')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface
# for the thinned process.
pdf('figures/dspatthinnedintensity.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          dspat_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_dist) +
          dspat_result$summary.fixed['intercept', 'mean'] +
          dspat_result$summary.fixed['river', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_river) +
          dspat_result$summary.fixed['habitat', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat) +
          dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat^2) +
          inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$mean)
        )),
        xrange = site_win$x,
        yrange = site_win$y),
     riblab = 'Observable Events per Square Unit', ribsep = 0.05,
     main = 'Posterior Thinned Intensity Function')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/dspatpost.pdf', width = 18, height = 6)
par(mfrow = c(1, 3), bty = 'n')
plot(inla.smarginal(dspat_result$marginals.hyperpar$Theta1), type = 'l',
     yaxt = 'n', xlab = expression(tau),
     main = expression(paste(bold('Posterior Distribution of '), tau)))
plot(inla.smarginal(dspat_result$marginals.hyperpar$Theta2), type = 'l',
     yaxt = 'n', xlab = expression(kappa),
     main = expression(paste(bold('Posterior Distribution of '), kappa)))
plot(inla.smarginal(dspat_result$marginals.fixed$dist), type = 'l',
     yaxt = 'n', xlab = expression(beta[4]),
     main = 'Posterior Distribution of Distance Coefficient')
dev.off()

# Plot posterior marginals of the coefficients.
pdf('figures/dspatpostcoefs.pdf', width = 12, height = 12)
par(mfrow = c(2, 2), bty = 'n')
plot(inla.smarginal(dspat_result$marginals.fixed$intercept), type = 'l',
     yaxt = 'n', xlab = expression(beta[0]),
     main = 'Posterior Distribution of the Intercept')
plot(inla.smarginal(dspat_result$marginals.fixed$river), type = 'l',
     yaxt = 'n', xlab = expression(beta[1]),
     main = 'Posterior Distribution of River Coefficient')
plot(inla.smarginal(dspat_result$marginals.fixed$habitat), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Habitat Linear Coefficient')
plot(inla.smarginal(dspat_result$marginals.fixed$`I(habitat^2)`), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Habitat Quadratic Coefficient')
dev.off()

# Plot the detection functions.
pdf('figures/dspatdetectcurves.pdf', width = 7, height = 6)
par(bty = 'n')
curve(exp(dspat_result$summary.fixed['dist', 'mean'] * x),
      from = 0, to = 10, ylim = 0:1, lty = 1,
      ylab = 'Detection Probability', xlab = 'Distance',
      main = 'Posterior Detection Function')
dev.off()


####################
## MODEL CHECKING ##
####################

# Residuls mostly positive. It seems lambda is biased low.

# Create grid counts for residual calculation.

# Number of grid cells.
NGRID_X <- 40
NGRID_Y <- 40

# Have spatstat find centers for the grid cells.
centers <- gridcenters(
  owin(range(dspat_boundary$loc[,1]), range(dspat_boundary$loc[,2])),
  NGRID_X, NGRID_Y)

# Compute the grid cell size.
dx <- sum(unique(centers$x)[1:2] * c(-1, 1)) / 2
dy <- sum(unique(centers$y)[1:2] * c(-1, 1)) / 2

# Initialize a data frame to store the counts.
dspat_resid_df <- data.frame(x = centers$x, y = centers$y,
                             count = NA_integer_, area = NA_real_)

# Loop through the cells, finding the event count and cell area.
for(r in seq_len(nrow(dspat_resid_df))){
  dspat_resid_df$count[r] <- sum(dspat_ppp$x >= dspat_resid_df$x[r] - dx &
                                 dspat_ppp$x < dspat_resid_df$x[r] + dx &
                                 dspat_ppp$y >= dspat_resid_df$y[r] - dy &
                                 dspat_ppp$y < dspat_resid_df$y[r] + dy)
  dspat_resid_df$area[r] <- area(obs_win[owin(c(dspat_resid_df$x[r] - dx,
                                                dspat_resid_df$x[r] + dx),
                                              c(dspat_resid_df$y[r] - dy,
                                                dspat_resid_df$y[r] + dy))])
}

# Set up a projection from the SPDE representation to the residual grid.
dspat_resid_proj <- inla.mesh.projector(dspat_mesh,
  loc = as.matrix(as.data.frame(centers)))

# Project the intensity surface to the residual grid.
dspat_resid_df$intensity <- as.vector(exp(
  dspat_result$summary.fixed['dist', 'mean'] *
    inla.mesh.project(dspat_resid_proj, dspat_mesh_dist) +
  dspat_result$summary.fixed['intercept', 'mean'] +
  dspat_result$summary.fixed['river', 'mean'] *
    inla.mesh.project(dspat_resid_proj, dspat_mesh_river) +
  dspat_result$summary.fixed['habitat', 'mean'] *
    inla.mesh.project(dspat_resid_proj, dspat_mesh_habitat) +
  dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
    inla.mesh.project(dspat_resid_proj, dspat_mesh_habitat^2) +
  inla.mesh.project(dspat_resid_proj, dspat_result$summary.random$idx$mean)
))

# Calculate the (approximate) expected counts and the Pearson residuals.
dspat_resid_df$expected <- dspat_resid_df$intensity * dspat_resid_df$area
dspat_resid_df$pearson <- (dspat_resid_df$count - dspat_resid_df$expected) /
  sqrt(dspat_resid_df$expected)

# Gridded Pearson residual plot.
pdf('figures/dspatresiduals.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(matrix(dspat_resid_df$pearson, nrow = length(unique(dspat_resid_df$x)))),
        unique(dspat_resid_df$x), unique(dspat_resid_df$y)),
     main = 'Gridded Pearson Residuals')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
points(dspat_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Set up a projection from the SPDE representation to the event locations.
dspat_event_proj <- inla.mesh.projector(dspat_mesh,
  loc = as.matrix(as.data.frame(dspat_ppp)))

# Create a copy of weeds_ppp and mark with 1/sqrt(lambda).
dspat_marked <- dspat_ppp
marks(dspat_marked) <- as.vector(1/sqrt(exp(
  dspat_result$summary.fixed['dist', 'mean'] *
    inla.mesh.project(dspat_event_proj, dspat_mesh_dist) +
  dspat_result$summary.fixed['intercept', 'mean'] +
  dspat_result$summary.fixed['river', 'mean'] *
    inla.mesh.project(dspat_event_proj, dspat_mesh_river) +
  dspat_result$summary.fixed['habitat', 'mean'] *
    inla.mesh.project(dspat_event_proj, dspat_mesh_habitat) +
  dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
    inla.mesh.project(dspat_event_proj, dspat_mesh_habitat^2) +
  inla.mesh.project(dspat_event_proj, dspat_result$summary.random$idx$mean)
)))

# Mark plot.
pdf('figures/dspatmarkplot.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(sqrt(exp(
          dspat_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_dist) +
          dspat_result$summary.fixed['intercept', 'mean'] +
          dspat_result$summary.fixed['river', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_river) +
          dspat_result$summary.fixed['habitat', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat) +
          dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
            inla.mesh.project(dspat_proj, dspat_mesh_habitat^2) +
          inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$mean)
        ))),
        xrange = site_win$x,
        yrange = site_win$y),
        riblab = expression(sqrt(hat(lambda))), ribsep = 0.05,
        main = 'Mark Plot')
plot(obs_win, border = NA, col = '#ffffff80', add = TRUE)
plot(dspat_marked, col = '#ffffff80', add = TRUE)
dev.off()

# Lurking variable plots.
lambda_im <- im(t(exp(
  dspat_result$summary.fixed['dist', 'mean'] *
    inla.mesh.project(dspat_proj, dspat_mesh_dist) +
  dspat_result$summary.fixed['intercept', 'mean'] +
  dspat_result$summary.fixed['river', 'mean'] *
    inla.mesh.project(dspat_proj, dspat_mesh_river) +
  dspat_result$summary.fixed['habitat', 'mean'] *
    inla.mesh.project(dspat_proj, dspat_mesh_habitat) +
  dspat_result$summary.fixed['I(habitat^2)', 'mean'] *
    inla.mesh.project(dspat_proj, dspat_mesh_habitat^2) +
  inla.mesh.project(dspat_proj, dspat_result$summary.random$idx$mean)
  )),
  xrange = site_win$x,
  yrange = site_win$y)

h_im <- 1/sqrt(lambda_im)

#TODO: lurking variable plots for river and habitat.

x_im <- im(matrix(seq(0, 100, 0.5), nrow = 201, ncol = 201, byrow = TRUE),
  xrange = site_win$x, yrange = site_win$y)

y_im <- im(matrix(seq(0, 100, 0.5), nrow = 201, ncol = 201, byrow = FALSE),
  xrange = site_win$x, yrange = site_win$y)

dist_im <- im(t(matrix(
  apply(sapply(unique(DSpat.obs$label), function(l){
      return(dist2line(
        as.ppp(expand.grid(x = seq(0, 100, 0.5), y = seq(0, 100, 0.5)), site_win),
        DSpat.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance
      )
    }), 1, min), nrow = 201)),
  xrange = site_win$x, yrange = site_win$y)

cum_pearson_x <- do.call(rbind, c(
  data.frame(x = 0, observed = 0, expected = 0, pearson = 0, v = 0, upper = 0, lower = 0),
  lapply(seq(0, 100, 0.5), function(x){
    sub_win <- obs_win[x_im <= x]
    sub_ppp <- dspat_ppp[sub_win]
    a <- area(sub_win)
    v <- a / area(obs_win)
    observed <- sub_ppp$n
    if(a > 0){
      expected <- mean(lambda_im[sub_win]) * a
      pearson <- (observed - expected) / sqrt(expected)
    } else {
      expected <- 0
      pearson <- 0
    }
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(x, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/dspatlurkx.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ x, data = cum_pearson_x, type = 'l',
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Horizontal Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ x, data = cum_pearson_x, type = 'l', lty = 3)
lines(upper ~ x, data = cum_pearson_x, type = 'l', lty = 3)
dev.off()

cum_pearson_y <- do.call(rbind, c(
  data.frame(y = 0, observed = 0, expected = 0, pearson = 0, v = 0, upper = 0, lower = 0),
  lapply(seq(0, 100, 0.5), function(y){
    sub_win <- obs_win[y_im <= y]
    sub_ppp <- dspat_ppp[sub_win]
    a <- area(sub_win)
    v <- a / area(obs_win)
    observed <- sub_ppp$n
    if(a > 0){
      expected <- mean(lambda_im[sub_win]) * a
      pearson <- (observed - expected) / sqrt(expected)
    } else {
      expected <- 0
      pearson <- 0
    }
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(y, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/dspatlurky.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ y, data = cum_pearson_y, type = 'l',
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Vertical Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ y, data = cum_pearson_y, type = 'l', lty = 3)
lines(upper ~ y, data = cum_pearson_y, type = 'l', lty = 3)
dev.off()

cum_pearson_dist <- do.call(rbind, c(
  data.frame(d = 0, observed = 0, expected = 0, pearson = 0, v = 0, upper = 0, lower = 0),
  lapply(seq(0, 7, 0.1), function(d){
    sub_win <- obs_win[dist_im <= d]
    sub_ppp <- dspat_ppp[sub_win]
    a <- area(sub_win)
    v <- a / area(obs_win)
    observed <- sub_ppp$n
    if(a > 0){
      expected <- mean(lambda_im[sub_win]) * a
      pearson <- (observed - expected) / sqrt(expected)
    } else {
      expected <- 0
      pearson <- 0
    }
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(d, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/dspatlurkdist.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ d, data = cum_pearson_dist, type = 'l',
     xlab = 'Distance from Transect', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Distance')
abline(h = 0, lty = 2)
lines(lower ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
lines(upper ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
dev.off()
