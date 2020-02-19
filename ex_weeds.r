#######################
## PACKAGES AND DATA ##
#######################

library(spatstat)
library(INLA) # For model fitting.
library(DSpat) # For dataset.

data(weeds.obs) # Incompletely observed (distance sampling) dataset.
data(weeds.lines) # Transect lines used for distance sampling.
data(weeds.covariates) # Lattice of predictor variables.

# Create spatstat objects for plotting.
weeds_win <- owin(c(0, 1200), c(0, 1200), unitname = c('meter', 'meters'))
weeds_psp <- psp(weeds.lines$x0, weeds.lines$y0, weeds.lines$x1, weeds.lines$y1, weeds_win)
weeds_ppp <- ppp(weeds.obs$x, weeds.obs$y, window = weeds_win)
weeds_covs_ppp <- as.ppp(weeds.covariates, weeds_win)
weeds_sheep_im <- im(
  matrix(weeds.covariates$sheep, ncol = length(unique(weeds.covariates$x))),
  xcol = unique(weeds.covariates$x), yrow = unique(weeds.covariates$y)
)
weeds_dist_im <- im(t(matrix(
  apply(sapply(unique(weeds.covariates$strip), function(l){
      return(dist2line(
        as.ppp(expand.grid(x = 0:1200, y = seq(0, 1200, 100)), weeds_win),
        weeds.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance
      )
    }), 1, min), nrow = 1201)),
  xrange = weeds_win$x, yrange = weeds_win$y)


# Plot the observed point pattern, its window, and the transects.
pdf('figures/weeds.pdf', width = 6, height = 6)
par(mar = c(3, 3, 2, 0))
plot(weeds_ppp, main = 'Devil\'s Claw Locations', border = 'grey')
plot(weeds_psp, add = TRUE)
axis(1, at = seq(0, 1200, 300))
mtext('Horizontal Coordinate', 1, 2)
axis(2, at = seq(0, 1200, 300))
mtext('Vertical Coordinate', 2, 2)
dev.off()

# Plot the sheep covariate and the transects.
pdf('figures/weedssheep.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(weeds_sheep_im, main = 'Presence of Sheep', border = 'grey',
     riblab = 'Sheep Indicator', ribsep = 0.05, ribn = 2,
     ribargs = list(at = c(300, 900), labels = 0:1))
plot(weeds_psp, add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the distance-to-transect surface.
pdf('figures/weedsdist.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(weeds_dist_im, riblab = 'Meters', ribsep = 0.05,
     main = 'Distance to Transect')
plot(weeds_win, border = '#80808080', add = TRUE)
plot(weeds_psp, add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


#######################
## MESH CONSTRUCTION ##
#######################

# Get the boundary polygon for the observation window
# and define mesh edge segments for INLA.
weeds_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(weeds_win)))

# Define edges for a finer mesh over the discontinuity in sheep presence.
# This is a 20m strip running the vertical length of the site. It should really
# be narrower, but we are keeping it wide enough so the continuitiy in the
# approximation is visible in the plots.
weeds_refine <- inla.mesh.segment(loc = cbind(
  c(590, 610, 610, 590), c(0, 0, 1200, 1200)
))

# Create the finer mesh.
weeds_refine_mesh <- inla.mesh.create(
  boundary = weeds_refine,
  refine = list(max.edge = 50)
)

# Create a Delaunay triangulation with maximum edge length of 50 meters to use
# as the mesh. Include the nodes from the finer mesh
weeds_mesh <- inla.mesh.create(
  loc = weeds_refine_mesh$loc[,1:2],
  boundary = weeds_boundary,
  refine = list(max.edge = 50)
)

# Define a SPDE representation of a GP with Matern covariance.
weeds_spde <- inla.spde2.matern(mesh = weeds_mesh)

# Set up a projection from the SPDE representation to a 400x400 grid.
weeds_proj <- inla.mesh.projector(weeds_mesh, dims = c(400, 400))


# Plot the mesh and point pattern.
pdf('figures/weedsmesh.pdf', width = 6, height = 6)
par(mar = c(0, 0, 2, 0))
plot(weeds_mesh, asp = 1, main = '')
plot(weeds_psp, add = TRUE, col = 'red')
points(weeds_ppp, pch = 20, col = 'red')
title('Mesh Over Devil\'s Claw Data')
dev.off()


################################
## PROJECT COVARIATES TO MESH ##
################################

# Create a ppp version of the mesh nodes to help with calculations.
weeds_mesh_ppp <- as.ppp(weeds_mesh$loc[,1:2], weeds_win)

# We are assuming the sheep indicator is exactly known and has the form
# sheep(x, y) = 1 if x > 600; 0 otherwise. Calculate this at each mesh node
# instead of interpolating or modeling.
weeds_mesh_sheep <- as.numeric(weeds_mesh_ppp$x > 600)

# Get the distance from each node to the nearest transect.
weeds_mesh_dist <- apply(sapply(unique(weeds.covariates$strip), function(l){
  return(dist2line(weeds_mesh_ppp, weeds.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance)
}), 1, min)

# Plot the piecewise linear approximation of the sheep presence surface.
pdf('figures/weedssheepmesh.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(weeds_proj, weeds_mesh_sheep)),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
        riblab = 'Sheep Indicator', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Sheep Presence')
plot(weeds_win, border = '#80808080', add = TRUE)
plot(weeds_psp, add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the distance surface.
pdf('figures/weedsdistmesh.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(weeds_proj, weeds_mesh_dist)),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
        riblab = 'Meters', ribsep = 0.05, zlim = c(0, 75),
        main = 'Piecewise Linear Approximation of Distance to Transect')
plot(weeds_win, border = '#80808080', add = TRUE)
plot(weeds_psp, add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


######################################
## SET UP SIMPSON METHOD PSEUDODATA ##
######################################

# Observed event locations.
weeds_pts <- cbind(weeds.obs$x, weeds.obs$y)

# Get the numbers of mesh nodes and real events.
# The sum of these will be the number of pseudodata points.
weeds_mesh_size <- weeds_mesh$n
weeds_n_events <- nrow(weeds_pts)

# Create the psuedodata. This is a vector giving the count of events at each
# pseudopoint location, that is 0 at each mesh node and 1 at each event.
weeds_pseudodata <- c(rep(0, weeds_mesh_size), rep(1, weeds_n_events))

# Get the numerical integration weights for the SPDE approach.
# Because the whole region was fully observed (constant sampling effort),
# these are left as is. If parts of the region were unobserved, these would be
# multiplied by the observed proportion of the area represented by the node.
weeds_int_weights <- diag(inla.mesh.fem(weeds_mesh)$c0)

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
weeds_pseudodata_exp <- c(weeds_int_weights, rep(0, weeds_n_events))

# Compute the barycentric coordinates of the observed events
# (i.e. project into the space spanned by the basis of mesh nodes).
weeds_bary <- inla.mesh.project(weeds_mesh, weeds_pts)$A

# Compute the barycentric coordinates of the nodes. Because the
# node coordinatess are the basis vectors, this is an identity matrix.
weeds_int_matrix <- sparseMatrix(
  i = seq_len(weeds_mesh_size),
  j = seq_len(weeds_mesh_size),
  x = rep(1, weeds_mesh_size)
)

# Bind the node and event coordinates into a single matrix of pseudodata
# locations in barycentric coordinates.
weeds_pseudopoints <- rbind(weeds_int_matrix, weeds_bary)


########################################
## FIT THE MODEL AND PLOT THE RESULTS ##
########################################

# Define the model formula.
# inla() requires a vector of 1s in the data argument for the intercept.
# The f() term specifies the GP with observations indexed by a variable called
# idx. The indices correspond to the indixes of the mesh nodes.
weeds_formula <- y ~ -1 + intercept + sheep * dist + f(idx, model = weeds_spde)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
weeds_inla_data <- list(
  y = weeds_pseudodata, # The whole pseudodata vector.
  sheep = weeds_mesh_sheep, # Sheep presence at nodes.
  dist = weeds_mesh_dist, # Distance at nodes.
  idx = seq_len(weeds_mesh_size), # Indices of the nodes.
  intercept = rep(1, weeds_mesh_size) # Intercept column.
)

# Fit the model as a Poisson GLM with exposures specified.
weeds_result <- inla(
  formula = weeds_formula,
  data = weeds_inla_data,
  family = 'poisson',
  control.predictor = list(A = weeds_pseudopoints),
  E = weeds_pseudodata_exp
)

# Plot the posterior mean of the latent surface.
pdf('figures/weedsmean.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$mean)),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Predicted Mean of Latent GP')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
pdf('figures/weedssd.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$sd)),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Prediction SD of Latent GP')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/weedsbetas.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(
          weeds_result$summary.fixed['intercept', 'mean'] +
          weeds_result$summary.fixed['sheep', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep)
        ),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(beta[0]*'|'*bold(x)) +
                         E(beta[1]*'|'*bold(x)) * z[1](u)), ribsep = 0.05,
     main = 'Posterior Mean of Fixed Effects')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the detection surface.
pdf('figures/weedsdetection.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          weeds_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_dist) +
          weeds_result$summary.fixed['sheep:dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep * weeds_mesh_dist)
        )),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
     zlim = 0:1, riblab = 'Detection Probability', ribsep = 0.05,
     main = 'Posterior Detection Surface')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/weedsintensity.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          weeds_result$summary.fixed['intercept', 'mean'] +
          weeds_result$summary.fixed['sheep', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep) +
          inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$mean)
        )),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')) * 1000000,
     riblab = 'Events per Square Kilometer', ribsep = 0.05,
     main = 'Posterior Intensity Function')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface
# for the thinned process.
pdf('figures/weedsthinnedintensity.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          weeds_result$summary.fixed['intercept', 'mean'] +
          weeds_result$summary.fixed['sheep', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep) +
          weeds_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_dist) +
          weeds_result$summary.fixed['sheep:dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep * weeds_mesh_dist) +
          inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$mean)
        )),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')) * 1000000,
     riblab = 'Observable Events per Square Kilometer', ribsep = 0.05,
     main = 'Posterior Thinned Intensity Function')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()


# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/weedspost.pdf', width = 18, height = 6)
par(mfrow = c(1, 3), bty = 'n')
plot(inla.smarginal(weeds_result$marginals.hyperpar$Theta1), type = 'l',
     yaxt = 'n', xlab = expression(tau),
     main = expression(paste(bold('Posterior Distribution of '), tau)))
plot(inla.smarginal(weeds_result$marginals.hyperpar$Theta2), type = 'l',
     yaxt = 'n', xlab = expression(kappa),
     main = expression(paste(bold('Posterior Distribution of '), kappa)))
plot(inla.smarginal(weeds_result$marginals.fixed$intercept), type = 'l',
     yaxt = 'n', xlab = expression(beta[0]),
     main = 'Posterior Distribution of the Intercept')
dev.off()

# Plot posterior marginals of the coefficients.
pdf('figures/weedspostcoefs.pdf', width = 18, height = 6)
par(mfrow = c(1, 3), bty = 'n')
plot(inla.smarginal(weeds_result$marginals.fixed$sheep), type = 'l',
     yaxt = 'n', xlab = expression(beta[1]),
     main = 'Posterior Distribution of Sheep Coefficient')
plot(inla.smarginal(weeds_result$marginals.fixed$dist), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Distance Coefficient')
plot(inla.smarginal(weeds_result$marginals.fixed$`sheep:dist`), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Sheep x Distance Coefficient')
dev.off()

# Plot the detection functions.
pdf('figures/weedsdetectcurves.pdf', width = 7, height = 6)
par(bty = 'n')
curve(exp(weeds_result$summary.fixed['dist', 'mean'] * x),
      from = 0, to = 75, ylim = 0:1, lty = 1, xaxt = 'n',
      ylab = 'Detection Probability', xlab = 'Distance',
      main = 'Posterior Detection Function')
curve(exp(sum(weeds_result$summary.fixed[c('dist', 'sheep:dist'), 'mean']) * x),
      from = 0, to = 75, add = TRUE, lty = 4)
axis(1, at = seq(0, 75, 15))
legend(37.5, 1, lty = c(1, 4), legend = c('No Sheep', 'Sheep'),
       horiz = TRUE, bty = 'n', xjust = 0.5, yjust = 0.3)
dev.off()


####################
## MODEL CHECKING ##
####################

# Create grid counts for residual calculation.

# Number of grid cells.
NGRID_X <- 20
NGRID_Y <- 20

# Have spatstat find centers for the grid cells.
centers <- gridcenters(
  owin(range(weeds_boundary$loc[,1]), range(weeds_boundary$loc[,2])),
  NGRID_X, NGRID_Y)

# Compute the grid cell size.
dx <- sum(unique(centers$x)[1:2] * c(-1, 1)) / 2
dy <- sum(unique(centers$y)[1:2] * c(-1, 1)) / 2

# Initialize a data frame to store the counts.
weeds_resid_df <- data.frame(x = centers$x, y = centers$y,
                             count = NA_integer_, area = NA_real_)

# Loop through the cells, finding the event count and cell area.
for(r in seq_len(nrow(weeds_resid_df))){
  weeds_resid_df$count[r] <- sum(weeds_ppp$x >= weeds_resid_df$x[r] - dx &
                                 weeds_ppp$x < weeds_resid_df$x[r] + dx &
                                 weeds_ppp$y >= weeds_resid_df$y[r] - dy &
                                 weeds_ppp$y < weeds_resid_df$y[r] + dy)
  weeds_resid_df$area[r] <- area(weeds_win[owin(c(weeds_resid_df$x[r] - dx,
                                                  weeds_resid_df$x[r] + dx),
                                                c(weeds_resid_df$y[r] - dy,
                                                  weeds_resid_df$y[r] + dy))])
}

# Set up a projection from the SPDE representation to the residual grid.
weeds_resid_proj <- inla.mesh.projector(weeds_mesh,
  loc = as.matrix(as.data.frame(centers)))

# Project the intensity surface to the residual grid.
weeds_resid_df$intensity <- as.vector(exp(
  weeds_result$summary.fixed['intercept', 'mean'] +
    weeds_result$summary.fixed['sheep', 'mean'] *
  inla.mesh.project(weeds_resid_proj, weeds_mesh_sheep) +
    weeds_result$summary.fixed['dist', 'mean'] *
  inla.mesh.project(weeds_resid_proj, weeds_mesh_dist) +
    weeds_result$summary.fixed['sheep:dist', 'mean'] *
  inla.mesh.project(weeds_resid_proj, weeds_mesh_sheep * weeds_mesh_dist) +
  inla.mesh.project(weeds_resid_proj, weeds_result$summary.random$idx$mean)
))

# Calculate the (approximate) expected counts and the Pearson residuals.
weeds_resid_df$expected <- weeds_resid_df$intensity * weeds_resid_df$area
weeds_resid_df$pearson <- (weeds_resid_df$count - weeds_resid_df$expected) /
  sqrt(weeds_resid_df$expected)

# Gridded Pearson residual plot.
pdf('figures/weedsresiduals.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(matrix(weeds_resid_df$pearson, nrow = length(unique(weeds_resid_df$x)))),
        unique(weeds_resid_df$x), unique(weeds_resid_df$y),
        unitname = c('meter', 'meters')),
     main = 'Gridded Pearson Residuals')
plot(weeds_win, border = 'white', add = TRUE)
points(weeds_ppp, pch = 21, col = '#00000080', bg = '#ffffff80')
dev.off()

# Set up a projection from the SPDE representation to the event locations.
weeds_event_proj <- inla.mesh.projector(weeds_mesh,
  loc = as.matrix(as.data.frame(weeds_ppp)))

# Create a copy of weeds_ppp and mark with 1/sqrt(lambda).
weeds_marked <- weeds_ppp
marks(weeds_marked) <- as.vector(1/sqrt(exp(
  weeds_result$summary.fixed['intercept', 'mean'] +
    weeds_result$summary.fixed['sheep', 'mean'] *
  inla.mesh.project(weeds_event_proj, weeds_mesh_sheep) +
    weeds_result$summary.fixed['dist', 'mean'] *
  inla.mesh.project(weeds_event_proj, weeds_mesh_dist) +
    weeds_result$summary.fixed['sheep:dist', 'mean'] *
  inla.mesh.project(weeds_event_proj, weeds_mesh_sheep * weeds_mesh_dist) +
  inla.mesh.project(weeds_event_proj, weeds_result$summary.random$idx$mean)
)))

# Mark plot.
pdf('figures/weedsmarkplot.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(sqrt(exp(
          weeds_result$summary.fixed['intercept', 'mean'] +
          weeds_result$summary.fixed['sheep', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep) +
          weeds_result$summary.fixed['dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_dist) +
          weeds_result$summary.fixed['sheep:dist', 'mean'] *
            inla.mesh.project(weeds_proj, weeds_mesh_sheep * weeds_mesh_dist) +
          inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$mean)
        ))),
        xrange = weeds_win$x,
        yrange = weeds_win$y,
        unitname = c('meter', 'meters')),
        riblab = expression(sqrt(hat(lambda))), ribsep = 0.05,
        main = 'Mark Plot')
plot(weeds_win, border = 'white', add = TRUE)
plot(weeds_marked, col = '#ffffff80', add = TRUE)
dev.off()

# Lurking variable plots.
lambda_im <- im(t(exp(
    weeds_result$summary.fixed['intercept', 'mean'] +
    weeds_result$summary.fixed['sheep', 'mean'] *
      inla.mesh.project(weeds_proj, weeds_mesh_sheep) +
    weeds_result$summary.fixed['dist', 'mean'] *
      inla.mesh.project(weeds_proj, weeds_mesh_dist) +
    weeds_result$summary.fixed['sheep:dist', 'mean'] *
      inla.mesh.project(weeds_proj, weeds_mesh_sheep * weeds_mesh_dist) +
    inla.mesh.project(weeds_proj, weeds_result$summary.random$idx$mean)
  )),
  xrange = weeds_win$x,
  yrange = weeds_win$y,
  unitname = c('meter', 'meters'))

h_im <- 1/sqrt(lambda_im)

x_im <- im(matrix(seq(0, 1200, 5), nrow = 241, ncol = 241, byrow = TRUE),
  xrange = weeds_win$x, yrange = weeds_win$y)

y_im <- im(matrix(seq(0, 1200, 5), nrow = 241, ncol = 241, byrow = FALSE),
  xrange = weeds_win$x, yrange = weeds_win$y)

dist_im <- im(t(matrix(
  apply(sapply(unique(weeds.covariates$strip), function(l){
      return(dist2line(
        as.ppp(expand.grid(x = seq(0, 1200, 5), y = seq(0, 1200, 5)), weeds_win),
        weeds.lines[l, c('x0', 'y0', 'x1', 'y1')])$distance
      )
    }), 1, min), nrow = 241)),
  xrange = weeds_win$x, yrange = weeds_win$y)

cum_pearson_x <- do.call(rbind, c(
  data.frame(x = 0, observed = 0, expected = 0, pearson = 0, v = 0, upper = 0, lower = 0),
  lapply(seq(0, 1200, 5), function(x){
    sub_ppp <- weeds_ppp[x_im <= x]
    a <- area(sub_ppp)
    v <- a / area(weeds_ppp)
    observed <- sub_ppp$n
    expected <- mean(lambda_im[x_im <= x]) * a
    pearson <- (observed - expected) / sqrt(expected)
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(x, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/weedslurkx.pdf', width = 7, height = 6)
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
  lapply(seq(0, 1200, 5), function(y){
    sub_ppp <- weeds_ppp[y_im <= y]
    a <- area(sub_ppp)
    v <- a / area(weeds_ppp)
    observed <- sub_ppp$n
    expected <- mean(lambda_im[y_im <= y]) * a
    pearson <- (observed - expected) / sqrt(expected)
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(y, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/weedslurky.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ y, data = cum_pearson_y, type = 'l', ylim = c(-2, 2),
     xlab = 'Vertical Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Vertical Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ y, data = cum_pearson_y, type = 'l', lty = 3)
lines(upper ~ y, data = cum_pearson_y, type = 'l', lty = 3)
dev.off()

cum_pearson_dist <- do.call(rbind, c(
  data.frame(d = 0, observed = 0, expected = 0, pearson = 0, v = 0, upper = 0, lower = 0),
  lapply(0:75, function(d){
    sub_ppp <- weeds_ppp[dist_im <= d]
    a <- area(sub_ppp)
    v <- a / area(weeds_ppp)
    observed <- sub_ppp$n
    expected <- mean(lambda_im[dist_im <= d]) * a
    pearson <- (observed - expected) / sqrt(expected)
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(d, observed, expected, pearson, v, upper, lower))
  })
))

pdf('figures/weedslurkdist.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ d, data = cum_pearson_dist, type = 'l',
     xlab = 'Distance from Transect', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Distance')
abline(h = 0, lty = 2)
lines(lower ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
lines(upper ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
dev.off()
