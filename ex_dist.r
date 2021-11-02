#######################
## PACKAGES AND DATA ##
#######################

library(spatstat) # For bei dataset.
library(INLA) # For model fitting.

# Utility function to compute the distance between each point in a ppp
# and each segment in a psp.
distfromxsect <- function(events, xsects){
  if(!is.ppp(events)){
    stop('"events" must be a ppp object.')
  }
  if(!is.psp(xsects)){
    stop('"xsects" must be a psp object.')
  }
  evts <- cbind(events$x, events$y)
  x0 <- xsects$ends$x0
  x1 <- xsects$ends$x1
  y0 <- xsects$ends$y0
  y1 <- xsects$ends$y1
  x1x0 <- x1 - x0
  y1y0 <- y1 - y0
  apply(evts, 1, function(xy){
    return(min(ifelse(x0 == x1 & y0 == y1,
      sqrt((xy[1] - x0)^2 + (xy[2] - y0)^2),
      abs(y1y0 * xy[1] - x1x0 * xy[2] + x1 * y0 - y1 * x0) /
        sqrt((x1 - x0)^2 + (y1 - y0)^2)
    )))
  })
}

# Simulate sampling along a systematic sample of 20 parallel transects.
# The transects run north-south with a 50 meter spacing.
# Events within 5 meters of a transect will be recorded.
set.seed(-3264)
bei_win <- Window(bei)
xsect_width <- 10
bei_lines <- data.frame(
  x0 = runif(1, xsect_width/2, 50 - xsect_width/2) + 50 * (0:19),
  y0 = min(bei_win$y),
  y1 = max(bei_win$y)
)
bei_lines$x1 <- bei_lines$x0
bei_psp <- as.psp(bei_lines, window = bei_win)

# Mark the points with the distance to the nearest transect.
marks(bei) <- distfromxsect(bei, bei_psp)

# Observe events with probability exp(-0.1 * distance^2).
bei_obs <- bei[as.logical(rbinom(bei$n, 1, exp(-0.1 * marks(bei)^2)))]

# Plot the observed point pattern and its window.
pdf('figures/bei-dist.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_win, border = 'grey', main = expression(bold(paste('Observed ', italic('Beilschmiedia pendula Lauraceae'), ' Locations'))))
plot(bei_psp, add = TRUE)
points(bei_obs, col = '#00000080')
dev.off()

# Plot the elevation surface.
pdf('figures/bei-dist_elev.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$elev, main = 'Elevation', riblab = 'Meters', ribsep = 0.05)
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the gradient surface.
pdf('figures/bei-dist_grad.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$grad, main = 'Gradient', rib.sep = 0.05)
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()


#######################
## MESH CONSTRUCTION ##
#######################

# Get the boundary polygon for the observation window
# and define mesh edge segments for INLA.
bei_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices(bei_win)))

# Create a Delaunay triangulation with maximum edge length of 25 meters to use
# as the mesh.
bei_mesh <- inla.mesh.create(
  boundary = bei_boundary,
  refine = list(max.edge = 25)
)

# Convert the mesh edges to a psp.
meshloc <- bei_mesh$loc[,-3]
meshadj <- bei_mesh$graph$vv
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
bei_spde <- inla.spde2.pcmatern(
  mesh = bei_mesh,
  alpha = 2,
  prior.range = c(5, 0.1), # Pr(range < 5) = 0.1
  prior.sigma = c(2, 0.1) # Pr(sd > 2) = 0.1
)

# Set up a projection from the SPDE representation to a 400x200 grid.
bei_proj <- inla.mesh.projector(bei_mesh, dims = c(400, 200))


# Plot the mesh and point pattern.
pdf('figures/bei-dist_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_mesh, asp = 1, main = '')
plot(bei_psp, add = TRUE, col = 'red', border = NA)
points(bei_mesh$loc[,1:2], pch = 20, col = 'black')
points(bei_obs, pch = 20, col = 'red')
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

# Compute distance from transect for each node.
bei_mesh_dist <- distfromxsect(as.ppp(bei_mesh$loc[,1:2], bei_win), bei_psp)

# Plot the piecewise linear approximation of the elevation surface.
pdf('figures/bei-dist_elev_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_elev)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        riblab = 'Meters', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Elevation')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the gradient surface.
pdf('figures/bei-dist_grad_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_grad)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Gradient')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the distance surface.
pdf('figures/bei-dist_dist_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_dist)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        zlim = c(0, 35), ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Distance to Transect')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()


######################################
## SET UP SIMPSON METHOD PSEUDODATA ##
######################################

# Observed event locations.
bei_pts <- cbind(bei_obs$x, bei_obs$y)

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
bei_formula <- y ~ -1 + intercept + elev + grad + I(dist^2) + f(idx, model = bei_spde)

# Define linear combinations to predict the posterior distribution of the linear
# predictor for the log-intensity surface at each node. We will use these for
# model checking, exponentiating them and projecting to different lattices as
# needed. A much more accurate overall approach is to specify a lincomb for
# every single combination of posterior quantity and projection before fitting
# the model (at added computational cost) but this quick-and-dirty approach is
# adequate for checking model fit.
bei_lcs <- inla.make.lincombs(
  intercept = rep(1, bei_mesh_size),
  elev = bei_mesh_elev,
  grad = bei_mesh_grad,
  `I(dist^2)` = bei_mesh_dist^2,
  idx = diag(bei_mesh_size)
)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
bei_inla_data <- list(
  y = bei_pseudodata, # The whole pseudodata vector.
  elev = bei_mesh_elev, # Elevation at nodes.
  grad = bei_mesh_grad, # Gradient at nodes.
  dist = bei_mesh_dist, # Distance at nodes.
  idx = seq_len(bei_mesh_size), # Indices of the nodes.
  intercept = rep(1, bei_mesh_size) # Intercept column.
)

# Fit the model as a Poisson GLM with exposures specified.
bei_result <- inla(
  formula = bei_formula,
  data = bei_inla_data,
  family = 'poisson',
  control.fixed = list(
    # Prior means and precisions for coefficients.
    mean.intercept = 0,
    prec.intercept = 0,
    mean = 0,
    prec = 0
  ),
  control.predictor = list(A = bei_pseudopoints),
  E = bei_pseudodata_exp,
  lincomb = bei_lcs
)

# Summarize the posterior marginals of the parameters.
print(bei_result$summary.fixed)
print(bei_result$summary.hyperpar)

# Plot the posterior mean of the latent surface.
pdf('figures/bei-dist_mean.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Predicted Mean of Latent GP')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
# A better way would be to use linear combinations to predict the
# latent surface at each pixel instead of interpolating the SD.
pdf('figures/bei-dist_sd.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$sd)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Prediction SD of Latent GP')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/bei-dist_betas.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad)
        ),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(beta[0]*'|'*bold(x)) +
                         E(beta[1]*'|'*bold(x)) * z[1](u) +
                         E(beta[2]*'|'*bold(x)) * z[2](u)), ribsep = 0.05,
     main = 'Posterior Mean of Fixed Effects')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the detecton function.
pdf('figures/bei-dist_detection.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          bei_result$summary.fixed['I(dist^2)', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_dist^2)
        )),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     zlim = 0:1, riblab = 'Detection Probability', ribsep = 0.05,
     main = 'Posterior Detection Surface')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/bei-dist_intensity.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(exp(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad) +
          inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)
        )),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = 'Events per Square Meter', ribsep = 0.05,
     main = 'Posterior Intensity Function')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface
# for the thinned process.
pdf('figures/bei-dist_thinned_intensity.pdf', width = 7, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj,
            sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp))
        ),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = 'Observable Events per Square Meter', ribsep = 0.05,
     main = 'Posterior Thinned Intensity Function')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()


# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/bei-dist_post.pdf', width = 12, height = 6)
par(mfrow = c(1, 2), bty = 'n')
plot(inla.smarginal(bei_result$marginals.hyperpar$`Stdev for idx`), type = 'l',
     yaxt = 'n', xlab = expression(sigma),
     main = 'Posterior Distribution of GP Standard Deviation')
plot(inla.smarginal(bei_result$marginals.hyperpar$`Range for idx`), type = 'l',
     yaxt = 'n', xlab = expression(rho),
     main = 'Posterior Distribution of the Range')
dev.off()

# Plot posterior marginals of the coefficients.
pdf('figures/bei-dist_post_coefs.pdf', width = 12, height = 12)
par(mfrow = c(2, 2), bty = 'n')
plot(inla.smarginal(bei_result$marginals.fixed$intercept), type = 'l',
     yaxt = 'n', xlab = expression(beta[0]),
     main = 'Posterior Distribution of the Intercept')
plot(inla.smarginal(bei_result$marginals.fixed$elev), type = 'l',
     yaxt = 'n', xlab = expression(beta[1]),
     main = 'Posterior Distribution of Elevation Coefficient')
plot(inla.smarginal(bei_result$marginals.fixed$grad), type = 'l',
     yaxt = 'n', xlab = expression(beta[2]),
     main = 'Posterior Distribution of Gradient Coefficient')
plot(inla.smarginal(bei_result$marginals.fixed$`I(dist^2)`), type = 'l',
     yaxt = 'n', xlab = expression(beta[3]),
     main = 'Posterior Distribution of the Distance Coefficient')
dev.off()

# Plot the detection function.
pdf('figures/bei-dist_detect_curve.pdf', width = 7, height = 6)
par(bty = 'n')
curve(exp(bei_result$summary.fixed['I(dist^2)', 'mean'] * x^2),
      from = 0, to = 25, ylim = 0:1,,
      ylab = 'Detection Probability', xlab = 'Distance',
      main = 'Posterior Detection Function')
dev.off()

####################
## MODEL CHECKING ##
####################

# Create grid counts for residual calculation.

# Number of grid cells.
NGRID_X <- 40
NGRID_Y <- 20

# Have spatstat produce quadrat counts.
bei_qcounts <- t(quadratcount(
  bei_obs,
  xbreaks = seq(min(bei_win$xrange), max(bei_win$xrange), length.out = NGRID_X + 1),
  ybreaks = seq(min(bei_win$yrange), max(bei_win$yrange), length.out = NGRID_Y + 1),
))

# Get the tesselation created by quadratcount().
bei_resid_tess <- as.tess(bei_qcounts)

# Get the centroids of the tesselation tiles and put them in a data frame.
bei_resid_df <- do.call(rbind, lapply(lapply(tiles(bei_resid_tess), centroid.owin), data.frame))

# Add the counts to the data frame.
bei_resid_df$count <- as.numeric(bei_qcounts)

# Find the observed area of each cell.
bei_resid_df$area <- tile.areas(bei_resid_tess)

# Set up a projection from the SPDE representation to the residual grid.
bei_resid_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(bei_resid_df[,c('x', 'y')]))

# Project the intensity surface to the residual grid.
bei_resid_df$intensity <- as.vector(inla.mesh.project(bei_resid_proj,
  sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp)))

# Calculate the (approximate) expected counts and the Pearson residuals.
bei_resid_df$expected <- bei_resid_df$intensity * bei_resid_df$area
bei_resid_df$pearson <- (bei_resid_df$count - bei_resid_df$expected) /
  sqrt(bei_resid_df$expected)

# Gridded Pearson residual plot.
pdf('figures/bei-dist_residuals.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(bei_resid_tess, border = NA, do.col = TRUE,
     values = as.hyperframe(bei_resid_df)$pearson,
     xlim = bei_win$xrange, ylim = bei_win$yrange,
     ribargs = list(lty = 1, box = TRUE),
     ribsep = 0.05, main = 'Gridded Pearson Residuals')
plot(bei_psp, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = '.', col = '#ffffff80')
dev.off()

# Set up a projection from the SPDE representation to the event locations.
bei_event_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(as.data.frame(bei_obs)))

# Create a copy of bei and mark with 1/sqrt(lambda).
bei_marked <- bei_obs
marks(bei_marked) <- as.vector(1/sqrt(inla.mesh.project(bei_event_proj,
  sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp))))

# Mark plot.
pdf('figures/bei-dist_markplot.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(sqrt(inla.mesh.project(bei_proj,
  sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp)))),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        ribsep = 0.05, main = 'Mark Plot')
plot(bei_psp, col = '#ffffff80', add = TRUE)
plot(bei_marked, col = '#ffffff80', add = TRUE)
dev.off()

# Lurking variable plots.
lambda_im <- im(t(inla.mesh.project(bei_proj, sapply(
    bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp
  ))),
  xrange = bei_win$x,
  yrange = bei_win$y)

h_im <- 1/sqrt(lambda_im)

x_im <- im(matrix(seq(0, 1000, 5), nrow = 101, ncol = 201, byrow = TRUE),
  xrange = bei_win$x, yrange = bei_win$y)

y_im <- im(matrix(seq(0, 500, 5), nrow = 101, ncol = 201, byrow = FALSE),
  xrange = bei_win$x, yrange = bei_win$y)

dist_im <- im(matrix(distfromxsect(as.ppp(expand.grid(x = 0:1000, y = 0:500), bei_win), bei_psp), ncol = 1001, byrow = TRUE),
  xrange = bei_win$x, yrange = bei_win$y)

cum_pearson_x <- do.call(rbind, c(
  list(data.frame(x = 0, observed = 0, expected = 0, pearson = 0, a = 0, v = 0, upper = 0, lower = 0)),
  lapply(seq(0, 1000, 5), function(x){
    sub_win <- bei_win[x_im <= x]
    a <- area(sub_win)
    v <- a / area(bei_win)
    observed <- sum(bei_obs$x <= x)
    if(a > 0){
      expected <- mean(lambda_im[sub_win]) * a
      pearson <- (observed - expected) / sqrt(expected)
    } else {
      expected <- 0
      pearson <- 0
    }
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(x, observed, expected, pearson, a, v, upper, lower))
  })
))

pdf('figures/bei-dist_lurk_x.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ x, data = cum_pearson_x, type = 'l', ylim = c(-4, 2),
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Horizontal Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ x, data = cum_pearson_x, type = 'l', lty = 3)
lines(upper ~ x, data = cum_pearson_x, type = 'l', lty = 3)
dev.off()

cum_pearson_y <- do.call(rbind, c(
  list(data.frame(y = 0, observed = 0, expected = 0, pearson = 0, a = 0, v = 0, upper = 0, lower = 0)),
  lapply(seq(0, 500, 5), function(y){
    sub_win <- bei_win[y_im <= y]
    a <- area(sub_win)
    v <- a / area(bei_win)
    observed <- sum(bei_obs$y <= y)
    if(a > 0){
      expected <- mean(lambda_im[sub_win]) * a
      pearson <- (observed - expected) / sqrt(expected)
    } else {
      expected <- 0
      pearson <- 0
    }
    upper <- 2 * sqrt(v)
    lower <- -2 * sqrt(v)
    return(data.frame(y, observed, expected, pearson, a, v, upper, lower))
  })
))

pdf('figures/bei-dist_lurk_y.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ y, data = cum_pearson_y, type = 'l', ylim = c(-3, 2),
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Vertical Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ y, data = cum_pearson_y, type = 'l', lty = 3)
lines(upper ~ y, data = cum_pearson_y, type = 'l', lty = 3)
dev.off()

cum_pearson_dist <- do.call(rbind,
  lapply(seq(0, 33, 0.5), function(d){
    sub_win <- bei_win[dist_im <= d]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_win)
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
    return(data.frame(d, observed, expected, pearson, a, v, upper, lower))
  })
)

pdf('figures/bei-dist_lurk_dist.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ d, data = cum_pearson_dist, type = 'l',
     xlab = 'Distance', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Distance')
abline(h = 0, lty = 2)
lines(lower ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
lines(upper ~ d, data = cum_pearson_dist, type = 'l', lty = 3)
dev.off()

cum_pearson_elev <- do.call(rbind,
  lapply(seq(118, 160, 2), function(elev){
    sub_win <- bei_win[bei.extra$elev <= elev]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_win)
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
    return(data.frame(elev, observed, expected, pearson, a, v, upper, lower))
  })
)

pdf('figures/bei-dist_lurk_elev.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ elev, data = cum_pearson_elev, type = 'l', ylim = c(-2, 2),
     xlab = 'Elevation', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Elevation')
abline(h = 0, lty = 2)
lines(lower ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
lines(upper ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
dev.off()

cum_pearson_grad <- do.call(rbind,
  lapply(seq(0, 0.35, 0.01), function(grad){
    sub_win <- bei_win[bei.extra$grad <= grad]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_win)
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
    return(data.frame(grad, observed, expected, pearson, a, v, upper, lower))
  })
)

pdf('figures/bei-dist_lurk_grad.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ grad, data = cum_pearson_grad, type = 'l', ylim = c(-2, 2),
     xlab = 'Gradient', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Gradient')
abline(h = 0, lty = 2)
lines(lower ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
lines(upper ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
dev.off()
