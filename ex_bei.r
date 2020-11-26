#######################
## PACKAGES AND DATA ##
#######################

library(spatstat) # For bei dataset.
library(INLA) # For model fitting.

bei_win <- Window(bei)

# Plot the observed point pattern and its window.
pdf('figures/bei.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei, main = expression(bold(paste(italic('Beilschmiedia pendula Lauraceae'), ' Locations'))))
dev.off()

# Plot the elevation surface.
pdf('figures/bei_elev.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(bei.extra$elev, main = 'Elevation', riblab = 'Meters', ribsep = 0.05)
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the gradient surface.
pdf('figures/bei_grad.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(bei.extra$grad, main = 'Gradient', riblab = 'm/m', ribsep = 0.05)
points(bei, pch = '.', col = '#ffffff80')
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

# Define a SPDE representation of a GP with Matern covariance.
bei_spde <- inla.spde2.pcmatern(
  mesh = bei_mesh,
  alpha = 2,
  prior.range = c(5, 0.1), # Pr(range < 5) = 0.1
  prior.sigma = c(2, 0.1) # Pr(sd > 2) = 0.1
)
# Note: alpha = 3/2 is exponential covariance. Only integer alpha are implemented.
# Moeller and Waagepetersen:
#   Exponential covariance
#   alpha is one-third of the range
#     log(alpha) ~ unif(log(1), log(235))
#   sigma is the sd
#     sigma ~ unif(0.001, inf)

# Set up a projection from the SPDE representation to a 400x200 grid.
bei_proj <- inla.mesh.projector(bei_mesh, dims = c(400, 200))


# Plot the mesh and point pattern.
pdf('figures/bei_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_mesh, asp = 1, main = '')
points(bei_mesh$loc[,1:2], pch = 20, col = 'black')
points(bei, pch = 20, col = 'red')
title('Mesh Over Bei Data')
dev.off()


################################
## PROJECT COVARIATES TO MESH ##
################################

# The covariates are measured on a lattice which does not match up to the mesh
# we'll use for the intensity function. INLA does support nested models where
# the prediction from one model (e.g. for a covariate) serves as a predictor
# in another model (e.g. the log-linear model for the point process intensity).
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
pdf('figures/bei_elev_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_elev)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        riblab = 'Meters', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Elevation')
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the gradient surface.
pdf('figures/bei_grad_mesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_grad)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        riblab = 'm/m',
        ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Gradient')
points(bei, pch = '.', col = '#ffffff80')
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
# node coordinates are the basis vectors, this is an identity matrix.
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
  idx = diag(bei_mesh_size)
)

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
pdf('figures/bei_mean.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Predicted Mean of Latent GP')
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
# A better way would be to use linear combinations to predict the
# latent surface at each pixel instead of interpolating the SD.
pdf('figures/bei_sd.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$sd)),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Prediction SD of Latent GP')
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/bei_betas.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
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
                         E(beta[2]*'|'*bold(x)) * z[2](u)), ribsep = 0.05,
     main = 'Posterior Mean of Fixed Effects')
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/bei_intensity.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
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
     riblab = 'Events per Square Meter', ribsep = 0.05,
     main = 'Posterior Intensity Function')
points(bei, pch = '.', col = '#ffffff80')
dev.off()

# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/bei_post.pdf', width = 18, height = 6)
par(mfrow = c(1, 3), bty = 'n')
plot(inla.smarginal(bei_result$marginals.hyperpar$`Stdev for idx`), type = 'l',
     yaxt = 'n', xlab = expression(sigma),
     main = 'Posterior Distribution of GP Standard Deviation')
plot(inla.smarginal(bei_result$marginals.hyperpar$`Range for idx`), type = 'l',
     yaxt = 'n', xlab = expression(rho),
     main = 'Posterior Distribution of the Range')
plot(inla.smarginal(bei_result$marginals.fixed$intercept), type = 'l',
     yaxt = 'n', xlab = expression(beta[0]),
     main = 'Posterior Distribution of the Intercept')
dev.off()

# Plot posterior marginals of the coefficients.
pdf('figures/bei_post_coefs.pdf', width = 12, height = 6)
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

# Have spatstat produce quadrat counts.
bei_qcounts <- t(quadratcount(
  bei,
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
pdf('figures/bei_residuals.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(bei_resid_tess, border = NA, do.col = TRUE,
     values = as.hyperframe(bei_resid_df)$pearson,
     ribargs = list(lty = 1, box = TRUE),
     ribsep = 0.05, main = 'Gridded Pearson Residuals')
points(bei, pch = '.', col = '#ffffff80')
mtext('Pearson Residual', 4)
dev.off()

# Set up a projection from the SPDE representation to the event locations.
bei_event_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(as.data.frame(bei)))

# Create a copy of bei and mark with 1/sqrt(lambda).
bei_marked <- bei
marks(bei_marked) <- as.vector(1/sqrt(inla.mesh.project(bei_event_proj,
  sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp))))

# Mark plot.
pdf('figures/bei_markplot.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(sqrt(inla.mesh.project(bei_proj,
  sapply(bei_result$marginals.lincomb.derived, inla.emarginal, fun = exp)))),
        xrange = Frame(bei)$x,
        yrange = Frame(bei)$y,
        unitname = c('meter', 'meters')),
        riblab = 'Square-Root Intensity',
        ribsep = 0.05, main = 'Mark Plot')
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

cum_pearson_x <- do.call(rbind, c(
  list(data.frame(x = 0, observed = 0, expected = 0, pearson = 0, a = 0, v = 0, upper = 0, lower = 0)),
  lapply(seq(0, 1000, 5), function(x){
    sub_win <- bei_win[x_im <= x]
    a <- area(sub_win)
    v <- a / area(bei_win)
    observed <- sum(bei$x <= x)
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

pdf('figures/bei_lurk_x.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ x, data = cum_pearson_x, type = 'l', ylim = c(-2, 2),
     xlab = 'Horizontal Coordinate (m)', ylab = 'Cumulative Pearson Residual',
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
    observed <- sum(bei$y <= y)
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

pdf('figures/bei_lurk_y.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ y, data = cum_pearson_y, type = 'l', ylim = c(-2, 2),
     xlab = 'Vertical Coordinate (m)', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Vertical Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ y, data = cum_pearson_y, type = 'l', lty = 3)
lines(upper ~ y, data = cum_pearson_y, type = 'l', lty = 3)
dev.off()

cum_pearson_elev <- do.call(rbind,
  lapply(seq(118, 160, 2), function(elev){
    sub_win <- bei_win[bei.extra$elev <= elev]
    sub_ppp <- bei[sub_win]
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

pdf('figures/bei_lurk_elev.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ elev, data = cum_pearson_elev, type = 'l', ylim = c(-2, 2),
     xlab = 'Elevation (m)', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Elevation')
abline(h = 0, lty = 2)
lines(lower ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
lines(upper ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
dev.off()

cum_pearson_grad <- do.call(rbind,
  lapply(seq(0, 0.35, 0.01), function(grad){
    sub_win <- bei_win[bei.extra$grad <= grad]
    sub_ppp <- bei[sub_win]
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

pdf('figures/bei_lurk_grad.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ grad, data = cum_pearson_grad, type = 'l', ylim = c(-3.5, 2),
     xlab = 'Gradient (m/m)', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Gradient')
abline(h = 0, lty = 2)
lines(lower ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
lines(upper ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
dev.off()
