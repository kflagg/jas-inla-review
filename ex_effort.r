#######################
## PACKAGES AND DATA ##
#######################

library(spatstat) # For bei dataset.
library(INLA) # For model fitting.

# Simulate sampling along a systematic sample of 15 parallel transects.
# The transects run north-south with a 68 meter spacing.
# Events within 5 meters of a transect will be recorded.
set.seed(-3264)
bei_win <- Window(bei)
bei_lines <- data.frame(
  x0 = runif(1, 0, 34) + 68 * (0:14),
  y0 = min(bei_win$y),
  y1 = max(bei_win$y)
)
bei_lines$x1 <- bei_lines$x0
bei_psp <- as.psp(bei_lines, window = bei_win)

# Create an window by giving the psp width.
bei_xsect <- dilation(bei_psp, 5)[bei_win]

# Subset the point pattern to the observation window.
bei_obs <- bei[bei_xsect]

# Plot the observed point pattern and its window.
pdf('figures/beieffort.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_win, border = 'grey', main = expression(paste('Observed ', italic('Beilschmiedia pendula Lauraceae'), ' Locations')))
plot(bei_xsect, add = TRUE)
points(bei_obs, col = '#00000040')
dev.off()

# Plot the elevation surface.
pdf('figures/beieffortelev.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$elev, main = 'Elevation', riblab = 'Meters', ribsep = 0.05)
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the gradient surface.
pdf('figures/beieffortgrad.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei.extra$grad, main = 'Gradient', rib.sep = 0.05)
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
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
bei_spde <- inla.spde2.matern(mesh = bei_mesh)

# Set up a projection from the SPDE representation to a 400x200 grid.
bei_proj <- inla.mesh.projector(bei_mesh, dims = c(400, 200))


# Plot the mesh and point pattern.
pdf('figures/beieffortmesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 0))
plot(bei_mesh, asp = 1, main = '')
plot(bei_xsect, add = TRUE, col = '#ff000080', border = NA)
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

# Plot the piecewise linear approximation of the elevation surface.
pdf('figures/beieffortelevmesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_elev)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        riblab = 'Meters', ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Elevation')
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the piecewise linear approximation of the gradient surface.
pdf('figures/beieffortgradmesh.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_mesh_grad)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        ribsep = 0.05,
        main = 'Piecewise Linear Approximation of Gradient')
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
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

# Take the intersection of obs_psp with each triangle.
mesh_tris <- apply(bei_mesh$graph$tv, 1, function(x){
  return(owin(poly = meshloc[x,]))
})
tri_areas <- sapply(mesh_tris, area)
lines_subsegs <- lapply(mesh_tris, function(x){return(bei_psp[x])})

# Calculate the area represented by each node.
mesh_area <- rep(0, bei_mesh$n)
for(i in seq_along(tri_areas)){
  mesh_area[bei_mesh$graph$tv[i,]] <- mesh_area[bei_mesh$graph$tv[i,]] +
    tri_areas[i] / 3
}

# Track which triangle each segment came from.
seg_tri_idx <- unlist(lapply(seq_along(lines_subsegs), function(i){return(rep(i, lines_subsegs[[i]]$n))}))

# Get the midpoints
lines_midpoints <- as.matrix(do.call(rbind, lapply(lapply(lines_subsegs, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
lines_bary <- inla.mesh.projector(bei_mesh, loc = lines_midpoints)$proj$A

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
bei_int_weights <- diag(inla.mesh.fem(bei_mesh)$c0)

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
bei_pseudodata_exp <- c(bei_int_weights * mesh_weights, rep(0, bei_n_events))

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
pdf('figures/beieffortmean.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(E(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Predicted Mean of Latent GP')
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior standard deviation of the latent surface.
pdf('figures/beieffortsd.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(inla.mesh.project(bei_proj, bei_result$summary.random$idx$sd)),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
     riblab = expression(SD(bold(e)(u)*'|'*bold(x))), ribsep = 0.05,
     main = 'Posterior Prediction SD of Latent GP')
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the posterior mean of the linear predictor.
pdf('figures/beieffortbetas.pdf', width = 12, height = 6)
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
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot the backtransformed posterior mean of the intensity surface.
pdf('figures/beieffortintensity.pdf', width = 12, height = 6)
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
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Plot posterior marginals of the covariance parameters and intercept.
pdf('figures/beieffortpost.pdf', width = 18, height = 6)
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
pdf('figures/beieffortpostcoefs.pdf', width = 12, height = 6)
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
  bei_resid_df$area[r] <- area(bei_win[owin(c(bei_resid_df$x[r] - dx,
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
pdf('figures/beieffortresiduals.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(matrix(bei_resid_df$pearson, nrow = length(unique(bei_resid_df$x)))),
        unique(bei_resid_df$x), unique(bei_resid_df$y),
        unitname = c('meter', 'meters')),
     ribseb = 0.05, main = 'Gridded Pearson Residuals')
plot(bei_xsect, border = NA, col = '#ffffff80', add = TRUE)
points(bei_obs, pch = 21, cex = 0.5, col = '#00000080', bg = '#ffffff80')
dev.off()

# Set up a projection from the SPDE representation to the event locations.
bei_event_proj <- inla.mesh.projector(bei_mesh,
  loc = as.matrix(as.data.frame(bei_obs)))

# Create a copy of bei and mark with 1/sqrt(lambda).
bei_marked <- bei_obs
marks(bei_marked) <- as.vector(1/sqrt(exp(
  bei_result$summary.fixed['intercept', 'mean'] +
    bei_result$summary.fixed['elev', 'mean'] *
  inla.mesh.project(bei_event_proj, bei_mesh_elev) +
    bei_result$summary.fixed['grad', 'mean'] *
  inla.mesh.project(bei_event_proj, bei_mesh_grad) +
  inla.mesh.project(bei_event_proj, bei_result$summary.random$idx$mean)
)))

# Mark plot.
pdf('figures/beieffortmarkplot.pdf', width = 12, height = 6)
par(mar = c(0, 0, 2, 2))
plot(im(t(sqrt(exp(
          bei_result$summary.fixed['intercept', 'mean'] +
          bei_result$summary.fixed['elev', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_elev) +
          bei_result$summary.fixed['grad', 'mean'] *
            inla.mesh.project(bei_proj, bei_mesh_grad) +
          inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)
        ))),
        xrange = bei_win$x,
        yrange = bei_win$y,
        unitname = c('meter', 'meters')),
        ribsep = 0.05, main = 'Mark Plot')
plot(bei_win, border = 'white', add = TRUE)
plot(bei_marked, col = '#ffffff80', add = TRUE)
dev.off()

# Lurking variable plots.
lambda_im <- im(t(exp(
  bei_result$summary.fixed['intercept', 'mean'] +
  bei_result$summary.fixed['elev', 'mean'] *
    inla.mesh.project(bei_proj, bei_mesh_elev) +
  bei_result$summary.fixed['grad', 'mean'] *
    inla.mesh.project(bei_proj, bei_mesh_grad) +
  inla.mesh.project(bei_proj, bei_result$summary.random$idx$mean)
  )),
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
    sub_win <- bei_xsect[x_im <= x]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_xsect)
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
    return(data.frame(x, observed, expected, pearson, a, v, upper, lower))
  })
))

pdf('figures/beieffortlurkx.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ x, data = cum_pearson_x, type = 'l',
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Horizontal Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ x, data = cum_pearson_x, type = 'l', lty = 3)
lines(upper ~ x, data = cum_pearson_x, type = 'l', lty = 3)
dev.off()

cum_pearson_y <- do.call(rbind, c(
  list(data.frame(y = 0, observed = 0, expected = 0, pearson = 0, a = 0, v = 0, upper = 0, lower = 0)),
  lapply(seq(0, 500, 5), function(y){
    sub_win <- bei_xsect[y_im <= y]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_xsect)
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
    return(data.frame(y, observed, expected, pearson, a, v, upper, lower))
  })
))

pdf('figures/beieffortlurky.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ y, data = cum_pearson_y, type = 'l',
     xlab = 'Horizontal Coordinate', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Vertical Coordinate')
abline(h = 0, lty = 2)
lines(lower ~ y, data = cum_pearson_y, type = 'l', lty = 3)
lines(upper ~ y, data = cum_pearson_y, type = 'l', lty = 3)
dev.off()

cum_pearson_elev <- do.call(rbind,
  lapply(seq(118, 160, 2), function(elev){
    sub_win <- bei_xsect[bei.extra$elev <= elev]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_xsect)
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

pdf('figures/beieffortlurkelev.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ elev, data = cum_pearson_elev, type = 'l',
     xlab = 'Elevation', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Elevation')
abline(h = 0, lty = 2)
lines(lower ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
lines(upper ~ elev, data = cum_pearson_elev, type = 'l', lty = 3)
dev.off()

cum_pearson_grad <- do.call(rbind,
  lapply(seq(0, 0.35, 0.01), function(grad){
    sub_win <- bei_xsect[bei.extra$grad <= grad]
    sub_ppp <- bei_obs[sub_win]
    a <- area(sub_win)
    v <- a / area(bei_xsect)
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

pdf('figures/beieffortlurkgrad.pdf', width = 7, height = 6)
par(bty = 'n')
plot(pearson ~ grad, data = cum_pearson_grad, type = 'l',
     xlab = 'Gradient', ylab = 'Cumulative Pearson Residual',
     main = 'Lurking Variable Plot for Gradient')
abline(h = 0, lty = 2)
lines(lower ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
lines(upper ~ grad, data = cum_pearson_grad, type = 'l', lty = 3)
dev.off()
