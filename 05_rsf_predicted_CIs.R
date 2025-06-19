library(MASS)

coefs <- readRDS("output/pigs_vcov.rds")

rasters <- cty_layers
#cov_stack is stacked rasters
# Create the squared covariate
x_sq <- rasters[["dist_water"]]^2
names(x_sq) <- "water2"  # name it to match the model formula

# Add it to the covariate stack
rasters <- c(rasters, x_sq)


rsf_CIs <- function(coefs, rasters) {
  
  # Extract coefficients and their variance-covariance matrix
  beta_hat <- coefs[[1]]
  vcov_mat <- coefs[[2]]
  
  # Prediction matrix (X_pred)
  # Get the raster values as a matrix (rows = cells, cols = covariates)
  X_pred <- values(rasters, mat = TRUE)

  # 2.Simulate RSF predictions
  set.seed(42)
  n_sims <- 1000  # number of draws
  beta_sim <- mvrnorm(n_sims, mu = beta_hat, Sigma = vcov_mat)  # matrix: n_sims x n_covs
  
  # Initialize matrix to hold predicted RSF values
  rsf_matrix <- matrix(NA_real_, nrow = nrow(X_pred), ncol = n_sims)
  
  # Compute RSF predictions: w(x) = exp(X %*% beta)
  eta_matrix <- X_pred %*% t(beta_sim)  # Matrix: pixels x sims
  
  # 3. Pixel-wise confidence intervals
  eta_lower <- apply(eta_matrix, 1, quantile, probs = 0.025, na.rm=TRUE)
  eta_upper <- apply(eta_matrix, 1, quantile, probs = 0.975, na.rm=TRUE)
  
  # Exponentiate
  rsf_lower <- exp(eta_lower)
  rsf_upper <- exp(eta_upper)

  # 4. Write back to raster layers
  # Use any covariate layer to copy geometry
  template <- rasters[[1]]
  
  r_lower <- template; values(r_lower) <- rsf_lower
  r_upper <- template; values(r_upper) <- rsf_upper

  # Log-transform
  r_lower <- r_lower + 1E-06
  r_upper <- r_upper + 1E-06
  r_lower <- log(r_lower)
  r_upper <- log(r_upper)
  
  # Combine into a multi-layer raster
  rsf_ci_stack <- c(r_lower, r_upper)
  names(rsf_ci_stack) <- c("rsf_lower", "rsf_upper")
  
  # Define the output file names. 
  filenames <- paste0("output/deer_county_preds/deer_pred_", split_counties[[i]]$COUNTYNAME,"_", names(rsf_ci_stack))
  
  # 5. Write each layer to a separate file
  writeRaster(rsf_ci_stack, filenames, overwrite=TRUE)
  
}

dim(X_pred)       # Should be (n_pixels, n_covariates)
dim(vcov_mat)     # Should be (n_covariates, n_covariates)


# Compute SE of linear predictor at each pixel
# Variance of each row: diag(X %*% V %*% t(X))
eta_var <- rowSums((X_pred %*% vcov_mat * X_pred))  # efficient shortcut
eta_se <- sqrt(eta_var)

template <- rasters[[1]]

r_se <- template; values(r_se) <- eta_se

r_se <- log(r_se)
plot(r_se)



apply(X_pred, 2, summary)

summary(diag(vcov_mat))  # Reasonable? Or are some >100?
