library(MASS)

final <- readRDS("output/deer_vcov.rds")

#cov_stack is stacked rasters

# Extract coefficients and their variance-covariance matrix
beta_hat <- final[[1]]
vcov_mat <- final[[2]]

# Prediction matrix (X_pred)
# Get the raster values as a matrix (rows = cells, cols = covariates)
X_pred <- terra::values(cov_stack, mat = TRUE)
X_pred <- cbind(Intercept = 1, X_pred)

# 2.Simulate RSF predictions
set.seed(42)
n_sims <- 1000  # number of draws
beta_sim <- MASS::mvrnorm(n_sims, mu = beta_avg, Sigma = vcov_avg)  # matrix: n_sims x n_covs

# Initialize matrix to hold predicted RSF values
rsf_matrix <- matrix(NA_real_, nrow = nrow(X_pred), ncol = n_sims)

# Compute RSF predictions: w(x) = exp(X %*% beta)
for (i in 1:n_sims) {
  eta <- X_pred %*% beta_sim[i, ]
  rsf_matrix[, i] <- exp(eta)
}

# 3. Pixel-wise confidence intervals
rsf_lower <- apply(rsf_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
rsf_upper <- apply(rsf_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
rsf_mean  <- apply(rsf_matrix, 1, mean, na.rm = TRUE)

# 4. Write back to raster layers
# Use any covariate layer to copy geometry
template <- cov_stack[[1]]

r_lower <- template; values(r_lower) <- rsf_lower
r_upper <- template; values(r_upper) <- rsf_upper
r_mean  <- template; values(r_mean)  <- rsf_mean

# Combine into a multi-layer raster
rsf_ci_stack <- c(r_mean, r_lower, r_upper)
names(rsf_ci_stack) <- c("rsf_mean", "rsf_lower", "rsf_upper")


