# Your Lefkovitch matrix
A <- matrix(c(
  # Example: replace with your actual matrix
  0.5, 0.7, 1.2,
  0.3, 0.0, 0.0,
  0.0, 0.6, 0.8
), nrow = 3, byrow = TRUE)

# Initial population vector
N0 <- c(100, 50, 30)

# Target harvest per year
H_target <- 50

# Fecundity indices (e.g., row 1, columns 2 and 3)
fecundity_indices <- c(1, 2, 3)  # adjust to match your matrix

# Stage from which harvest is taken (e.g., adults)
harvest_stage <- 3

# Carrying capacity for fecundity
Kf <- 500

simulate_theta <- function(theta, A, N0, harvest_stage, fecundity_indices, Kf, H_target, years = 50) {
  N <- N0
  total_harvest <- 0
  
  for (t in 1:years) {
    Nf_t <- sum(N[fecundity_indices])
    fecundity_modifier <- 1 / (1 + (Nf_t / Kf)^theta)
    
    A_mod <- A
    A_mod[fecundity_indices] <- A[fecundity_indices] * fecundity_modifier
    
    N <- A_mod %*% N
    h <- min(N[harvest_stage], H_target)
    N[harvest_stage] <- pmax(N[harvest_stage] - h, 0)
    total_harvest <- total_harvest + h
  }
  
  avg_harvest <- total_harvest / years
  return(abs(avg_harvest - H_target))  # minimize difference from target
}

optimize_theta <- optimize(
  f = simulate_theta,
  interval = c(0.01, 10),
  A = A,
  N0 = N0,
  harvest_stage = harvest_stage,
  fecundity_indices = fecundity_indices,
  Kf = Kf,
  H_target = H_target
)

# Best theta value
optimize_theta$minimum


theta_vals <- seq(0.1, 10, length.out = 100)
harvest_diff <- sapply(theta_vals, function(th) simulate_theta(th, A, N0, harvest_stage, fecundity_indices, Kf, H_target))

plot(theta_vals, harvest_diff, type = "l", col = "blue", lwd = 2,
     xlab = expression(theta), ylab = "Harvest Error", main = "Sensitivity of Harvest to Theta")
abline(h = 0, col = "red", lty = 2)





