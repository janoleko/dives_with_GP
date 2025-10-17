library(LaMa)
library(splines2)
library(RTMB)
library(fmesher)

# very wasteful right now
sim_profile <- function(beta, sigma, rho, T_dive = 100, s = function(T) T_dive, y0 = 0) {
  u <- seq(0, 1, length.out = T_dive)  # relative time through dive in [0,1]

  df <- length(beta)
  degree <- 3  # spline degree (cubic)

  # Velocity basis: B(u)
  B <- bSpline(u, df = df, degree = degree, intercept = TRUE)
  # Integrated basis: A(u) = ∫_0^u B(s) ds
  A <- bSpline(u, df = df, degree = degree, intercept = TRUE, integral = TRUE)

  # we need the constraint y(1) = 0, i.e. A(1) %*% beta = 0
  A_end <- A[nrow(A), , drop=FALSE]
  N <- diag(df) - matrix(A_end, ncol=1) %*% (A_end / sum(A_end^2)) # maps beta into nullspace of A_end
  A <- A %*% N # design matrices w constraint that y(T) = 0
  B <- B %*% N

  # Mean depth profile from the ODE solution
  y_mean <- y0 + s(T_dive) * (A %*% beta)

  # creating mesh and finite element matrices
  mesh <- fm_mesh_1d(u * T_dive, boundary = "dirichlet")
  spde <- fm_fem(mesh)

  # transform marginal sd and range to SPDE parameters
  kappa <- 1 / rho
  tau   <- 1 / (2 * sqrt(pi * kappa) * sigma)

  Q <- tau^2*(kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2) # sparse precision
  field <- RTMB:::rgmrf0(1, Q)
  y_sim <- y_mean + c(0, field, 0) # boundary conditions -> fixed at 0

  list(
    depth = as.numeric(y_sim),
    mean = as.numeric(y_mean)
  )
}

sim_depth_hmm <- function(beta, sigma, rho,
                          n_dives = 20, T_dive = 100,
                          n_states = 2,
                          Gamma = tpm(rep(qlogis(0.2), 2))
                          ) {
  delta <- stationary(Gamma)

  # sim states
  s <- rep(NA, n_dives)
  s[1] <- sample(1:n_states, 1, prob = delta)
  for(d in 2:n_dives) {
    s[d] <- sample(1:n_states, 1, prob = Gamma[s[d-1], ])
  }

  # sample profiles conditional on dives
  depth <- means <- list()
  for(d in 1:n_dives) {
    sim <- sim_profile(T_dive = T_dive,
                       beta = beta[s[d], ],
                       sigma = sigma[s[d]],
                       rho = rho[s[d]])
    depth[[d]] <- sim$depth
    means[[d]] <- sim$mean
  }

  list(
    depth = depth,
    mean = means,
    state = s
  )
}

beta <- matrix(
  c(-2, -0.2, -0.2, -1, 0, 0, 0.1, 0.1, 0.5, 3,
    -6, -4, -2, -3, -0.5, 0.5, 2, 0, 0.5, 4),
  nrow = 2, byrow = TRUE)

sigma <- c(0.1, 0.6)         # marginal standard deviation in depth units
rho <- c(10, 10)

Gamma <- tpm(c(qlogis(0.2), qlogis(0.1)))


data <- sim_depth_hmm(beta, sigma, rho)

depth <- simplify2array(data$depth)
depth <- as.numeric(depth)
plot(depth, type = "l")

mean <- simplify2array(data$mean)
mean <- as.numeric(mean)
plot(mean, type = "l")
