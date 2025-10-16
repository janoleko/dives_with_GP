library(splines2)
library(RTMB)
library(fmesher)
library(MASS)

# ============================================================
#  Gaussian process mean model based on a varying-coefficient ODE
#  (Depth profile model for one dive)
#
#  ODE formulation:
#    dy/dt = s(T) * B(u) %*% beta      where  u = t / T,  t ∈ [0, T]
#
#  Solution:
#    y(t) = y0 + s(T) * ∫_0^{t} B(s/T) %*% beta ds
#         = y0 + s(T) * T * ∫_0^{u} B(r) %*% beta dr
#         = y0 + s(T) * A(u) %*% beta
#
#  where:
#    - B(u) : B-spline basis functions (velocity profile over relative time)
#    - A(u) : integrated B-spline basis (depth profile over relative time)
#    - beta : coefficients defining the shape of the velocity curve
#    - s(T) : scaling function controlling how depth scales with duration
# ============================================================

# ---------------------------
# Dive setup
# ---------------------------
T_dive <- 200  # dive duration in seconds (e.g., 10 minutes at 1 Hz sampling)
u <- seq(0, 1, length.out = T_dive)  # relative time through dive in [0,1]

# ---------------------------
# Basis setup
# ---------------------------
df <- 10      # number of B-spline basis functions
degree <- 3  # spline degree (cubic)

# Velocity basis: B(u)
B <- bSpline(u, df = df, degree = degree, intercept = TRUE)

# Integrated basis: A(u) = ∫_0^u B(s) ds
# 'integral = TRUE' makes bSpline return the integrated basis directly
A <- bSpline(u, df = df, degree = degree, intercept = TRUE, integral = TRUE)

# constraint at the end of the dive
A_end <- A[nrow(A), , drop=FALSE]

# A_end: vector of length df
N <- diag(df) - matrix(A_end, ncol=1) %*% (A_end / sum(A_end^2))

A <- A %*% N
B <- B %*% N

# ---------------------------
# Scaling function s(T)
# ---------------------------
# Controls how overall depth scales with dive duration.
# Examples:
#   s <- function(T) 1          # shape fixed on relative-time scale
#   s <- function(T) T          # linear scaling with duration (velocity model)
#   s <- function(T) T^gamma    # flexible scaling, estimate gamma
s <- function(T) T  # default: linear scaling with duration

# ---------------------------
# Simulate one dive
# ---------------------------
# Beta defines the shape of the velocity curve over relative time.
# Negative coefficients = descent; positive = ascent.
beta <- c(0, -4, -2, -3, -0.5, 0.5, 2, 0, 0.5, 2)
y0 <- 0  # starting depth (surface)

# Mean depth profile from the ODE solution
y_mean <- y0 + s(T_dive) * (A %*% beta)

# ---------------------------
# Plot depth and velocity
# ---------------------------
plot(u * T_dive, y_mean, type = "l", lwd = 2, col = "steelblue",
     xlab = "Time through dive",
     ylab = "Depth (arbitrary units)",
     main = "Mean dive profile from ODE with scaling s(T)")
abline(h = 0, col = "#00000020")

# Velocity (for intuition): v(u) = B(u) %*% beta
v <- B %*% beta
lines(u * T_dive, v, col = "grey50", lty = 2)
legend("bottomright", legend = c("Depth profile", "Velocity (relative)"),
       col = c("steelblue", "grey50"), lty = c(1, 2), bty = "n")



# Now sample GP with that mean --------------------------------------------

# creating mesh and finite element matrices
mesh <- fm_mesh_1d(u, boundary = "dirichlet")
spde <- fm_fem(mesh)

kappa <- 0.1
tau <- 0.01

Q <- tau^2*(kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2) # sparse precision

plot(u * T_dive, y_mean, type = "l", lwd = 1, col = "steelblue",
     ylim = c(min(y_mean) * 1.2, 0),
     xlab = "Time through dive",
     ylab = "Depth (arbitrary units)",
     main = "Mean dive profile from ODE with scaling s(T)")

set.seed(123)
for(i in 1:100) {
  field <- RTMB:::rgmrf0(1, Q)
  y_sim <- y_mean + c(0, field, 0)
  lines(u * T_dive, y_sim, col = "#00000010")
}

lines(u * T_dive, y_mean, type = "l", lwd = 3, col = "orange")

