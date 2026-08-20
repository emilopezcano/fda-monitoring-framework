## Setup for simulations

## Parallel processing setup
if (!exists("cl") || !inherits(cl, "cluster")) {
  ncores <- parallel::detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

## Functional mean constant shift
deltas <- c(0, 0.3, 0.5, 0.7, 0.9)

## Smoothing parameter $\gamma$
smo <- 0.05

## False alarm rate $\alpha$
alpha <- 0.05

## Number of Monte Carlo Simulations
mc <- mc_reps <- 1000
mc_chart <- 500

## Number of bootstrap samples
B <- 500

## Number of permutations
P <- 500

## Number of subgroups for monitoring
K <- 20

## Tolerance for eigenvalues
tol <- 1e-6

## Grid of time points and means for functional data
tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3 / 2)

## Theoretical variance 
var.teor <- 1

## Theoretical correlation structure (Scenarios A)
corr.teor.A <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2)
)

## Theoretical correlation structure (Scenarios B)
corr.teor.B <- outer(
  tt,
  tt,
  function(s, t)
    exp(-2 * (s - t)^2) * (s + 0.5) * (t + 0.5)
)

## Shape change parameters scenarios 2 (L1, L2)
etas.L <- c(0, 0.3, 0.5, 0.6)

## Shape change parameters scenarios 2 (T2)
etas.T <- c(0, 0.02, 0.03, 0.05)

## Trimming parameter
rho <- 0
