# Simulation code for Scenario 2A:
#   * Standardized T2 statistic
#   * Monte Carlo method
#   * n_1 = 40 (calibration sample size)
#   * n_2 = 20 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.A

## Values for etas according to the statistic used (L1, L2, T2)
etas <- etas.T

## Set seed for reproducibility
set.seed(123)

## Calibration sample size
n1 <- 40

## Monitoring sample size
n2 <- 20

## Phase I function generator
f0 <- func.sim.set(
  tt,
  var.teor,
  trend.teor = mu0,
  corr.teor,
  rho
)

## Object for saving simulation result (power)
potencia_t2_montecarlo_40_20 <- numeric(length(etas))

## Object for saving simulation result (out of control signal)
senal_t2_montecarlo_40_20 <- vector("list", length(etas))

cat("--- Scenario 2A simulation for T2, Montecarlo,", n1, n2, "\n")

for (i in seq_along(etas)) {
  start <- Sys.time()
  cat(
    "[",
    format(start, "%HH:%MM"),
    "] Running simulation for eta = ",
    eta,
    "\n",
    sep = ""
  )
  eta <- etas[i]
  senal_eta <- vector("list", mc_chart)
  f1 <- func.sim.set(
    t = tt,
    var.teor = var.teor,
    trend.teor = (1 - eta) *
      30 *
      tt *
      (1 - tt)^(3 / 2) +
      eta * 30 * tt^(3 / 2) * (1 - tt),
    corr.teor = corr.teor,
    rho = rho
  )

  for (g in seq_len(mc_chart)) {
    estadistico_h0 <- numeric(mc_reps)
    calibrado_h0 <- t(f0(n1))
    mu <- rep(0, ncol(calibrado_h0))

    for (b in seq_len(mc_reps)) {
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_hotelling_impl(
          x = as.matrix(calibrado_h0),
          y = as.matrix(matriz_grupo_mc),
          mu = mu,
          paired = FALSE,
          step_size = 0.02,
          use_correction = FALSE,
          tolerance = tol
        )
    }

    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    senal_grafico_t2_montecarlo_40_20 <- logical(K)

    for (k in seq_len(K)) {
      if (eta == 0) {
        matriz_grupo <- t(f0(n2))
      } else {
        matriz_grupo <- t(f1(n2))
      }

      T2 <-
        fdahotelling:::stat_hotelling_impl(
          x = calibrado_h0,
          y = matriz_grupo,
          mu = mu,
          paired = FALSE,
          step_size = 0.02,
          use_correction = FALSE,
          tolerance = tol
        )

      senal_grafico_t2_montecarlo_40_20[k] <- T2 > UCL
    }

    senal_eta[[g]] <- senal_grafico_t2_montecarlo_40_20
  }

  senal_t2_montecarlo_40_20[[i]] <- unlist(senal_eta)

  potencia_t2_montecarlo_40_20[i] <- mean(senal_t2_montecarlo_40_20[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}


end <- Sys.time()
cat(
    "[",
    format(end, "%HH:%MM"),
    "] END simulation Scenario 2A for T2, Montecarlo,", n1, n2, "\n",
    format(end - start, digits = 3),
    sep = ""
  )

save(
  potencia_t2_montecarlo_40_20,
  senal_t2_montecarlo_40_20,
  file = "results/simulations/t2_montecarlo_40_20_2A.RData"
)
