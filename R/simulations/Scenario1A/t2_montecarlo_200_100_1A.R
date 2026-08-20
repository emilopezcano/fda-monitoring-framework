# Simulation code for Scenario 1A:
#   * Hotelling T^2 statistic
#   * Montecarlo method
#   * n_1 = 200 (calibration sample size)
#   * n_2 = 100 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.A

## Set seed for reproducibility
set.seed(123)

## Calibration sample size
n1 <- 200

## Monitoring sample size
n2 <- 100

## Phase I function generator
f0 <- func.sim.set(
  tt,
  var.teor,
  trend.teor = mu0,
  corr.teor,
  rho
)

## Object for saving simulation result (power)
potencia_t2_montecarlo_200_100 <- numeric(length(deltas))

## Object for saving simulation result (out of control signal)
senal_t2_montecarlo_200_100 <- vector("list", length(deltas))

for (i in seq_along(deltas)) {
  delta <- deltas[i]
  start <- Sys.time()
  cat("Running simulation for delta =", delta, "\n")

  senal_delta <- vector("list", mc_chart)

  f1 <- func.sim.set(
    t = tt,
    var.teor = var.teor,
    trend.teor = mu0 + delta,
    corr.teor = corr.teor,
    rho = rho
  )

  for (g in seq_len(mc_chart)) {
    # 1. UCL Estimation
    estadistico_h0 <- numeric(mc_reps)

    calibrado_h0 <- t(f0(n1))

    mu <- rep(0, ncol(calibrado_h0))

    for (b in seq_len(mc_reps)) {
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_hotelling_impl(
          x = calibrado_h0,
          y = matriz_grupo_mc,
          mu = mu,
          paired = FALSE,
          step_size = 0.02,
          use_correction = FALSE,
          tolerance = tol
        )
    }
    # UCL

    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )
    # 2. Monitoring Phase II

    senal_grafico_t2_montecarlo_200_100 <- logical(K)

    for (k in seq_len(K)) {
      if (delta == 0) {
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

      senal_grafico_t2_montecarlo_200_100[k] <-
        T2 > UCL
    }

    senal_delta[[g]] <-
      senal_grafico_t2_montecarlo_200_100
  }

  # Resultado final para delta

  senal_t2_montecarlo_200_100[[i]] <-
    unlist(senal_delta)

  potencia_t2_montecarlo_200_100[i] <-
    mean(senal_t2_montecarlo_200_100[[i]])

  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

# Save

save(
  potencia_t2_montecarlo_200_100,
  senal_t2_montecarlo_200_100,
  file = "results/simulations/t2_montecarlo_200_100_1A.RData"
)
