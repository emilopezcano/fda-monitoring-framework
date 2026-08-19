# Simulation code for Scenario 1A:
#   * Standardized L^1 statistic
#   * Montecdarlo method
#   * n_1 = 100 (calibration sample size)
#   * n_2 = 50 (monitoring sample size)

## Set seed for reproducibility
set.seed(123)

## Calibration sample size
n1 <- 100

## Monitoring sample size
n2 <- 50

## Phase I function generator
f0 <- func.sim.set(
  tt,
  var.teor,
  trend.teor = mu0,
  corr.teor,
  rho
)
## Object for saving simulation result (out of control signal)
senal_l1std_montecarlo_100_50 <- vector("list", length(deltas))

## Object for saving simulation result (power)
potencia_l1std_montecarlo_100_50 <- numeric(length(deltas))

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

    for (b in seq_len(mc_reps)) {
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo_mc
        )
    }

    # UCL específico del gráfico
    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    # 2. Monitoring Phase II
    senal_grafico_l1std_montecarlo_100_50 <- logical(K)

    for (k in seq_len(K)) {
      if (delta == 0) {
        matriz_grupo <- t(f0(n2))
      } else {
        matriz_grupo <- t(f1(n2))
      }

      l1std <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo
        )

      senal_grafico_l1std_montecarlo_100_50[k] <-
        l1std > UCL
    }

    # Result for the chart
    senal_delta[[g]] <-
      senal_grafico_l1std_montecarlo_100_50
  }

  # Final result for delta
  senal_l1std_montecarlo_100_50[[i]] <-
    unlist(senal_delta)

  potencia_l1std_montecarlo_100_50[i] <-
    mean(senal_l1std_montecarlo_100_50[[i]])
  
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

## Save ----
save(
  potencia_l1std_montecarlo_100_50,
  senal_l1std_montecarlo_100_50,
  file = "results/simulations/l1std_montecarlo_100_50.RData"
)
