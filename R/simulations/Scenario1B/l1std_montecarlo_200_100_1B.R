# Simulation code for Scenario 1B:
#   * Standardized L^1 statistic
#   * Montecarlo method
#   * n_1 = 200 (calibration sample size)
#   * n_2 = 100 (monitoring sample size)

## Theoretical correlation structure (Scenarios B)
corr.teor <- corr.teor.B

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
potencia_l1std_montecarlo_200_100 <- numeric(length(deltas))

## Object for saving simulation result (out of control signal)
senal_l1std_montecarlo_200_100 <- vector("list", length(deltas))





cat("--- Scenario 1B simulation for L1, Montecarlo,", n1, "/", n2, "\n")
start0 <- Sys.time()

for (i in seq_along(deltas)) {
  delta <- deltas[i]
  start <- Sys.time()
  cat(
    "[",
    format(start, "%HH:%MM"),
    "] Running simulation for delta = ",
    delta,
    "\n",
    sep = ""
  )


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

    # Phase I sample under H0, fixed for this control chart
    calibrado_h0 <- t(f0(n1))

    for (b in seq_len(mc_reps)) {
      # Generate a Phase II sample under H0
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo_mc
        )
    }

    # Chart-specific UCL
    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    # 2. Phase II Monitoring
    senal_grafico_l1std_montecarlo_200_100 <- logical(K)

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

      senal_grafico_l1std_montecarlo_200_100[k] <-
        l1std > UCL
    }

    # Control chart result
    senal_delta[[g]] <-
      senal_grafico_l1std_montecarlo_200_100
  }

  # Final result for delta
  senal_l1std_montecarlo_200_100[[i]] <-
    unlist(senal_delta)

  potencia_l1std_montecarlo_200_100[i] <-
    mean(senal_l1std_montecarlo_200_100[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

## Save ----
save(
  potencia_l1std_montecarlo_200_100,
  senal_l1std_montecarlo_200_100,
  file = "results/simulations/l1std_montecarlo_200_100_1B.RData"
)

end <- Sys.time()
cat(
  "[",
  format(end, "%HH:%MM"),
  "] END simulation Scenario 1B for L1, Montecarlo,", n1, "/", n2, "\n",
  format(end - start0, digits = 3),
  sep = ""
)
