# Simulation code for Scenario 1A:
#   * Standardized L^2 statistic
#   * Montecarlo method
#   * n_1 = 100 (calibration sample size)
#   * n_2 = 50 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.A

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

## Object for saving simulation result (power)
potencia_l2std_montecarlo_100_50 <- numeric(length(deltas))

## Object for saving simulation result (out of control signal)
senal_l2std_montecarlo_100_50 <- vector("list", length(deltas))




cat("--- Scenario 1A simulation for L2, Montecarlo,", n1, "/", n2, "\n")
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

    calibrado_h0 <- t(f0(n1))

    for (b in seq_len(mc_reps)) {
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_L2_std(
          x = calibrado_h0,
          y = matriz_grupo_mc
        )
    }
    # UCL

    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    # 2. Monitoring Phase II
    senal_grafico_l2std_montecarlo_100_50 <- logical(K)

    for (k in seq_len(K)) {
      if (delta == 0) {
        matriz_grupo <- t(f0(n2))
      } else {
        matriz_grupo <- t(f1(n2))
      }

      l2std <-
        fdahotelling:::stat_L2_std(
          x = calibrado_h0,
          y = matriz_grupo
        )

      senal_grafico_l2std_montecarlo_100_50[k] <-
        l2std > UCL
    }

    senal_delta[[g]] <-
      senal_grafico_l2std_montecarlo_100_50
  }

  senal_l2std_montecarlo_100_50[[i]] <-
    unlist(senal_delta)

  potencia_l2std_montecarlo_100_50[i] <-
    mean(senal_l2std_montecarlo_100_50[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

# Save

save(
  potencia_l2std_montecarlo_100_50,
  senal_l2std_montecarlo_100_50,
  file = "results/simulations/l2std_montecarlo_100_50_1A.RData"
)

end <- Sys.time()
cat(
  "[",
  format(end, "%HH:%MM"),
  "] END simulation Scenario 1A for L2, Montecarlo,", n1, "/", n2, "\n",
  format(end - start0, digits = 3),
  sep = ""
)
