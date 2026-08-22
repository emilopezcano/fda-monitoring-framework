# Simulation code for Scenario 2B:
#   * Standardized L1 statistic
#   * Monte Carlo method
#   * n_1 = 40 (calibration sample size)
#   * n_2 = 20 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.B

## Values for etas according to the statistic used (L1, L2, T2)
etas <- etas.L

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
potencia_l1std_montecarlo_40_20 <- numeric(length(etas))

## Object for saving simulation result (out of control signal)
senal_l1std_montecarlo_40_20 <- vector("list", length(etas))

cat("--- Scenario 2B simulation for L1, Montecarlo,", n1, n2, "\n")

for (i in seq_along(etas)) {
  eta <- etas[i]
  start <- Sys.time()
  cat(
    "[",
    format(start, "%HH:%MM"),
    "] Running simulation for eta = ",
    eta,
    "\n",
    sep = ""
  )

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

    for (b in seq_len(mc_reps)) {
      matriz_grupo_mc <- t(f0(n2))

      estadistico_h0[b] <-
        fdahotelling:::stat_L1_std(
          x = as.matrix(calibrado_h0),
          y = as.matrix(matriz_grupo_mc)
        )
    }

    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    senal_grafico_l1std_montecarlo_40_20 <- logical(K)

    for (k in seq_len(K)) {
      if (eta == 0) {
        matriz_grupo <- t(f0(n2))
      } else {
        matriz_grupo <- t(f1(n2))
      }

      L1std <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo
        )

      senal_grafico_l1std_montecarlo_40_20[k] <- L1std > UCL
    }

    senal_eta[[g]] <- senal_grafico_l1std_montecarlo_40_20
  }

  senal_l1std_montecarlo_40_20[[i]] <- unlist(senal_eta)

  potencia_l1std_montecarlo_40_20[i] <- mean(senal_l1std_montecarlo_40_20[[i]])

  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

end <- Sys.time()
cat(
    "[",
    format(end, "%HH:%MM"),
    "] END simulation Scenario 2B for L1, Montecarlo,", n1, n2, "\n",
    format(end - start, digits = 3),
    sep = ""
  )


save(
  potencia_l1std_montecarlo_40_20,
  senal_l1std_montecarlo_40_20,
  file = "results/simulations/l1std_montecarlo_40_20_2B.RData"
)
