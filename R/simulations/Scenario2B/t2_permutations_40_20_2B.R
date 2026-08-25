# Simulation code for Scenario 2B:
#   * Standardized T2 statistic
#   * Permutation method
#   * n_1 = 40 (calibration sample size)
#   * n_2 = 20 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.B

## Values for etas according to the statistic used (L1, L2, T2)
etas <- etas.T

## Object for saving simulation result (power)
potencia_t2_perms_40_20 <- numeric(length(etas))

## Object for saving simulation result (out of control signal)
senal_t2_perms_40_20 <- vector("list", length(etas))

## Initialize seeds for reproducibility
seeds_perm <- 1000 + seq_along(etas)

## Calibration sample size
n1 <- 40

## Monitoring sample size
n2 <- 20



cat("--- Scenario 2B simulation for T2, Permutations,", n1, "/", n2, "\n")
start0 <- Sys.time()

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

  senal <- foreach(
    g = seq_len(mc),
    .packages = c("fda", "fda.usc", "fdahotelling", "qcr"),
    .options.RNG = seeds_perm[i]
  ) %dorng%
    {
      f0 <- func.sim.set(
        t = tt,
        var.teor = var.teor,
        trend.teor = mu0,
        corr.teor = corr.teor,
        rho = rho
      )

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

      eps <- 0.05

      n_out <- max(1, ceiling(eps * n1))

      calibrado <- rbind(
        t(f0(n1 - n_out)),
        t(f1(n_out))
      )

      calibrado <- calibrado[
        sample(seq_len(nrow(calibrado))),
        ,
        drop = FALSE
      ]

      if (eta != 0) {
        fdchart <- fdqcd(as.data.frame(calibrado))

        fddep <- fdqcs.depth(
          fdchart,
          depth = fda.usc::depth.mode,
          nb = 2000,
          plot = FALSE,
          ns = 0.05
        )

        if (length(fddep$out) > 0) {
          calibrado <- calibrado[-fddep$out, , drop = FALSE]
        }
      }

      n_cal <- nrow(calibrado)

      n1_perm <- n_cal
      n2_perm <- n2
      n_total <- n1_perm + n2_perm

      perm_estadistico <- numeric(P)

      mu <- rep(0, ncol(calibrado))

      ind <- sample(
        seq_len(n_cal),
        size = n_total,
        replace = TRUE
      )

      data_perm <- calibrado[ind, , drop = FALSE]

      ##########################################################################
      # Permutations
      ##########################################################################

      for (p in seq_len(P)) {
        perm <- sample(n_total)

        D1 <- data_perm[
          perm[1:n1_perm],
          ,
          drop = FALSE
        ]

        D2 <- data_perm[
          perm[(n1_perm + 1):n_total],
          ,
          drop = FALSE
        ]

        perm_estadistico[p] <-
          fdahotelling:::stat_hotelling_impl(
            x = D1,
            y = D2,
            mu = mu,
            paired = FALSE,
            step_size = 0.02,
            use_correction = FALSE,
            tolerance = tol
          )
      }

      UCL <- quantile(
        perm_estadistico,
        probs = 1 - alpha,
        na.rm = TRUE
      )

      estad_monitor <- numeric(K)

      for (k in seq_len(K)) {
        if (eta == 0) {
          f0 <- func.sim.set(
            t = tt,
            var.teor = var.teor,
            trend.teor = mu0,
            corr.teor = corr.teor,
            rho = rho
          )

          matriz_grupo <- t(f0(n2))
        } else {
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

          matriz_grupo <- t(f1(n2))
        }

        estad_monitor[k] <-
          fdahotelling:::stat_hotelling_impl(
            x = calibrado,
            y = matriz_grupo,
            mu = mu,
            paired = FALSE,
            step_size = 0.02,
            use_correction = FALSE,
            tolerance = tol
          )
      }

      senal_grafico <- estad_monitor > UCL

      list(
        s = senal_grafico,
        estad_monitor = estad_monitor,
        UCL = UCL,
        per = perm_estadistico
      )
    }

  senal_t2_perms_40_20[[i]] <- unlist(lapply(senal, function(x) x$s))

  potencia_t2_perms_40_20[i] <- mean(senal_t2_perms_40_20[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}



save(
  potencia_t2_perms_40_20,
  senal_t2_perms_40_20,
  file = "results/simulations/t2_perms_40_20_2B.RData"
)

end <- Sys.time()
cat(
  "[",
  format(end, "%HH:%MM"),
  "] END simulation Scenario 2B for T2, Permutations,", n1, "/", n2, "\n",
  format(end - start0, digits = 3), "\n",
  sep = ""
)
