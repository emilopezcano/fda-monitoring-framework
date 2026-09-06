# Simulation code for Scenario 2B:
#   * Standardized L1 statistic
#   * Permutation method
#   * n_1 = 60 (calibration sample size)
#   * n_2 = 40 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.B

## Values for etas according to the statistic used (L1, L2, T2)
etas <- etas.L

## Object for saving simulation result (power)
potencia_l1std_perms_60_40 <- numeric(length(etas))

## Object for saving simulation result (out of control signal)
senal_l1std_perms_60_40 <- vector("list", length(etas))

## Initialize seeds for reproducibility
seeds_perm <- 1000 + seq_along(etas)

## Calibration sample size
n1 <- 60

## Monitoring sample size
n2 <- 40



cat("--- Scenario 2B simulation for L1, Permutations,", n1, "/", n2, "\n")
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
          fdahotelling:::stat_L1_std(
            x = D1,
            y = D2
          )
      }

      UCL <- quantile(
        perm_estadistico,
        probs = 1 - alpha,
        na.rm = TRUE
      )

      ######################################################################
      # 5. Monitoreo
      ######################################################################

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
          fdahotelling:::stat_L1_std(
            x = calibrado,
            y = matriz_grupo
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

  senal_l1std_perms_60_40[[i]] <- unlist(lapply(senal, function(x) x$s))

  potencia_l1std_perms_60_40[i] <- mean(senal_l1std_perms_60_40[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}



save(
  potencia_l1std_perms_60_40,
  senal_l1std_perms_60_40,
  file = "results/simulations/l1std_perms_60_40_2B.RData"
)

end <- Sys.time()
cat(
  "[",
  format(end, "%HH:%MM"),
  "] END simulation Scenario 2B for L1, Permutations,", n1, "/", n2, "\n",
  format(end - start0, digits = 3), "\n",
  sep = ""
)
