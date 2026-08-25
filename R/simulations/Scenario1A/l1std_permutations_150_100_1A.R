# Simulation code for Scenario 1A:
#   * Standardized L^1 statistic
#   * Permutations method
#   * n_1 = 150 (calibration sample size)
#   * n_2 = 100 (monitoring sample size)

## Theoretical correlation structure (Scenarios A)
corr.teor <- corr.teor.A

## Object for saving simulation result (power)
potencia_l1std_perms_150_100 <- numeric(length(deltas))

## Object for saving simulation result (out of control signal)
senal_l1std_perms_150_100 <- vector("list", length(deltas))

## Initialize seeds for reproducibility
seeds_perm <- 1000 + seq_along(deltas)

## Calibration sample size
n1 <- 150

## Monitoring sample size
n2 <- 100




cat("--- Scenario 1A simulation for L1, Permutations,", n1, "/", n2, "\n")
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
        trend.teor = mu0 + delta,
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

      if (delta != 0) {
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
        if (delta == 0) {
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
            trend.teor = mu0 + delta,
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

  senal_l1std_perms_150_100[[i]] <- unlist(lapply(senal, function(x) x$s))

  potencia_l1std_perms_150_100[i] <- mean(senal_l1std_perms_150_100[[i]])
  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}


save(
  potencia_l1std_perms_150_100,
  senal_l1std_perms_150_100,
  file = "results/simulations/l1std_perms_150_100_1A.RData"
)

end <- Sys.time()
cat(
  "[",
  format(end, "%HH:%MM"),
  "] END simulation Scenario 1A for L1, Permutations,", n1, "/", n2, "\n",
  format(end - start0, digits = 3),
  sep = ""
)
