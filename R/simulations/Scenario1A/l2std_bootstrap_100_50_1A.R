# Simulation code for Scenario 1A:
#   * Standardized L^2 statistic
#   * Bootstrap method
#   * n_1 = 100 (calibration sample size)
#   * n_2 = 50 (monitoring sample size)

## Object for saving simulation result (out of control signal)
senal_l2std_boot_100_50 <- vector("list", length(deltas))

## Object for saving simulation result (power)
potencia_l2std_boot_100_50 <- numeric(length(deltas))

## Initialize seeds for reproducibility
seeds <- 345 + seq_along(deltas)

## Calibration sample size
n1 <- 100

## Monitoring sample size
n2 <- 50

for (i in seq_along(deltas)) {
  delta <- deltas[i]
  start <- Sys.time()
  cat("Running simulation for delta =", delta, "\n")

  senal <- foreach(
    g = seq_len(mc),
    .packages = c("fda", "fda.usc", "fdahotelling", "qcr"),
    .options.RNG = seeds[i]
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
          depth = depth.mode,
          nb = 2000,
          plot = FALSE,
          ns = 0.05
        )

        if (length(fddep$out) > 0) {
          calibrado <- calibrado[-fddep$out, , drop = FALSE]
        }
      }

      n1_boot <- nrow(calibrado)
      n2_boot <- n2
      n_boot <- n1_boot + n2_boot

      boot_estadistico <- rep(NA_real_, B)

      mu <- rep(0, ncol(calibrado))

      for (b in seq_len(B)) {
        ind <- sample(
          seq_len(n1_boot),
          size = n_boot,
          replace = TRUE
        )

        data_boot <- calibrado[ind, , drop = FALSE]

        cov.est <- var(data_boot)

        eig <- eigen(cov.est, symmetric = TRUE)

        lambda <- eig$values
        V <- eig$vectors

        if (!all(lambda >= -tol * abs(lambda[1]))) {
          next
        }

        L.boot <- V %*%
          diag(sqrt(pmax(lambda * smo, 0)))

        dat <- data_boot +
          t(
            L.boot %*%
              matrix(
                rnorm(n_boot * ncol(data_boot)),
                nrow = ncol(data_boot)
              )
          )

        boot_1 <- dat[
          1:n1_boot,
          ,
          drop = FALSE
        ]

        boot_2 <- dat[
          (n1_boot + 1):n_boot,
          ,
          drop = FALSE
        ]

        boot_estadistico[b] <-
          fdahotelling:::stat_L2_std(
            x = boot_1,
            y = boot_2
          )
      }

      UCL <- quantile(
        boot_estadistico,
        probs = 1 - alpha,
        na.rm = TRUE
      )

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

          matriz_grupo <- t(f0(n2_boot))
        } else {
          f1 <- func.sim.set(
            t = tt,
            var.teor = var.teor,
            trend.teor = mu0 + delta,
            corr.teor = corr.teor,
            rho = rho
          )

          matriz_grupo <- t(f1(n2_boot))
        }

        estad_monitor[k] <-
          fdahotelling:::stat_L2_std(
            x = matriz_grupo,
            y = calibrado
          )
      }

      senal_grafico <- estad_monitor > UCL

      list(
        s = senal_grafico,
        estad_monitor = estad_monitor,
        UCL = UCL,
        boot = boot_estadistico
      )
    }

  senal_l2std_boot_100_50[[i]] <- unlist(lapply(senal, function(x) x$s))
  potencia_l2std_boot_100_50[i] <- mean(senal_l2std_boot_100_50[[i]])

  cat("\t", format(Sys.time() - start, digits = 3), "\n")
}

save(
  potencia_l2std_boot_100_50,
  senal_l2std_boot_100_50,
  file = "results/simulations/l2std_boot_100_50_1A.RData"
)
