# Simulation code for Scenario 1A: 
#   * Functional mean constant shift with deltas 0, 0.3, 0.5, 0.7 and 0.9
#   * Bootstrap method
#   * n_1 = 100 (calibration sample size)
#   * n_2 = 50 (monitoring sample size)
#   * Standardized L^1 statistic

## Parallel processing setup

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoSNOW(cl)

## Functional mean constant shift
deltas <- c(0, 0.3, 0.5, 0.7, 0.9)

## Object for saving simulation result
senal_l1std_boot_100_50 <- vector("list", length(deltas))

## Smoothing parameter $\gamma$
smo <- 0.05

## Object for saving simulation result (power)
potencia_l1std_boot_100_50 <- numeric(length(deltas))

## False alarm rate $\alpha$
alpha <- 0.05

## Random seeds for reproducibility
seeds <- 345 + seq_along(deltas)

## Number of Monte Carlo Simulations
mc <- 1000

## Number of bootstrap samples
B <- 500

## Number of subgroups for monitoring
K <- 20

## Tolerance for eigenvalues
tol <- 1e-6

## Calibration sample size
n1 <- 100

## Monitoring sample size
n2 <- 50

## Grid of time points and means for functional data
tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3 / 2)

## Theoretical variance and correlation structure for the functional data
var.teor <- 1
corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2)
)

## Trimming parameter
rho <- 0


for (i in seq_along(deltas)) {
  delta <- deltas[i]
  cat("Running simulation for delta =", delta, "\n")

  senal <- foreach(
    g = seq_len(mc),
    .packages = c("fda", "fda.usc", "fdahotelling", "qcr"),
    # .export = c(
    #   "fdqcd",
    #   "fdqcs.depth",
    #   "fdqcs.depth.default",
    #   "fdqcs.depth.fdqcd",
    #   "func.sim.set",
    #   "mu0"
    # ),
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
          fdahotelling:::stat_L1_std(
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
          fdahotelling:::stat_L1_std(
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

  senal_l1std_boot_100_50[[i]] <- unlist(lapply(senal, function(x) x$s))
  potencia_l1std_boot_100_50[i] <- mean(senal_l1std_boot_100_50[[i]])
}


stopCluster(cl)

save(
  potencia_l1std_boot_100_50,
  senal_l1std_boot_100_50,
  file = "l1std_boot_100_50_1A.RData"
)
