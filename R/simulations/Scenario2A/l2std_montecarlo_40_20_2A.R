################################################################################
# Global configuration
################################################################################

message(paste(rep("-", 80), collapse = ""))
message("\tGlobal configuration")

options(
  scipen = 99,
  stringsAsFactors = FALSE
)

set.seed(123)

################################################################################
# Simulation parameters
################################################################################

mc_chart <- 500
mc_reps <- 1000
n1 <- 40
n2 <- 20
K <- 20
alpha <- 0.05
rho <- 0

################################################################################
# 2. Simulation scenarios: 2A
################################################################################

tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3 / 2)
var.teor <- 1
corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2)
)


################################################################################
# Functional data generators
################################################################################

message("\tCreating functional data generators")

func.sim.set <- function(
  t,
  var.teor = 1,
  trend.teor,
  corr.teor,
  rho = 0
) {
  mdata <- length(t)

  sd.teor <- sqrt(as.numeric(var.teor))

  if (rho != 0) {
    corr.teor <- corr.teor *
      sqrt((1 + rho) / (1 - rho))
  }

  C <- svd(t(corr.teor))

  L.corr.teor <- C$u %*% diag(sqrt(C$d))

  func.sim <- function(rep) {
    err.norm <- matrix(
      rnorm(mdata * rep),
      nrow = mdata
    )

    data.err <- L.corr.teor %*% err.norm

    if (rho != 0) {
      data.err[, 1] <-
        data.err[, 1] *
        sqrt((1 - rho) / (1 + rho))

      for (i in 2:rep) {
        data.err[, i] <-
          rho * data.err[, i - 1] + (1 - rho) * data.err[, i]
      }
    }

    res <- as.numeric(trend.teor) + data.err

    return(res)
  }

  return(func.sim)
}


f0 <- func.sim.set(
  tt,
  var.teor,
  trend.teor = mu0,
  corr.teor,
  rho
)

################################################################################
# Monte Carlo
################################################################################

#etas <- c(0, 0.3, 0.5, 0.7, 0.9)
etas <- c(0, 0.02, 0.03, 0.05)

potencia_l2std_montecarlo_40_20 <- numeric(length(etas))
senal_l2std_montecarlo_40_20 <- vector("list", length(etas))

for (i in seq_along(etas)) {
  eta <- etas[i]
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
        fdahotelling:::stat_L2_std(
          x = as.matrix(calibrado_h0),
          y = as.matrix(matriz_grupo_mc)
        )
    }

    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )

    senal_grafico_l2std_montecarlo_40_20 <- logical(K)

    for (k in seq_len(K)) {
      if (eta == 0) {
        matriz_grupo <- t(f0(n2))
      } else {
        matriz_grupo <- t(f1(n2))
      }

      l2std <-
        fdahotelling:::stat_L2_std(
          x = calibrado_h0,
          y = matriz_grupo
        )

      senal_grafico_l2std_montecarlo_40_20[k] <- l2std > UCL
    }

    senal_eta[[g]] <- senal_grafico_l2std_montecarlo_40_20
  }

  senal_l2std_montecarlo_40_20[[i]] <- unlist(senal_eta)

  potencia_l2std_montecarlo_40_20[i] <- mean(senal_l2std_montecarlo_40_20[[i]])
}


save(
  potencia_l2std_montecarlo_40_20,
  senal_l2std_montecarlo_40_20,
  file = "l2std_montecarlo_40_20_2A.RData"
)
