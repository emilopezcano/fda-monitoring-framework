############################################################
#   ARL SIMULATION FOR FUNCTIONAL HOTELLING T2 WITH PCA
############################################################

library(MASS)

# 1. Function to generate functional data ----

generate.sample <- function(n, shift = 0){
  grid <- seq(0, 1, length = 100)
  mu <- sin(2 * pi * grid) + shift
  
  cov.mat <- outer(grid, grid, function(s,t) exp(-abs(s-t)/0.2))
  L <- chol(cov.mat + 1e-6 * diag(100))
  
  data <- matrix(0, nrow = n, ncol = length(grid))
  for(i in 1:n){
    data[i, ] <- mu + L %*% rnorm(length(grid), sd = 0.3)
  }
  
  return(data)
}

# 2. Funtcion to generate Hotelling T2 data in PCA score space ---

T2.fun <- function(Xs, Ys){
  n1 <- nrow(Xs)
  n2 <- nrow(Ys)
  
  m1 <- colMeans(Xs)
  m2 <- colMeans(Ys)
  
  S1 <- cov(Xs)
  S2 <- cov(Ys)
  
  Spooled <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  
  T2 <- t(m1 - m2) %*% solve(Spooled) %*% (m1 - m2)
  return(as.numeric(T2))
}

# 3. PHASE I CALIBRATION ----

set.seed(123)

n_phase1 <- 80
phase1 <- generate.sample(n_phase1)

# PCA on Phase I curves
pca.obj <- prcomp(phase1, scale. = TRUE)
scores1 <- pca.obj$x[, 1:5]   # first 5 components

B.boot <- 1500
boot.stat <- numeric(B.boot)

for(b in 1:B.boot){
  X <- scores1[sample(1:n_phase1, n_phase1, replace = TRUE), ]
  Y <- scores1[sample(1:n_phase1, 10, replace = TRUE), ]
  
  boot.stat[b] <- T2.fun(X, Y)
}

UCL <- quantile(boot.stat, 0.95)
cat("UCL =", UCL, "\n")

# 4. ARL0 SIMULATION ----

B <- 1000
RL0 <- numeric(B)

for(r in 1:B){
  rl <- 0
  
  repeat{
    rl <- rl + 1
    
    subgroup <- generate.sample(10, shift = 0)
    scores2 <- predict(pca.obj, newdata = subgroup)[, 1:5]
    
    stat <- T2.fun(scores1, scores2)
    
    if(stat > UCL) break
    if(rl > 5000) break
  }
  
  RL0[r] <- rl
}

ARL0 <- mean(RL0)
SDRL0 <- sd(RL0)

cat("ARL0 =", ARL0, "\n")
cat("SDRL0 =", SDRL0, "\n")

# 5. ARL1 SIMULATION ----

deltas <- c(0.1, 0.2, 0.3, 0.5)
results <- data.frame()

for(delta in deltas){
  RL1 <- numeric(B)
  
  for(r in 1:B){
    rl <- 0
    
    repeat{
      rl <- rl + 1
      
      subgroup <- generate.sample(10, shift = delta)
      scores2 <- predict(pca.obj, newdata = subgroup)[, 1:5]
      
      stat <- T2.fun(scores1, scores2)
      
      if(stat > UCL) break
      if(rl > 5000) break
    }
    
    RL1[r] <- rl
  }
  
  results <- rbind(
    results,
    data.frame(
      delta = delta,
      ARL1 = mean(RL1),
      SDRL1 = sd(RL1)
    )
  )
}


# 6. PRINT ALL RESULTS ----

cat("\n====================\n")
cat("ARL0 =", ARL0, "\n")
cat("SDRL0 =", SDRL0, "\n")
cat("====================\n\n")

print(results)

