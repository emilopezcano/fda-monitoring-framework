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
mc_reps<-1000
n1 <- 80   
n2 <- 40    
K <- 20   
tol <- 1e-6
alpha <- 0.05
rho <- 0

################################################################################
# 2. Simulation scenarios: 2B
################################################################################

tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3/2)
var.teor <- 1
corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2) * (s + 0.5) * (t + 0.5)
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
){
  
  mdata <- length(t)
  
  sd.teor <- sqrt(as.numeric(var.teor))
  
  if(rho != 0)
    corr.teor <- corr.teor *
    sqrt((1 + rho)/(1 - rho))
  
  C <- svd(t(corr.teor))
  
  L.corr.teor <- C$u %*% diag(sqrt(C$d))
  
  func.sim <- function(rep){
    
    err.norm <- matrix(
      rnorm(mdata * rep),
      nrow = mdata
    )
    
    data.err <- L.corr.teor %*% err.norm
    
    if(rho != 0){
      
      data.err[,1] <-
        data.err[,1] *
        sqrt((1-rho)/(1+rho))
      
      for(i in 2:rep){
        
        data.err[,i] <-
          rho * data.err[,i-1] +
          (1-rho) * data.err[,i]
        
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
# Monte Carlo T2
################################################################################

#etas <- c(0, 0.3, 0.5, 0.7, 0.9)
etas <- c(0, 0.02, 0.03, 0.05)

potencia_t2_montecarlo_80_40 <- numeric(length(etas))
senal_t2_montecarlo_80_40 <- vector("list", length(etas))

for(i in seq_along(etas)) 
{
  eta <- etas[i]
  senal_eta <- vector("list", mc_chart)
  f1 <- func.sim.set(
    t          = tt,
    var.teor   = var.teor,
    trend.teor = (1-eta)*30 * tt * (1 - tt)^(3/2)+eta*30*tt^(3/2)*(1-tt),
    corr.teor  = corr.teor,
    rho        = rho
  )
  
  for(g in seq_len(mc_chart)) 
  {
    
    
    estadistico_h0 <- numeric(mc_reps)
    calibrado_h0 <- t(f0(n1))
    mu<-rep(0,ncol(calibrado_h0))
    
    for(b in seq_len(mc_reps))
    {
      
      matriz_grupo_mc <- t(f0(n2))
      
      estadistico_h0[b] <-
        fdahotelling:::stat_hotelling_impl(
          x = as.matrix(calibrado_h0),
          y = as.matrix(matriz_grupo_mc),
          mu = mu,
          paired = FALSE,
          step_size = 0.02,
          use_correction = FALSE,
          tolerance = tol
        )
    }
    
    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )
    
    senal_grafico_t2_montecarlo_80_40 <- logical(K)
    
    for(k in seq_len(K)){
      
      if(eta == 0){
        
        matriz_grupo <- t(f0(n2))
        
      }else{
        
        matriz_grupo <- t(f1(n2))
        
      }
      
      T2 <-
        fdahotelling:::stat_hotelling_impl(
          x = calibrado_h0,
          y = matriz_grupo,
          mu = mu,
          paired = FALSE,
          step_size = 0.02,
          use_correction = FALSE,
          tolerance = tol
        )
      
      senal_grafico_t2_montecarlo_80_40[k] <- T2 > UCL
      
    }
    
    
    senal_eta[[g]] <- senal_grafico_t2_montecarlo_80_40
    
  }
  
  senal_t2_montecarlo_80_40[[i]] <- unlist(senal_eta)
  
  potencia_t2_montecarlo_80_40[i] <- mean(senal_t2_montecarlo_80_40[[i]])
}


save(
  potencia_t2_montecarlo_80_40,
  senal_t2_montecarlo_80_40,
  file = "t2_montecarlo_80_40_2B.RData"
)
