set.seed(123)

################################################################################
# Simulation parameters
################################################################################

mc_chart <- 500   
mc_reps<-1000
n1 <- 150   
n2 <- 100   
K <- 20
tol <- 1e-6
alpha <- 0.05
rho <- 0

################################################################################
# 2. Simulation scenarios: 1A
################################################################################

tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3/2)
var.teor <- 1
corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2)
)



f0 <- func.sim.set(
  tt,
  var.teor,
  trend.teor = mu0,
  corr.teor,
  rho
)


################################################################################
# Monte Carlo l1std
################################################################################

deltas <- c(0, 0.3, 0.5, 0.7, 0.9)

potencia_l1std_montecarlo_150_100 <- numeric(length(deltas))
senal_l1std_montecarlo_150_100 <- vector("list", length(deltas))

for(i in seq_along(deltas)) 
{
  delta <- deltas[i]
  
  senal_delta <- vector("list", mc_chart)
  
  f1 <- func.sim.set(
    t          = tt,
    var.teor   = var.teor,
    trend.teor = mu0 + delta,
    corr.teor  = corr.teor,
    rho        = rho
  )
  
  for(g in seq_len(mc_chart)) 
  {
    
    ##########################################################################
    # 1. UCL Estimation
    ##########################################################################
    
    estadistico_h0 <- numeric(mc_reps)
    
    # Phase I bajo H0: fija para este gráfico
    calibrado_h0 <- t(f0(n1))
    
    for(b in seq_len(mc_reps))
    {
      
      # Phase II bajo H0
      matriz_grupo_mc <- t(f0(n2))
      
      estadistico_h0[b] <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo_mc
        )
    }
    
    ##########################################################################
    # UCL específico del gráfico
    ##########################################################################
    
    UCL <- quantile(
      estadistico_h0,
      probs = 1 - alpha
    )
    
    ##########################################################################
    # 2. Monitoring Phase II
    ##########################################################################
    
    senal_grafico_l1std_montecarlo_150_100 <- logical(K)
    
    for(k in seq_len(K))
    {
      
      if(delta == 0){
        
        matriz_grupo <- t(f0(n2))
        
      }else{
        
        matriz_grupo <- t(f1(n2))
        
      }
      
      l1std <-
        fdahotelling:::stat_L1_std(
          x = calibrado_h0,
          y = matriz_grupo
        )
      
      senal_grafico_l1std_montecarlo_150_100[k] <-
        l1std > UCL
    }
    

    
    senal_delta[[g]] <-
      senal_grafico_l1std_montecarlo_150_100
  }
  
  
  senal_l1std_montecarlo_150_100[[i]] <-
    unlist(senal_delta)
  
  potencia_l1std_montecarlo_150_100[i] <-
    mean(senal_l1std_montecarlo_150_100[[i]])
}

################################################################################
# Save
################################################################################

save(
  potencia_l1std_montecarlo_150_100,
  senal_l1std_montecarlo_150_100,
  file = "l1std_montecarlo_150_100.RData"
)