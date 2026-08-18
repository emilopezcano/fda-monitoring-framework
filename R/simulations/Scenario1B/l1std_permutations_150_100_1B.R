ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoSNOW(cl)


deltas <- c(0,0.3,0.5,0.7,0.9)
potencia_l1std_perms_150_100 <- numeric(length(deltas))
senal_l1std_perms_150_100 <- vector("list", length(deltas))
alpha<-0.05

seeds_perm <- 1000 + seq_along(deltas)
mc <- 1000
P  <- 500
K  <- 20
n1<-150
n2<-100



tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3/2)
var.teor <- 1
corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2) * (s + 0.5) * (t + 0.5)
)


rho <- 0


for(i in seq_along(deltas)) 
{
  
  delta <- deltas[i]
  
  
  senal <- foreach(
    g = seq_len(mc),
    .packages = c("fda","fda.usc","fdahotelling"),
    .export = c(
      "fdqcd",
      "fdqcs.depth",
      "fdqcs.depth.default",
      "fdqcs.depth.fdqcd",
      "func.sim.set",
      "mu0"
    ),
    .options.RNG = seeds_perm[i]
  ) %dorng% {
    
    
    f0 <- func.sim.set(
      t          = tt,
      var.teor   = var.teor,
      trend.teor = mu0,
      corr.teor  = corr.teor,
      rho        = rho
    )
    
    f1 <- func.sim.set(
      t          = tt,
      var.teor   = var.teor,
      trend.teor = mu0 + delta,
      corr.teor  = corr.teor,
      rho        = rho
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
    
    if(delta != 0){
      
      fdchart <- fdqcd(as.data.frame(calibrado))
      
      fddep <- fdqcs.depth(
        fdchart,
        depth = fda.usc::depth.mode,
        nb = 2000,
        plot = FALSE,
        ns = 0.05
      )
      
      if(length(fddep$out) > 0){
        
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
    
    for(p in seq_len(P)){
      
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
    
    
    
    estad_monitor <- numeric(K)
    
    
    for(k in seq_len(K)){
      
      if(delta == 0){
        
        f0 <- func.sim.set(
          t          = tt,
          var.teor   = var.teor,
          trend.teor = mu0,
          corr.teor  = corr.teor,
          rho        = rho
        )
        
        matriz_grupo <- t(f0(n2))
        
      }else{
        
        f1 <- func.sim.set(
          t          = tt,
          var.teor   = var.teor,
          trend.teor = mu0 + delta,
          corr.teor  = corr.teor,
          rho        = rho
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
      per =  perm_estadistico
    )
    
    
  } 
  
  senal_l1std_perms_150_100[[i]]<- unlist(lapply(senal, function(x) x$s))
  
  
  potencia_l1std_perms_150_100[i] <- mean(senal_l1std_perms_150_100[[i]])
}


save(
  potencia_l1std_perms_150_100,
  senal_l1std_perms_150_100,
  file = "l1std_perms_150_100_1B.RData"
)

stopCluster(cl)


