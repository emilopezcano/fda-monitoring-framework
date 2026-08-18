ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoSNOW(cl)

deltas <- c(0,0.3,0.5,0.7,0.9)
senal_l2std_boot_150_100 <- vector("list", length(deltas))
smo<-0.05
potencia_l2std_boot_150_100 <- numeric(length(deltas))
alpha<-0.05
seeds <- 345 + seq_along(deltas)

mc <- 1000
B  <- 150
K  <- 100
tol <- 1e-6
n1 <-150
n2<-100

tt <- seq(0, 1, length.out = 25)
mu0 <- 30 * tt * (1 - tt)^(3/2)

var.teor <- 1

corr.teor <- outer(
  tt,
  tt,
  function(s, t) exp(-2 * (s - t)^2)
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
    .options.RNG = seeds[i]
  ) %dorng%{
    
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
        depth = depth.mode,
        nb = 2000,
        plot = FALSE,
        ns = 0.05
      )
      
      if(length(fddep$out) > 0){
        
        calibrado <- calibrado[-fddep$out, , drop = FALSE]
        
      }
    }
    
    
    n1_boot <- nrow(calibrado)  
    n2_boot <- n2                
    n_boot  <- n1_boot + n2_boot
    
    boot_estadistico <- rep(NA_real_, B)
    
    mu <- rep(0, ncol(calibrado))
    
    
    
    for(b in seq_len(B))
    {
      
      
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
      
      if(!all(lambda >= -tol * abs(lambda[1])))
        next
      
      
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
    
    
    for(k in seq_len(K)){
      
      if(delta == 0){
        
        f0 <- func.sim.set(
          t          = tt,
          var.teor   = var.teor,
          trend.teor = mu0,
          corr.teor  = corr.teor,
          rho        = rho
        )
        
        matriz_grupo <- t(f0(n2_boot))
        
      }else{
        
        f1 <- func.sim.set(
          t          = tt,
          var.teor   = var.teor,
          trend.teor = mu0 + delta,
          corr.teor  = corr.teor,
          rho        = rho
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
  
  senal_l2std_boot_150_100[[i]]  <- unlist(lapply(senal, function(x) x$s))
  potencia_l2std_boot_150_100[i] <-  mean(senal_l2std_boot_150_100[[i]])
}


stopCluster(cl)

save(
  potencia_l2std_boot_150_100,
  senal_l2std_boot_150_100,
  file = "l2std_boot_150_100_1A.RData"
)

