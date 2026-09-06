# Daily residuals with respect to the clean Phase I reference
residuos_phase1_clean <- lapply(
  matrices_phase1_clean,
  function(dia) {
    
    # Check that the stations are exactly the same
    # and in the same order as in the clean reference
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    dia - ref_phase1_clean
  }
)

# Stack all Phase I residuals
residuos_todos_clean <- do.call(
  rbind,
  residuos_phase1_clean
)


#===============================================================================
# COVARIANCE MATRIX OF PHASE I RESIDUALS
#===============================================================================

cov_res <- var(
  residuos_todos_clean
)

dim(cov_res)

eig_res <- eigen(
  cov_res,
  symmetric = TRUE
)

lambda_res <- eig_res$values
V_res <- eig_res$vectors


#===============================================================================
# SMOOTHING MATRIX
#===============================================================================

# Madrid application
gamma <- 0.02

L_res <- V_res %*%
  diag(
    sqrt(lambda_res * gamma),
    nrow = length(lambda_res)
  )


#===============================================================================
# GENERAL PARAMETERS
#===============================================================================

nb <- 5000
tol <- 1e-6
alpha <- 0.05

n_est <- nrow(
  ref_phase1_clean
)

p_dim <- ncol(
  ref_phase1_clean
)

mu <- rep(
  0,
  p_dim
)


#===============================================================================
# PHASE II: DATA
#===============================================================================

data_monitoreo <- data_interior_m30[
  fecha >= as.Date("2020-03-02") &
    fecha <= as.Date("2020-04-30") &
    as.POSIXlt(fecha)$wday %in% 1:5
]

data_monitoreo <- as.data.table(
  data_monitoreo
)

fecha_monitoreo <- sort(
  unique(data_monitoreo$fecha)
)


#===============================================================================
# PHASE II: CREATE DAILY MATRICES
#===============================================================================

matrices_monitoreo <- lapply(
  fecha_monitoreo,
  function(f) {
    
    datos_dia <- data_monitoreo[
      fecha == f
    ]
    
    create_functional_matrix(
      datos_dia,
      orden_estaciones
    )
  }
)

names(matrices_monitoreo) <- as.character(
  fecha_monitoreo
)


#===============================================================================
# PHASE II: INTERPOLATION
#===============================================================================

matrices_monitoreo_interp <- lapply(
  matrices_monitoreo,
  interpolate_matrix
)


#===============================================================================
# PHASE II: KEEP EXACTLY THE SAME STATIONS AS IN CLEAN PHASE I
#===============================================================================

matrices_monitoreo_interp <- lapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    #---------------------------------------------------------------------------
    # Check that all retained Phase I stations are present
    #---------------------------------------------------------------------------
    
    stopifnot(
      all(
        estaciones_clean %in% rownames(dia)
      )
    )
    
    
    #---------------------------------------------------------------------------
    # Keep only retained Phase I stations
    #---------------------------------------------------------------------------
    
    dia <- dia[
      estaciones_clean,
      ,
      drop = FALSE
    ]
    
    
    #---------------------------------------------------------------------------
    # Check exact station order
    #---------------------------------------------------------------------------
    
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    
    #---------------------------------------------------------------------------
    # Check that Phase I and Phase II have the same number of hours
    #---------------------------------------------------------------------------
    
    stopifnot(
      ncol(dia) ==
        ncol(ref_phase1_clean)
    )
    
    
    #---------------------------------------------------------------------------
    # Check column names/order when available
    #---------------------------------------------------------------------------
    
    if (
      !is.null(colnames(dia)) &&
      !is.null(colnames(ref_phase1_clean))
    ) {
      
      stopifnot(
        identical(
          colnames(dia),
          colnames(ref_phase1_clean)
        )
      )
    }
    
    
    #---------------------------------------------------------------------------
    # Check that there are no remaining missing values
    #---------------------------------------------------------------------------
    
    stopifnot(
      !anyNA(dia)
    )
    
    dia
  }
)


#===============================================================================
# FINAL CONSISTENCY CHECK BEFORE MONITORING
#===============================================================================

stopifnot(
  length(matrices_monitoreo_interp) ==
    length(fecha_monitoreo)
)

stopifnot(
  all(
    vapply(
      matrices_monitoreo_interp,
      function(dia) {
        identical(
          dim(dia),
          dim(ref_phase1_clean)
        )
      },
      logical(1)
    )
  )
)

stopifnot(
  !anyNA(ref_phase1_clean)
)


#===============================================================================
# HOTELLING T2: SMOOTHED BOOTSTRAP
#===============================================================================

boot_t2 <- rep(
  NA_real_,
  nb
)

set.seed(
  params$seed$case_study_1$hotelling
)

for (b in seq_len(nb)) {
  
  #---------------------------------------------------------------------------
  # 1. Select one complete Phase I residual matrix
  #---------------------------------------------------------------------------
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  #---------------------------------------------------------------------------
  # 2. Generate smooth functional noise
  #---------------------------------------------------------------------------
  
  Z <- matrix(
    rnorm(p_dim * n_est),
    nrow = p_dim,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(
    ref_phase1_clean
  )
  
  colnames(ruido) <- colnames(
    ref_phase1_clean
  )
  
  
  #---------------------------------------------------------------------------
  # 3. Smoothed residual
  #---------------------------------------------------------------------------
  
  res_star <- res_boot + ruido
  
  
  #---------------------------------------------------------------------------
  # 4. Generate hypothetical in-control day
  #---------------------------------------------------------------------------
  
  dia_star <- ref_phase1_clean +
    res_star
  
  
  #---------------------------------------------------------------------------
  # 5. Check dimensions/order
  #---------------------------------------------------------------------------
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  
  
  #---------------------------------------------------------------------------
  # 6. Calculate T2
  #---------------------------------------------------------------------------
  
  boot_t2[b] <-
    fdahotelling:::stat_hotelling_impl(
      x = as.matrix(dia_star),
      y = as.matrix(ref_phase1_clean),
      mu = mu,
      paired = TRUE,
      step_size = 0.02,
      use_correction = FALSE,
      tolerance = tol
    )
}


#===============================================================================
# HOTELLING T2: UCL
#===============================================================================

UCL_t2 <- unname(
  quantile(
    boot_t2,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL_t2


#===============================================================================
# HOTELLING T2: PHASE II
#===============================================================================

estadistico_t2 <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    fdahotelling:::stat_hotelling_impl(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1_clean),
      mu = mu,
      paired = TRUE,
      step_size = 0.02,
      use_correction = FALSE,
      tolerance = tol
    )
  }
)

estadistico_t2


#===============================================================================
# HOTELLING T2: PLOT DATA
#===============================================================================

data_plot_t2 <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(
    estadistico_t2
  )
)


#===============================================================================
# HOTELLING T2: CONTROL CHART
#===============================================================================

g_t2 <- ggplot(
  data_plot_t2,
  aes(
    x = fecha,
    y = estadistico,
    group = 1
  )
) +
  
  geom_line(
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      color = estadistico > UCL_t2
    ),
    size = 1.6
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "#E41A1C",
      "FALSE" = "black"
    ),
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = UCL_t2,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot_t2$fecha),
    y = UCL_t2 +
      0.02 *
      max(
        data_plot_t2$estadistico,
        na.rm = TRUE
      ),
    label = paste0(
      "UCL = ",
      round(UCL_t2, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(
      T^2 ~ "Statistic"
    )
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot_t2$fecha,
    expand = c(0.01, 0.01)
  ) +
  
  theme_bw(
    base_family = "Arial"
  ) +
  
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0,
      vjust = 2,
      face = "plain"
    ),
    legend.position = "none",
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 8.5
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor =
      element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g_t2


ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "T2_phaseII_monitoring_p.pdf"
  ),
  plot = g_t2,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


#===============================================================================
# L1 STD: SMOOTHED BOOTSTRAP
#===============================================================================

boot_l1 <- rep(
  NA_real_,
  nb
)

set.seed(
  params$seed$case_study_1$l1_std
)

for (b in seq_len(nb)) {
  
  #---------------------------------------------------------------------------
  # 1. Select one complete Phase I residual matrix
  #---------------------------------------------------------------------------
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  #---------------------------------------------------------------------------
  # 2. Generate smooth functional noise
  #---------------------------------------------------------------------------
  
  Z <- matrix(
    rnorm(p_dim * n_est),
    nrow = p_dim,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(
    ref_phase1_clean
  )
  
  colnames(ruido) <- colnames(
    ref_phase1_clean
  )
  
  
  #---------------------------------------------------------------------------
  # 3. Smoothed residual
  #---------------------------------------------------------------------------
  
  res_star <- res_boot + ruido
  
  
  #---------------------------------------------------------------------------
  # 4. Generate hypothetical in-control day
  #---------------------------------------------------------------------------
  
  dia_star <- ref_phase1_clean +
    res_star
  
  
  #---------------------------------------------------------------------------
  # 5. Check pairing/order
  #---------------------------------------------------------------------------
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  
  
  #---------------------------------------------------------------------------
  # 6. Calculate standardized L1
  #---------------------------------------------------------------------------
  
  boot_l1[b] <-
    fdahotelling:::stat_L1_std(
      x = as.matrix(dia_star),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
}


#===============================================================================
# L1 STD: UCL
#===============================================================================

UCL_l1 <- unname(
  quantile(
    boot_l1,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL_l1


#===============================================================================
# L1 STD: PHASE II
#===============================================================================

estadistico_l1 <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    fdahotelling:::stat_L1_std(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
  }
)

estadistico_l1


#===============================================================================
# L1 STD: PLOT DATA
#===============================================================================

data_plot_l1 <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(
    estadistico_l1
  )
)


#===============================================================================
# L1 STD: CONTROL CHART
#===============================================================================

g_l1 <- ggplot(
  data_plot_l1,
  aes(
    x = fecha,
    y = estadistico,
    group = 1
  )
) +
  
  geom_line(
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      color = estadistico > UCL_l1
    ),
    size = 1.6
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "#E41A1C",
      "FALSE" = "black"
    ),
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = UCL_l1,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot_l1$fecha),
    y = UCL_l1 +
      0.02 *
      max(
        data_plot_l1$estadistico,
        na.rm = TRUE
      ),
    label = paste0(
      "UCL = ",
      round(UCL_l1, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(
      L^1 ~ "Statistic"
    )
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot_l1$fecha,
    expand = c(0.01, 0.01)
  ) +
  
  theme_bw(
    base_family = "Arial"
  ) +
  
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0,
      vjust = 2,
      face = "plain"
    ),
    legend.position = "none",
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 8.5
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor =
      element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g_l1


ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "L1std_phaseII_monitoring_p.pdf"
  ),
  plot = g_l1,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


#===============================================================================
# L2 STD: SMOOTHED BOOTSTRAP
#===============================================================================

boot_l2 <- rep(
  NA_real_,
  nb
)

set.seed(
  params$seed$case_study_1$l2_std
)

for (b in seq_len(nb)) {
  
  #---------------------------------------------------------------------------
  # 1. Select one complete Phase I residual matrix
  #---------------------------------------------------------------------------
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  #---------------------------------------------------------------------------
  # 2. Generate smooth functional noise
  #---------------------------------------------------------------------------
  
  Z <- matrix(
    rnorm(p_dim * n_est),
    nrow = p_dim,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(
    ref_phase1_clean
  )
  
  colnames(ruido) <- colnames(
    ref_phase1_clean
  )
  
  
  #---------------------------------------------------------------------------
  # 3. Smoothed residual
  #---------------------------------------------------------------------------
  
  res_star <- res_boot + ruido
  
  
  #---------------------------------------------------------------------------
  # 4. Generate hypothetical in-control day
  #---------------------------------------------------------------------------
  
  dia_star <- ref_phase1_clean +
    res_star
  
  
  #---------------------------------------------------------------------------
  # 5. Check pairing/order
  #---------------------------------------------------------------------------
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  
  
  #---------------------------------------------------------------------------
  # 6. Calculate standardized L2
  #---------------------------------------------------------------------------
  
  boot_l2[b] <-
    fdahotelling:::stat_L2_std(
      x = as.matrix(dia_star),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
}


#===============================================================================
# L2 STD: UCL
#===============================================================================

UCL_l2 <- unname(
  quantile(
    boot_l2,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL_l2


#===============================================================================
# L2 STD: PHASE II
#===============================================================================

estadistico_l2 <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    fdahotelling:::stat_L2_std(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
  }
)

estadistico_l2


#===============================================================================
# L2 STD: PLOT DATA
#===============================================================================

data_plot_l2 <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(
    estadistico_l2
  )
)


#===============================================================================
# L2 STD: CONTROL CHART
#===============================================================================

g_l2 <- ggplot(
  data_plot_l2,
  aes(
    x = fecha,
    y = estadistico,
    group = 1
  )
) +
  
  geom_line(
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      color = estadistico > UCL_l2
    ),
    size = 1.6
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "#E41A1C",
      "FALSE" = "black"
    ),
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = UCL_l2,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot_l2$fecha),
    y = UCL_l2 +
      0.02 *
      max(
        data_plot_l2$estadistico,
        na.rm = TRUE
      ),
    label = paste0(
      "UCL = ",
      round(UCL_l2, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(
      L^2 ~ "Statistic"
    )
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot_l2$fecha,
    expand = c(0.01, 0.01)
  ) +
  
  theme_bw(
    base_family = "Arial"
  ) +
  
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0,
      vjust = 2,
      face = "plain"
    ),
    legend.position = "none",
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 8.5
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor =
      element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g_l2


ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "L2std_phaseII_monitoring_p.pdf"
  ),
  plot = g_l2,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


#===============================================================================
# Q CHART
#===============================================================================

set.seed(
  params$seed$case_study_1$q
)


#===============================================================================
# Q: PHASE I REFERENCE
#===============================================================================

fdata_mon <- fdata(
  as.matrix(ref_phase1_clean)
)


#===============================================================================
# Q: PHASE II
#===============================================================================

estadistico_q <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    #---------------------------------------------------------------------------
    # Check exact station pairing/order
    #---------------------------------------------------------------------------
    
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    
    #---------------------------------------------------------------------------
    # Check dimensions
    #---------------------------------------------------------------------------
    
    stopifnot(
      identical(
        dim(dia),
        dim(ref_phase1_clean)
      )
    )
    
    
    #---------------------------------------------------------------------------
    # Convert Phase II matrix to functional data
    #---------------------------------------------------------------------------
    
    fdata_ref <- fdata(
      as.matrix(dia)
    )
    
    
    #---------------------------------------------------------------------------
    # Functional depth
    #---------------------------------------------------------------------------
    
    if (nrow(fdata_mon$data) > 1) {
      
      profs <- depth.FM(
        fdata_mon,
        fdata_ref
      )$dep
      
      mean(profs)
      
    } else {
      
      NA_real_
    }
  }
)

estadistico_q


#===============================================================================
# Q: LOWER CONTROL LIMIT
#===============================================================================

n <- nrow(
  ref_phase1_clean
)

LCL_q <- (
  factorial(n) *
    alpha /
    n^n
)^(1 / n)

LCL_q


#===============================================================================
# Q: PLOT DATA
#===============================================================================

data_plot_q <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(
    estadistico_q
  )
)


#===============================================================================
# Q: CONTROL CHART
#===============================================================================

g_q <- ggplot(
  data_plot_q,
  aes(
    x = fecha,
    y = estadistico,
    group = 1
  )
) +
  
  geom_line(
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      color = estadistico < LCL_q
    ),
    size = 1.6
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "#E41A1C",
      "FALSE" = "black"
    ),
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = LCL_q,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot_q$fecha),
    y = LCL_q +
      0.02 *
      max(
        data_plot_q$estadistico,
        na.rm = TRUE
      ),
    label = paste0(
      "LCL = ",
      round(LCL_q, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(
      Q ~ "Statistic"
    )
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot_q$fecha,
    expand = c(0.01, 0.01)
  ) +
  
  theme_bw(
    base_family = "Arial"
  ) +
  
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0
    ),
    legend.position = "none",
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 8.5
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor =
      element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g_q


ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "Q_phaseII_monitoring.pdf"
  ),
  plot = g_q,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


#===============================================================================
# FINAL CONSISTENCY SUMMARY
#===============================================================================

# cat(
#   "\nT2 UCL:",
#   UCL_t2,
#   "\nL1 UCL:",
#   UCL_l1,
#   "\nL2 UCL:",
#   UCL_l2,
#   "\nQ LCL:",
#   LCL_q,
#   "\n"
# )
# 
# cat(
#   "\nNumber of Phase II monitoring days:",
#   length(fecha_monitoreo),
#   "\nNumber of retained stations:",
#   nrow(ref_phase1_clean),
#   "\nNumber of hourly points:",
#   ncol(ref_phase1_clean),
#   "\n"
# )