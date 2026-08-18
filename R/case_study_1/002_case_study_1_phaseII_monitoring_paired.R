# Phase I daily residuals

## fore each Phase I day, how much each station deviates from its reference profile.
residuos_phase1_clean <- lapply(
  matrices_phase1_clean,
  function(dia) {
    
    # Check station pairing
    stopifnot(
      identical(
        rownames(dia),
        rownames(ref_phase1_clean)
      )
    )
    
    dia - ref_phase1_clean
  }
)

residuos_todos_clean <- do.call(
  rbind,
  residuos_phase1_clean
)

# Covariance matrix

cov_res <- var(residuos_todos_clean)

dim(cov_res)

eig_res <- eigen(
  cov_res,
  symmetric = TRUE
)

lambda_res <- eig_res$values
V_res <- eig_res$vectors


smo <- 0.02

# Smoothing matrix

L_res <- V_res %*%
  diag(
    sqrt(lambda_res * smo),
    nrow = length(lambda_res)
  )


# Smoothed bootstrap

nb <- 5000
tol <- 1e-6
smo <- 0.02

n_est <- nrow(ref_phase1_clean)
p <- ncol(ref_phase1_clean)

mu <- rep(0, p)

boot_estadistico <- rep(NA_real_, nb)

set.seed(params$seed$case_study_1$hotelling)

for (b in seq_len(nb)) {
  
  # ----------------------------------------------------------
  # 1. Select one complete Phase I day
  # ----------------------------------------------------------
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  # ----------------------------------------------------------
  # 2. Generate smooth functional noise
  # ----------------------------------------------------------
  
  Z <- matrix(
    rnorm(p * n_est),
    nrow = p,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(ref_phase1_clean)
  colnames(ruido) <- colnames(ref_phase1_clean)
  
  
  # ----------------------------------------------------------
  # 3. Smoothed bootstrap residual
  # ----------------------------------------------------------
  
  res_star <- res_boot + ruido
  
  
  #----------------------------------------------------------------------------
  # 4. Generate hypothetical in-control day
  #-----------------------------------------------------------------------------
  
  dia_star <- ref_phase1_clean + res_star
  
  
  #-----------------------------------------------------------------------------
  # 5. Check pairing
  #-----------------------------------------------------------------------------
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  
  
  #----------------------------------------------------------------------------
  # 6. Paired T2
  #-----------------------------------------------------------------------------
  
  boot_estadistico[b] <-
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


# Upper Control Limit

alpha <- 0.05

UCL <- unname(
  quantile(
    boot_estadistico,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL

# PHASE II - MONITORING DATA

data_monitoreo <- data_interior_m30[
  fecha >= as.Date("2020-03-02") &
    fecha <= as.Date("2020-04-30") &
    as.POSIXlt(fecha)$wday %in% 1:5
]

data_monitoreo <- as.data.table(data_monitoreo)

fecha_monitoreo <- sort(
  unique(data_monitoreo$fecha)
)

# DAILY MONITORING MATRICES

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

names(matrices_monitoreo) <- as.character(fecha_monitoreo)

matrices_monitoreo_interp <- lapply(
  matrices_monitoreo,
  interpolate_matrix
)

# PHASE II T2

mu <- rep(
  0,
  ncol(ref_phase1)
)

estadistico <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    fdahotelling:::stat_hotelling_impl(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1),
      mu = mu,
      paired = TRUE,
      step_size = 0.02,
      use_correction = FALSE,
      tolerance = tol
    )
  }
)

estadistico


## Data for Phase II control chart

data_plot <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(estadistico)
)


## Phase II control chart

g <- ggplot(
  data_plot,
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
    aes(color = estadistico > UCL),
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
    yintercept = UCL,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot$fecha),
    y = UCL + 0.02 * max(data_plot$estadistico, na.rm = TRUE),
    label = paste0(
      "UCL = ",
      round(UCL, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(T^2 ~ "Statistic")
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot$fecha,
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
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8.5),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g

## Save figure

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "T2_phaseII_monitoring_p.pdf"
  ),
  plot = g,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)

# Phase II L1 std control chart

## Smoothed bootstrap

nb <- 5000
tol <- 1e-6
smo <- 0.02

n_est <- nrow(ref_phase1_clean)
p <- ncol(ref_phase1_clean)

mu <- rep(0, p)

boot_estadistico <- rep(NA_real_, nb)

set.seed(params$seed$case_study_1$l1_std)

for (b in seq_len(nb)) {
  
  ### 1. Select one complete Phase I day
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  ### 2. Generate smooth functional noise
  
  Z <- matrix(
    rnorm(p * n_est),
    nrow = p,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(ref_phase1_clean)
  colnames(ruido) <- colnames(ref_phase1_clean)
  
  
  ### 3. Smoothed bootstrap residual
  
  res_star <- res_boot + ruido
  
  
  ### 4. Generate hypothetical in-control day
  
  dia_star <- ref_phase1_clean + res_star
  
  
  ### 5. Check pairing
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  
  
  ### 6. L1std
  
  boot_estadistico[b] <-
    fdahotelling:::stat_L1_std(
      x = as.matrix(dia_star),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
}



## Upper Control Limit
alpha <- 0.05

UCL <- unname(
  quantile(
    boot_estadistico,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL


# PHASE II L1 std

mu <- rep(
  0,
  ncol(ref_phase1)
)

estadistico <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    fdahotelling:::stat_L1_std(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1),
      mu = mu
    )
  }
)

estadistico


# Data for Phase II control chart

data_plot <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(estadistico)
)


# Phase II control chart

g <- ggplot(
  data_plot,
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
    aes(color = estadistico > UCL),
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
    yintercept = UCL,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot$fecha),
    y = UCL + 0.02 * max(data_plot$estadistico, na.rm = TRUE),
    label = paste0(
      "UCL = ",
      round(UCL, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(L^1 ~ "Statistic")
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot$fecha,
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
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8.5),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g

## Save figure

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "L1std_phaseII_monitoring_p.pdf"
  ),
  plot = g,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)

# Phase II L2std control chart

## Smoothed bootstrap

nb <- 5000
tol <- 1e-6
smo <- 0.02

n_est <- nrow(ref_phase1_clean)
p <- ncol(ref_phase1_clean)

mu <- rep(0, p)

boot_estadistico <- rep(NA_real_, nb)

set.seed(params$seed$case_study_1$l2_std)

for (b in seq_len(nb)) {
  
  ### 1. Select one complete Phase I day
  
  j <- sample(
    seq_along(residuos_phase1_clean),
    size = 1
  )
  
  res_boot <- residuos_phase1_clean[[j]]
  
  
  ### 2. Generate smooth functional noise
  
  Z <- matrix(
    rnorm(p * n_est),
    nrow = p,
    ncol = n_est
  )
  
  ruido <- t(
    L_res %*% Z
  )
  
  rownames(ruido) <- rownames(ref_phase1_clean)
  colnames(ruido) <- colnames(ref_phase1_clean)
  ### 3. Smoothed bootstrap residual
  
  res_star <- res_boot + ruido
  
  ### 4. Generate hypothetical in-control day
  
  dia_star <- ref_phase1_clean + res_star
  
  
  ### 5. Check pairing
  
  stopifnot(
    identical(
      rownames(dia_star),
      rownames(ref_phase1_clean)
    )
  )
  ### 6. Paired L1std
  
  boot_estadistico[b] <-
    fdahotelling:::stat_L2_std(
      x = as.matrix(dia_star),
      y = as.matrix(ref_phase1_clean),
      mu = mu
    )
}


## Upper Control Limit

alpha <- 0.05

UCL <- unname(
  quantile(
    boot_estadistico,
    probs = 1 - alpha,
    na.rm = TRUE
  )
)

UCL


## PHASE II l2 std

mu <- rep(
  0,
  ncol(ref_phase1)
)

estadistico <- sapply(
  matrices_monitoreo_interp,
  function(dia) {
    
    fdahotelling:::stat_L2_std(
      x = as.matrix(dia),
      y = as.matrix(ref_phase1),
      mu = mu

    )
  }
)

estadistico


## Data for Phase II control chart

data_plot <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = as.numeric(estadistico)
)


## Phase II control chart

g <- ggplot(
  data_plot,
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
    aes(color = estadistico > UCL),
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
    yintercept = UCL,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot$fecha),
    y = UCL + 0.02 * max(data_plot$estadistico, na.rm = TRUE),
    label = paste0(
      "UCL = ",
      round(UCL, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(L^2 ~ "Statistic")
  ) +
  
  scale_x_date(
    date_labels = "%b-%d",
    breaks = data_plot$fecha,
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
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8.5),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

g

## Save figure

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "L2std_phaseII_monitoring_p.pdf"
  ),
  plot = g,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


# Phase II Q control chart

set.seed(params$seed$case_study_1$q)


## Phase I reference

fdata_mon <- fdata(
  as.matrix(ref_phase1_clean)
)


## Phase II monitoring

estadistico_q <- rep(
  NA_real_,
  length(fecha_monitoreo)
)

for (i in seq_along(fecha_monitoreo)) {
  
  # Data for monitoring day i
  aux1 <- data_monitoreo[
    fecha == fecha_monitoreo[i],
    c("estacion", "hora", "concentracion_horaria")
  ]
  
  # Daily station x hour matrix
  aux2 <- dcast(
    aux1,
    estacion ~ hora,
    value.var = "concentracion_horaria",
    fun.aggregate = mean
  )
  
  
  #---------------------------------------------------------------------------
  # Keep the same clean stations and the same Phase I order
  #---------------------------------------------------------------------------
  
  aux2 <- aux2[
    match(estaciones_clean, aux2$estacion),
  ]
  
  stopifnot(
    identical(
      as.character(aux2$estacion),
      as.character(estaciones_clean)
    )
  )
  
  
  #---------------------------------------------------------------------------
  # Interpolation
  #---------------------------------------------------------------------------
  
  data_monitoreo_interp <- t(
    apply(
      aux2[, -1],
      1,
      na_interpolation
    )
  )
  
  rownames(data_monitoreo_interp) <- aux2$estacion
  
  
  # Check pairing/order with Phase I reference
  stopifnot(
    identical(
      rownames(data_monitoreo_interp),
      rownames(ref_phase1_clean)
    )
  )
  
  
  #---------------------------------------------------------------------------
  # Functional depth
  #---------------------------------------------------------------------------
  
  fdata_ref <- fdata(
    as.matrix(data_monitoreo_interp)
  )
  
  if (nrow(fdata_mon$data) > 1) {
    
    profs <- depth.FM(
      fdata_mon,
      fdata_ref
    )$dep
    
    estadistico_q[i] <- mean(profs)
    
  } else {
    
    estadistico_q[i] <- NA_real_
  }
}

## Lower Control Limit

n <- nrow(ref_phase1_clean)
alpha <- 0.05

LCL <- (
  factorial(n) * alpha / n^n
)^(1 / n)

## Data for plot

data_plot_q <- data.frame(
  fecha = as.Date(fecha_monitoreo),
  estadistico = estadistico_q
)


## Figure

p <- ggplot(
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
    aes(color = estadistico < LCL),
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
    yintercept = LCL,
    linetype = "dashed",
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = max(data_plot_q$fecha),
    y = LCL +
      0.02 * max(data_plot_q$estadistico, na.rm = TRUE),
    label = paste0(
      "LCL = ",
      round(LCL, 4)
    ),
    hjust = 1,
    vjust = 0,
    color = "#E41A1C",
    size = 3
  ) +
  
  labs(
    title = "Phase II Monitoring",
    x = "Date",
    y = expression(Q ~ "Statistic")
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
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8.5),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

## Save figure

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "Q_phaseII_monitoring.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)
