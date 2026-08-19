## This script generates phase II monitoring statistics and plots

# Phase II Hotelling T2 control chart ----

## Parameters ----

n1 <- nrow(calibrado)
n2 <- 5
n <- n1 + n2
mu <- rep(0, ncol(calibrado))
nb <- 5000
smo <- 0.01
tol <- 1e-3


## Bootstrap-based control limit ----

set.seed(params$seed$case_study_2$hotelling)

boot_estadistico <- rep(NA_real_, nb)

for (i in seq_len(nb)) {
  # 1. Resample from the single Phase I sample

  ind <- sample(
    seq_len(nrow(calibrado)),
    size = n,
    replace = TRUE
  )

  data_boot <- as.matrix(
    calibrado[ind, -1]
  )

  # 2. Estimate covariance

  cov.est <- var(data_boot)

  eig <- eigen(
    cov.est,
    symmetric = TRUE
  )

  lambda <- eig$values
  V <- eig$vectors

  if (!all(lambda >= -tol * abs(lambda[1]))) {
    next
  }

  # 3. Smoothed bootstrap

  L.boot <- V %*%
    diag(
      sqrt(pmax(lambda * smo, 0))
    )

  dat <- data_boot +
    t(
      L.boot %*%
        matrix(
          rnorm(n * ncol(data_boot)),
          nrow = ncol(data_boot)
        )
    )

  # 4. Split into two independent bootstrap samples

  boot_1 <- dat[
    seq_len(n1),
    ,
    drop = FALSE
  ]

  boot_2 <- dat[
    (n1 + 1):n,
    ,
    drop = FALSE
  ]

  # 5. Hotelling T2

  boot_estadistico[i] <-
    fdahotelling:::stat_hotelling_impl(
      x = as.matrix(boot_1),
      y = as.matrix(boot_2),
      mu = mu,
      paired = FALSE,
      step_size = 0.02,
      use_correction = FALSE,
      tolerance = tol
    )
}


## Upper control limit ----

UCL <- quantile(
  boot_estadistico,
  probs = 0.95,
  na.rm = TRUE
)

UCL

## Prepare data for the plot ----

estad_semana <- sapply(7:14, function(sem) {
  matriz_semana <- datos_combinados %>%
    filter(semana == sem) %>%
    group_by(fecha, time) %>%
    summarise(
      valor = mean(.data[[nombre_variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = valor
    ) %>%
    arrange(fecha)

  fdahotelling:::stat_hotelling_impl(
    x = as.matrix(calibrado[, -1]),
    y = as.matrix(matriz_semana[, -1]),
    mu = mu,
    paired = FALSE,
    step_size = 0.02,
    use_correction = FALSE,
    tolerance = tol
  )
})

df_plot <- tibble(
  Semana = 7:14,
  Estadistico = estad_semana
)

## Create plot ----
p <- ggplot(df_plot, aes(Semana, Estadistico)) +

  geom_line(linewidth = 0.5) +

  geom_point(
    aes(color = Estadistico >= UCL),
    size = 1.8
  ) +

  geom_hline(
    yintercept = UCL,
    colour = "#E41A1C",
    linetype = "dashed",
    linewidth = 0.6
  ) +

  annotate(
    "text",
    x = max(df_plot$Semana),
    y = UCL + 0.02 * max(df_plot$Estadistico, na.rm = TRUE),
    label = paste0("UCL = ", round(UCL, 4)),
    hjust = 1,
    vjust = 0,
    colour = "#E41A1C",
    size = 3
  ) +

  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#E41A1C"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = 7:14
  ) +

  labs(
    title = "Phase II Monitoring",
    x = "Week",
    y = expression(T^2 ~ "Statistic")
  ) +

  theme_bw(base_family = "Arial") +

  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

p

## Save figure ----

ggsave(
  filename = file.path(
    "results/case_study_2/figures",
    "T2_phaseII_monitoring_case_study2_2.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


# Phase II L1 control chart ----

## Parameters ----

n1 <- nrow(calibrado)
n2 <- 5
n <- n1 + n2
mu <- rep(0, ncol(calibrado))
nb <- 5000
smo <- 0.01

## Bootstrap-based control limit ----

set.seed(params$seed$case_study_2$hotelling)

boot_estadistico <- rep(NA_real_, nb)


for (i in seq_len(nb)) {
  # 1. Resample from the single Phase I sample
  ind <- sample(
    seq_len(nrow(calibrado)),
    size = n,
    replace = TRUE
  )

  data_boot <- as.matrix(
    calibrado[ind, -1]
  )

  # 2. Estimate covariance
  cov.est <- var(data_boot)

  eig <- eigen(
    cov.est,
    symmetric = TRUE
  )

  lambda <- eig$values
  V <- eig$vectors

  if (!all(lambda >= -tol * abs(lambda[1]))) {
    next
  }

  # 3. Smoothed bootstrap
  L.boot <- V %*%
    diag(
      sqrt(pmax(lambda * smo, 0))
    )

  dat <- data_boot +
    t(
      L.boot %*%
        matrix(
          rnorm(n * ncol(data_boot)),
          nrow = ncol(data_boot)
        )
    )
  # 4. Split into two independent bootstrap samples
  boot_1 <- dat[
    seq_len(n1),
    ,
    drop = FALSE
  ]

  boot_2 <- dat[
    (n1 + 1):n,
    ,
    drop = FALSE
  ]

  # 5. L1 standardized

  boot_estadistico[i] <-
    fdahotelling:::stat_L1_std(
      x = as.matrix(boot_1),
      y = as.matrix(boot_2)
    )
}


## Upper control limit ----

UCL <- quantile(
  boot_estadistico,
  probs = 0.95,
  na.rm = TRUE
)

UCL

## Prepare data for plots ----

estad_semana <- sapply(7:14, function(sem) {
  matriz_semana <- datos_combinados %>%
    filter(semana == sem) %>%
    group_by(fecha, time) %>%
    summarise(
      valor = mean(.data[[nombre_variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = valor
    ) %>%
    arrange(fecha)

  fdahotelling:::stat_L2_std(
    x = as.matrix(calibrado[, -1]),
    y = as.matrix(matriz_semana[, -1]),
  )
})

df_plot <- tibble(
  Semana = 7:14,
  Estadistico = estad_semana
)

## Create plot ----

p <- ggplot(df_plot, aes(Semana, Estadistico)) +

  geom_line(linewidth = 0.5) +

  geom_point(
    aes(color = Estadistico >= UCL),
    size = 1.8
  ) +

  geom_hline(
    yintercept = UCL,
    colour = "#E41A1C",
    linetype = "dashed",
    linewidth = 0.6
  ) +

  annotate(
    "text",
    x = max(df_plot$Semana),
    y = UCL + 0.02 * max(df_plot$Estadistico, na.rm = TRUE),
    label = paste0("UCL = ", round(UCL, 4)),
    hjust = 1,
    vjust = 0,
    colour = "#E41A1C",
    size = 3
  ) +

  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#E41A1C"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = 7:14
  ) +

  labs(
    title = "Phase II Monitoring",
    x = "Week",
    y = expression(L^1 ~ "Statistic")
  ) +

  theme_bw(base_family = "Arial") +

  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

p


## Save figure ----

ggsave(
  filename = file.path(
    "results/case_study_2/figures/",
    "L1std_phaseII_monitoring_case_study2_2.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)

# Phase II L2 standardized ----

## Bootstrap-based control limit ----

set.seed(params$seed$case_study_2$l2_std)

boot_estadistico <- rep(NA_real_, nb)


for (i in seq_len(nb)) {
  # 1. Resample from the single Phase I sample

  ind <- sample(
    seq_len(nrow(calibrado)),
    size = n,
    replace = TRUE
  )

  data_boot <- as.matrix(
    calibrado[ind, -1]
  )

  # 2. Estimate covariance

  cov.est <- var(data_boot)

  eig <- eigen(
    cov.est,
    symmetric = TRUE
  )

  lambda <- eig$values
  V <- eig$vectors

  if (!all(lambda >= -tol * abs(lambda[1]))) {
    next
  }

  # 3. Smoothed bootstrap

  L.boot <- V %*%
    diag(
      sqrt(pmax(lambda * smo, 0))
    )

  dat <- data_boot +
    t(
      L.boot %*%
        matrix(
          rnorm(n * ncol(data_boot)),
          nrow = ncol(data_boot)
        )
    )

  # 4. Split into two independent bootstrap samples

  boot_1 <- dat[
    seq_len(n1),
    ,
    drop = FALSE
  ]

  boot_2 <- dat[
    (n1 + 1):n,
    ,
    drop = FALSE
  ]

  # 5. L2std

  boot_estadistico[i] <-
    fdahotelling:::stat_L2_std(
      x = as.matrix(boot_1),
      y = as.matrix(boot_2)
    )
}


## Upper control limit ----

UCL <- quantile(
  boot_estadistico,
  probs = 0.95,
  na.rm = TRUE
)

UCL

## prepare data for the plot ----

estad_semana <- sapply(7:14, function(sem) {
  matriz_semana <- datos_combinados %>%
    filter(semana == sem) %>%
    group_by(fecha, time) %>%
    summarise(
      valor = mean(.data[[nombre_variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = valor
    ) %>%
    arrange(fecha)

  fdahotelling:::stat_L2_std(
    x = as.matrix(calibrado[, -1]),
    y = as.matrix(matriz_semana[, -1]),
  )
})


df_plot <- tibble(
  Semana = 7:14,
  Estadistico = estad_semana
)

## Create plot ----
p <- ggplot(df_plot, aes(Semana, Estadistico)) +

  geom_line(linewidth = 0.5) +

  geom_point(
    aes(color = Estadistico >= UCL),
    size = 1.8
  ) +

  geom_hline(
    yintercept = UCL,
    colour = "#E41A1C",
    linetype = "dashed",
    linewidth = 0.6
  ) +

  annotate(
    "text",
    x = max(df_plot$Semana),
    y = UCL + 0.02 * max(df_plot$Estadistico, na.rm = TRUE),
    label = paste0("UCL = ", round(UCL, 4)),
    hjust = 1,
    vjust = 0,
    colour = "#E41A1C",
    size = 3
  ) +

  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#E41A1C"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = 7:14
  ) +

  labs(
    title = "Phase II Monitoring",
    x = "Week",
    y = expression(L^2 ~ "Statistic")
  ) +

  theme_bw(base_family = "Arial") +

  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

p

## Save figure ----

ggsave(
  filename = file.path(
    "results/case_study_2/figures",
    "L2std_phaseII_monitoring_case_study2_2.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)


#  Phase II  Q control chart ----

set.seed(params$seed$case_study_1$q)

## Prepare data to compute statistic ----
nombre_variable <- names(datos_combinados)[4]

fdata_ref <- fdata(
  as.matrix(calibrado[, -1]),
  argvals = seq_len(ncol(calibrado[, -1]))
)

semanas <- 7:14

Q_vals <- sapply(semanas, function(sem) {
  mon_df <- datos_combinados %>%
    filter(semana == sem) %>%
    group_by(fecha, time) %>%
    summarise(
      valor = mean(.data[[nombre_variable]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = time,
      values_from = valor
    )

  if (nrow(mon_df) <= 1) {
    return(NA)
  }

  fdata_mon <- fdata(
    as.matrix(mon_df[, -1]),
    argvals = seq_len(ncol(mon_df[, -1]))
  )

  mean(
    depth.FM(
      fdata_mon,
      fdata_ref
    )$dep
  )
})

## Compute limit ----
n <- 5
alpha <- 0.05
LCL <- (factorial(n) * alpha / n^n)^(1 / n)
LCL


## Prepare data for the plot ----

df_plot <- tibble(
  Week = semanas,
  Q = Q_vals,
  OutOfControl = Q_vals < LCL
)

## Create plot ----
p <- ggplot(
  df_plot,
  aes(Week, Q)
) +
  geom_line(linewidth = 0.6) +
  geom_point(
    aes(color = OutOfControl),
    size = 2
  ) +
  geom_hline(
    yintercept = LCL,
    colour = "#E41A1C",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = max(df_plot$Week),
    y = LCL + 0.02 * max(df_plot$Q, na.rm = TRUE),
    label = paste0("LCL = ", round(LCL, 4)),
    hjust = 1,
    colour = "#E41A1C",
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#E41A1C"
    ),
    guide = "none"
  ) +
  labs(
    title = "Phase II Monitoring",
    x = "Week",
    y = "Q statistic"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.border = element_rect(
      colour = "grey70",
      fill = NA,
      linewidth = 0.5
    )
  )

p

ggsave(
  filename = file.path(
    "results/case_study_2/figures",
    "Q_phaseII_monitoring_case_study2_2.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)
