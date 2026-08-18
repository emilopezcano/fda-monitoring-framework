# Read data files

## List of files to be read
ifiles <- list.files(
  path = "data/case_study_1",
  pattern = "\\.xlsx$",
  full.names = TRUE
)

## Read files and combine into a single data.table
ldata <- lapply(
  ifiles,
  readxl::read_excel,
  sheet = 1,
  col_names = TRUE,
  guess_max = 24000
)

data_contaminacion <- data.table::rbindlist(
  ldata,
  fill = TRUE
)

## Set appropriate names (in Spanish)
setnames(
  data_contaminacion,
  c(
    "fecha",
    "hora",
    "magnitud",
    "tipo_estacion",
    "zona",
    "estacion",
    "concentracion_horaria"
  )
)

## Convert date column to Date type
data_contaminacion[,
  fecha := as.Date(fecha)
]

# Data Processing

## Parameters for data processing

### Concentration limit for filtering observations
limite_concentracion <- 200

### Study area
zona_estudio <- "Interior M-30"

### Excluded stations
estaciones_excluidas <- c(
  "Plaza de España",
  "Barrio del Pilar"
)

### Days of January and February 2020
dias_enero <- as.Date(c(
  "2020-01-01",
  "2020-01-02",
  "2020-01-03",
  "2020-01-06",
  "2020-01-07",
  "2020-01-08",
  "2020-01-09",
  "2020-01-10",
  "2020-01-13",
  "2020-01-14",
  "2020-01-15",
  "2020-01-16",
  "2020-01-17",
  "2020-01-20",
  "2020-01-21",
  "2020-01-22",
  "2020-01-23",
  "2020-01-24",
  "2020-01-27",
  "2020-01-28",
  "2020-01-29",
  "2020-01-30",
  "2020-01-31"
))
dias_febrero <- as.Date(c(
  "2020-02-03",
  "2020-02-04",
  "2020-02-05",
  "2020-02-06",
  "2020-02-07",
  "2020-02-10",
  "2020-02-11",
  "2020-02-12",
  "2020-02-13",
  "2020-02-14",
  "2020-02-17",
  "2020-02-18",
  "2020-02-19",
  "2020-02-20",
  "2020-02-21",
  "2020-02-24",
  "2020-02-25",
  "2020-02-26",
  "2020-02-27",
  "2020-02-28"
))


## Filtering observations

data_tratamiento <- copy(data_contaminacion)

data_tratamiento <- data_tratamiento[
  concentracion_horaria < limite_concentracion |
    is.na(concentracion_horaria)
]
data_interior_m30 <- data_tratamiento[
  zona == zona_estudio &
    !estacion %in% estaciones_excluidas
]
interior_m30_enero <- data_interior_m30[
  fecha %in% dias_enero
]
interior_m30_febrero <- data_interior_m30[
  fecha %in% dias_febrero
]

## Get matrices for each month

### January
fechas_enero <- sort(unique(interior_m30_enero$fecha))

orden_estaciones <- sort(
  unique(interior_m30_enero$estacion)
)

matrices_enero <- lapply(
  fechas_enero,
  function(f) {
    datos_dia <- interior_m30_enero[
      fecha == f
    ]

    create_functional_matrix(
      datos_dia,
      orden_estaciones
    )
  }
)

names(matrices_enero) <- as.character(fechas_enero)

matrices_enero_interp <- lapply(
  matrices_enero,
  interpolate_matrix
)

### February

fechas_febrero <- sort(
  unique(interior_m30_febrero$fecha)
)

matrices_febrero <- lapply(
  fechas_febrero,
  function(f) {
    datos_dia <- interior_m30_febrero[
      fecha == f
    ]
    create_functional_matrix(
      datos_dia,
      orden_estaciones
    )
  }
)

names(matrices_febrero) <- as.character(fechas_febrero)

matrices_febrero_interp <- lapply(
  matrices_febrero,
  interpolate_matrix
)

matrices_phase1 <- c(
  matrices_enero_interp,
  matrices_febrero_interp
)


matrices_phase1 <- c(
  matrices_enero_interp,
  matrices_febrero_interp
)


prototype <- tibble(
  hora = 1:24,
  proto = colMeans(
    do.call(rbind, matrices_phase1)
  )
)


data_monitoreo <- data_interior_m30[
  fecha >= "2020-03-02" &
    fecha <= "2020-04-30" &
    fecha != "2020-03-21" &
    fecha != "2020-03-14" &
    fecha != "2020-03-15" &
    fecha != "2020-03-21" &
    fecha != "2020-03-22" &
    fecha != "2020-03-28" &
    fecha != "2020-03-29" &
    fecha != "2020-04-04" &
    fecha != "2020-04-05" &
    fecha != "2020-04-11" &
    fecha != "2020-04-12" &
    fecha != "2020-04-18" &
    fecha != "2020-04-19" &
    fecha != "2020-04-25" &
    fecha != "2020-04-26" &
    fecha != "2020-03-01" &
    fecha != "2020-03-07" &
    fecha != "2020-03-08"
]


gdata <- data_monitoreo |>

  arrange(fecha, estacion, hora) |>

  dplyr::select(
    fecha,
    estacion,
    hora,
    concentracion_horaria
  ) |>

  group_by(fecha, estacion) |>

  mutate(
    concentracion_horaria = na_interpolation(concentracion_horaria)
  ) |>

  ungroup() |>

  group_by(fecha, hora) |>

  summarise(
    concentracion_horaria_med = mean(concentracion_horaria),
    .groups = "drop"
  ) |>

  mutate(
    semana = c(
      rep(1:8, each = 24 * 5),
      rep(9, 4 * 24)
    )
  )

lfacets_df <- gdata |>
  group_by(semana) |>
  summarise(
    periodo = paste(range(fecha), collapse = " to "),
    .groups = "drop"
  )

lfacets <- setNames(
  lfacets_df$periodo,
  lfacets_df$semana
)

dias <- weekdays(gdata$fecha)

if (any(dias == "lunes")) {
  gdata <- gdata |>
    mutate(
      wday = factor(
        dias,
        levels = c(
          "lunes",
          "martes",
          "miércoles",
          "jueves",
          "viernes"
        ),
        labels = c(
          "Monday",
          "Tuesday",
          "Wednesday",
          "Thursday",
          "Friday"
        )
      )
    )
} else {
  gdata <- gdata |>
    mutate(
      wday = factor(
        dias,
        levels = c(
          "Monday",
          "Tuesday",
          "Wednesday",
          "Thursday",
          "Friday"
        )
      )
    )
}

p <- ggplot(
  gdata,
  aes(
    x = hora,
    y = concentracion_horaria_med,
    colour = wday,
    group = fecha
  )
) +

  geom_line(linewidth = 0.6) +

  geom_point(size = 0.2) +

  facet_wrap(
    ~semana,
    ncol = 2,
    labeller = as_labeller(lfacets)
  ) +

  geom_line(
    data = prototype,
    aes(
      x = hora,
      y = proto,
      linewidth = "Phase I mean"
    ),
    inherit.aes = FALSE,
    colour = "black",
    linetype = 2
  ) +

  labs(
    title = "Overview of Phase II case study data",
    x = "Hour",
    y = expression(NO[2] ~ "(" * mu * "g/" * m^3 * ")"),
    colour = "Weekday",
    linewidth = ""
  ) +

  scale_linewidth_manual(
    values = 1,
    breaks = "Phase I mean"
  ) +

  theme_bw() +

  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),

    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),

    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),

    strip.text = element_text(size = 16, face = "bold")
  )


#-------------------------------------------------------------------------------
# Save figure
#-------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "Weekly grouping of working‐day NO2 profiles, with the Phase I mean profile as reference (dashed line).pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 11,
  height = 10,
  units = "in",
  dpi = 300
) #-------------------------------------------------------------------------------
# Save figure
#-------------------------------------------------------------------------------

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
