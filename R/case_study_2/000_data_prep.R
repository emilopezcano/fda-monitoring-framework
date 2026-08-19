# Reading data ----

data21 <- read_xlsx(
  file.path("data/case_study_2/data21.xlsx")
)

data21 <- data.table(data21)

data22 <- read_xlsx(
  file.path("data/case_study_2/data22.xlsx"),
  sheet = "analisis_semanas"
)

# Data processing ----

## Clean names
data21 <- data21 %>%
  clean_names()

data22 <- data22 %>%
  clean_names()

names(data21)[1] <- "fecha"

## Combine data
datos_combinados <- data21 %>%
  left_join(
    data22 %>% dplyr::select(fecha, semana),
    by = "fecha"
  ) %>%
  filter(!is.na(semana))


# Figure 2. Mean power by hour and week ----

## Prepare data for plot
nombre_variable <- names(datos_combinados)[4]
datos_semana_media <- datos_combinados %>%
  group_by(semana, time) %>%
  summarise(
    media = mean(.data[[nombre_variable]], na.rm = TRUE),
    .groups = "drop"
  )

primeras <- datos_semana_media %>%
  filter(semana %in% 1:5)

otras <- datos_semana_media %>%
  filter(semana %in% 7:14)

media_1_5 <- primeras %>%
  group_by(time) %>%
  summarise(
    media = mean(media),
    .groups = "drop"
  )

## Set colors
colores_okabe <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7",
  "#9467BD"
)

names(colores_okabe) <- as.character(7:14)

## Create plot
g2 <- ggplot() +
  
  geom_line(
    data = primeras,
    aes(time, media,
        group = semana,
        color = "Weeks 1–5"),
    linewidth = 0.8,
    alpha = 0.9
  ) +
  
  geom_line(
    data = otras,
    aes(time, media,
        group = semana,
        color = as.factor(semana)),
    linewidth = 0.6,
    alpha = 0.95
  ) +
  
  geom_line(
    data = media_1_5,
    aes(time, media,
        color = "Mean (1–5)"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = media_1_5,
    aes(time, media,
        color = "Mean (1–5)"),
    size = 1.8
  ) +
  
  scale_color_manual(
    values = c(
      "Weeks 1–5" = "#8C8C8C",
      colores_okabe,
      "Mean (1–5)" = "#E41A1C"
    ),
    breaks = c(
      "Weeks 1–5",
      as.character(7:14),
      "Mean (1–5)"
    ),
    labels = c(
      "Weeks 1–5",
      paste("Week", 7:14),
      "Mean (1–5)"
    ),
    name = NULL
  ) +
  
  labs(
    x = "Hour of the day",
    y = "Power"
  ) +
  
  theme_bw(base_family = "Arial") +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.35, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
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


## Save plot
ggsave(
  filename = file.path("results/case_study_2/figures/power_vs_hour_of_the_day.pdf"),
  plot = g2,
  device = cairo_pdf,
  width = 7,
  height = 3.5,
  units = "in",
  dpi = 300
)



