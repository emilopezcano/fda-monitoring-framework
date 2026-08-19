# Calibration sample: daily functional observations from Weeks 1--5

nombre_variable <- names(datos_combinados)[4]

calibrado <- datos_combinados %>%
  filter(semana %in% 1:5) %>%
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
