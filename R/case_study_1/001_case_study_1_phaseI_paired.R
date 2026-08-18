## Phase I: Plot Figure 3 of the paper and get the reference curve for Phase I

# Phase I reference
ref_phase1 <- Reduce(
  "+",
  matrices_phase1
) /
  length(matrices_phase1)


# Outlier detection
fdchart <- fdqcd(ref_phase1)

fddep <- fdqcs.depth(
  fdchart,
  depth = depth.mode,
  nb = 2000,
  plot = FALSE,
  ns = 0.05
)

p <- plot_depth(fddep)
p

# Save figure

ggsave(
  filename = file.path(
    "results/case_study_1/figures",
    "Phase I depth‐based selection of in‐control curves for the Madrid NO2 data.pdf"
  ),
  plot = p,
  device = cairo_pdf,
  width = 10,
  height = 3.5,
  units = "in",
  dpi = 300
)


# Identify outlying stations, if any
if (length(fddep$out) > 0) {
  estaciones_out <- rownames(ref_phase1)[fddep$out]
} else {
  estaciones_out <- character(0)
}

estaciones_out

# Stations retained

estaciones_clean <- setdiff(
  rownames(ref_phase1),
  estaciones_out
)

estaciones_clean


# Clean Phase I daily matrices

matrices_phase1_clean <- lapply(
  matrices_phase1,
  function(dia) {
    dia[
      estaciones_clean,
      ,
      drop = FALSE
    ]
  }
)

# Clean Phase I reference

ref_phase1_clean <- Reduce(
  "+",
  matrices_phase1_clean
) /
  length(matrices_phase1_clean)
