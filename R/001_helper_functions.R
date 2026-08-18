# Helper Functions

#' Create a functional matrix from the data
#'
#' @param .data data.frame with the input data 
#' @param orden_estaciones vector with the order of the stations
#'
#' @returns a data.frame with the aggregated functional matrix
#'
#' @export
#' @examples
#' 
create_functional_matrix <- function(.data, orden_estaciones) {
  matriz <- dcast(
    .data,
    estacion ~ hora,
    value.var = "concentracion_horaria",
    fun.aggregate = mean
  )
  matriz <- matriz[
    match(orden_estaciones, matriz$estacion),
  ]
  return(matriz)
}


#' Interpolate missing values in a matrix
#'
#' @param .data data.frame with the input data
#'
#' @returns a data.frame with the interpolated matrix
#'
#' @export
#' @examples
interpolate_matrix <- function(.data) {
  int_matrix <- t(
    apply(
      .data[, -1],
      1,
      na_interpolation
    )
  )
  rownames(int_matrix) <- .data$estacion
  return(int_matrix)
}

#' Create a function to simulate functional data with a given trend, variance, 
#' and correlation structure
#'
#' @param t vector of time points
#' @param var.teor variance of the functional data
#' @param trend.teor trend of the functional data
#' @param corr.teor correlation structure of the functional data
#' @param rho autocorrelation parameter
#'
#' @returns a function that generates simulated functional data
#'
#' @export
#' @examples
func.sim.set <- function(
  t,
  var.teor = 1,
  trend.teor,
  corr.teor,
  rho = 0
) {
  mdata <- length(t)

  sd.teor <- sqrt(as.numeric(var.teor))

  if (rho != 0) {
    corr.teor <- corr.teor *
      sqrt((1 + rho) / (1 - rho))
  }

  C <- svd(t(corr.teor))

  L.corr.teor <- C$u %*% diag(sqrt(C$d))

  func.sim <- function(rep) {
    err.norm <- matrix(
      rnorm(mdata * rep),
      nrow = mdata
    )

    data.err <- L.corr.teor %*% err.norm

    if (rho != 0) {
      data.err[, 1] <-
        data.err[, 1] *
        sqrt((1 - rho) / (1 + rho))

      for (i in 2:rep) {
        data.err[, i] <-
          rho * data.err[, i - 1] + (1 - rho) * data.err[, i]
      }
    }

    res <- as.numeric(trend.teor) + data.err

    return(res)
  }

  return(func.sim)
}

#' Plot the depth chart and the FDA chart
#'
#' @param res list with the results of the FDA analysis
#'
#' @returns a ggplot object with the depth chart and the FDA chart
#'
#' @export
#' @examples
plot_depth <- function(res) {
  df <- as.data.frame(res$fdata$data)

  df$Curve <- factor(seq_len(nrow(df)))

  df_long <-
    pivot_longer(
      df,
      -Curve,
      names_to = "Hour",
      values_to = "Value"
    )

  df_long$Hour <- as.numeric(df_long$Hour)

  med <-
    data.frame(
      Hour = res$fdata$argvals,
      Value = as.numeric(res$fmed$data)
    )

  env <-
    data.frame(
      Hour = res$fdata$argvals,
      Lower = as.numeric(res$fmin$data),
      Upper = as.numeric(res$fmax$data)
    )

  if (length(res$out) > 0) {
    out_df <-
      df_long %>%
      filter(Curve %in% res$out)
  } else {
    out_df <- NULL
  }

  p1 <-
    ggplot() +

    geom_line(
      data = df_long,
      aes(Hour, Value, group = Curve),
      colour = "grey70",
      linewidth = .4
    ) +

    geom_line(
      data = env,
      aes(Hour, Lower),
      colour = "red",
      linetype = 2,
      linewidth = .8
    ) +

    geom_line(
      data = env,
      aes(Hour, Upper),
      colour = "red",
      linetype = 2,
      linewidth = .8
    ) +

    geom_line(
      data = med,
      aes(Hour, Value),
      colour = "blue",
      linewidth = 1
    ) +

    {
      if (length(res$out) > 0) {
        geom_line(
          data = out_df,
          aes(Hour, Value, group = Curve),
          colour = "black",
          linetype = 3,
          linewidth = .8
        )
      }
    } +

    labs(
      title = "Phase I: FDA Chart",
      x = expression(t),
      y = expression(X(t))
    ) +

    theme_bw() +

    theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none"
    )
  p1 <- p1 +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(10, 10, 10, 10)
    )

  p1 <- p1 +

    # Calibration
    annotate(
      "segment",
      x = 2,
      xend = 3,
      y = -1,
      yend = -1,
      colour = "grey70",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 3.3,
      y = -1,
      label = "Calibration",
      hjust = 0,
      size = 3
    ) +

    annotate(
      "segment",
      x = 8,
      xend = 9,
      y = -1,
      yend = -1,
      colour = "blue",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = 9.3,
      y = -1,
      label = "Median (deepest)",
      hjust = 0,
      size = 3
    ) +

    # Envelope
    annotate(
      "segment",
      x = 16,
      xend = 17,
      y = -1,
      yend = -1,
      colour = "red",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = 17.3,
      y = -1,
      label = "Envelope 99%",
      hjust = 0,
      size = 3
    )

  depth <-
    data.frame(
      Curve = seq_along(res$Depth),

      Depth = res$Depth
    )

  p2 <-
    ggplot(
      depth,
      aes(Curve, Depth)
    ) +

    geom_line() +

    geom_point(size = 2) +

    geom_hline(
      yintercept = res$LCL,
      colour = "red",
      linetype = 2
    ) +

    labs(
      title = "Phase I: Depth Chart",
      x = "Curve index",
      y = "Depth"
    ) +

    theme_bw() +

    theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )

  p1 + p2
}
