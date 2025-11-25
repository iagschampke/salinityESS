## Leer archivo .nc en el mes especificado. Creación data array, vectores lon/lat y level
read_sal_nc <- function(filename, month) {
  nc <- open.nc(filename)
  data <- read.nc(nc)
  close.nc(nc)
  
  dat <- data$salt[,,,month]
  lon0 <- data$lon
  lat0 <- data$lat
  
  coord <- expand.grid(lon0, lat0)
  lon <- coord[,1]
  lat <- coord[,2]
  
  out <- array(dim = c(length(lon0) * length(lat0), length(data$level)))
  for (j in seq_along(data$level)) out[, j] <- c(dat[,, j])
  
  list(array = out, lon = lon, lat = lat, levels = data$level)
}

## Quitar todas las filas de la matriz donde la salinidad sea NA en todos los niveles. También
## se pasa longitud de 0,360 a -180,180
remove_na_rows <- function(df, matrix_data) {
  d <- rowSums(is.finite(matrix_data)) > 0
  matrix_data <- matrix_data[d, , drop = FALSE]
  df <- data.frame(lon = df$lon[d], lat = df$lat[d], matrix_data)
  df$lon <- ifelse(df$lon > 180, df$lon - 360, df$lon)
  df
}

## Se extrae de la matriz de información el área cuadrada especificada. También
## para cada nivel se centra en media 0.
area_filter <- function(xmin, xmax, ymin, ymax, df) {
  df <- dplyr::filter(df, lon >= xmin, lon <= xmax, lat >= ymin, lat <= ymax)
  sal_col_names <- grep("^X", names(df), value = TRUE)
  
  # If no rows, skip scaling entirely
  if (nrow(df) == 0L) {
    return(df)
  }
  
  # If there are no salinity columns, just return what we have
  if (length(sal_col_names) == 0L) {
    return(df)
  }
  
  # Safe scaling: convert scaled matrix to data.frame before assignment
  df[sal_col_names] <- as.data.frame(
    scale(df[sal_col_names], center = TRUE, scale = FALSE)
  )
  
  return(df)
}

## Proyección sinusoidal.
sinusoidal_projection <- function(df, lon_0 = 0, radius = 6378.137) {
  df_proj <- df 
  phi <- df_proj$lat * pi / 180
  lambda <- df_proj$lon * pi / 180
  lambda_0 <- lon_0 * pi / 180
  df_proj$lon <- radius * (lambda - lambda_0) * cos(phi)
  df_proj$lat <- radius * phi
  return(df_proj)
}

## Eliminación de todas las columnas tales que exista al menos un NA
filter_na_depths <- function(df, levels) {
  sal_col_names <- grep("^X", names(df), value = TRUE)
  
  # No salinity columns at all -> nothing to use
  if (length(sal_col_names) == 0L) {
    return(list(df_filtered = df, levels_filtered = numeric(0)))
  }
  
  sal_df <- df[, sal_col_names, drop = FALSE]
  na_col_indexes <- which(sapply(sal_df, function(x) any(is.na(x))))
  
  if (length(na_col_indexes) > 0L) {
    cols_to_remove <- sal_col_names[na_col_indexes]
    df <- df[, !(names(df) %in% cols_to_remove), drop = FALSE]
    levels <- levels[-na_col_indexes]
  }
  
  # If all salinity columns were removed, treat as empty
  remaining_sal <- grep("^X", names(df), value = TRUE)
  if (length(remaining_sal) == 0L) {
    return(list(df_filtered = df, levels_filtered = numeric(0)))
  }
  
  return(list(df_filtered = df, levels_filtered = levels))
}

## Obtención gráficos de salinidad para el área y nivel seleccionados
area_graphs <- function(xmin, xmax, ymin, ymax, df_full, df_zoom, level) {
  lvl_p <- df_full[,level + 2]
  lvl_t <- df_zoom[,level + 2]
  
  p_full <- levelplot(lvl_p ~ df_full$lon * df_full$lat, col.regions = viridis(100))
  map <- p_full + latticeExtra::layer(panel.rect(xleft = xmin, ybottom = ymin, xright = xmax, ytop = ymax, border = "red", lwd = 1), 
                                      data = list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
  zoom <- levelplot(lvl_t ~ df_zoom$lon * df_zoom$lat, col.regions = viridis(100))
  list(map = map, zoom = zoom)
}

## Creación functional data object. Lambda se obtiene como aquel que minimiza el error por GCV
## Se utiliza base de esplines cúbicos
make_fdobj <- function(df, v) {
  sal_col <- grep("^X", names(df))
  sal_mat <- t(as.matrix(df[, sal_col]))
  
  base_spline <- create.bspline.basis(rangeval = range(v), breaks = v, norder = 4)
  
  lambdas <- 10^seq(-10, -5, length.out = 100)
  gcvs <- numeric(length(lambdas))
  
  for (i in seq_along(lambdas)) {
    fdPar_obj <- fdPar(base_spline, Lfdobj = 2, lambda = lambdas[i])
    fit <- smooth.basis(argvals = v, y = sal_mat, fdPar = fdPar_obj)
    gcvs[i] <- mean(fit$gcv)
  }
  lambda <- lambdas[which.min(gcvs)]
  
  fdParObj <- fdPar(base_spline, Lfdobj = 2, lambda = lambda)
  fd_obj <- smooth.basis(argvals = v, y = sal_mat, fdPar = fdParObj)
  list(sal_mat = sal_mat, base_spline = base_spline, fd_obj = fd_obj)
}

## Creación matriz con funciones que interpolan nuestros datos. Sólo es utilizada para
## gráfico
smooth_cols <- function(matrix_data, v) {
  out <- matrix(NA_real_, nrow = nrow(matrix_data), ncol = ncol(matrix_data))
  for (j in seq_len(ncol(matrix_data))) {
    fit <- smooth.spline(x = v, y = matrix_data[, j], spar = 0.5)
    out[, j] <- predict(fit, x = v)$y
  }
  out
}

## Obtención data frame con las anomalías (año a evaluar - media histórica)
anomaly_df <- function(path_a, path_b, month) {
  A <- read_sal_nc(path_a, month = month)
  B <- read_sal_nc(path_b, month = month)
  diff_array <- A$array - B$array
  dfr <- remove_na_rows(data.frame(lon = A$lon, lat = A$lat, diff_array), diff_array)
  list(df = dfr, levels = A$levels)
}

## Vector de profundidad. Se pasa a logaritmo
depths <- function(levels, df_like) {
  k <- ncol(df_like) - 2
  v <- log(as.vector(levels)[1:k])
  return(v)
}

## Obtención del variograma traza empírico.
empirical_tracevariog <- function(coords, fd_obj, base_spline, scale_factor = 1) {
  M <- bsplinepen(base_spline, Lfdobj = 2)
  s <- ncol(fd_obj$fd$coefs)
  L2 <- l2.norm(s = s, datafd = fd_obj$fd, M = M)
  L2_scaled <- L2 * scale_factor
  emp_tv <- trace.variog(coords = coords, L2norm = L2_scaled)
  emp_tv_bin <- trace.variog(coords = coords, L2norm = L2_scaled, bin = TRUE)
  
  list(emp_tv = emp_tv, emp_tv_bin = emp_tv_bin)
}

## Función obsoleta.
init_variog_params <- function(emp_tv_bin) {
  est_nugget <- emp_tv_bin$v[1]
  est_sigma2 <- max(emp_tv_bin$v) - est_nugget
  est_phi <- emp_tv_bin$u[which.max(emp_tv_bin$v)]
  list(nugget = est_nugget, sigma2 = est_sigma2, phi = est_phi)
}

## Parámetros de inicialización ajuste variograma traza. Problema: se puede hacer al ojo pero es bastante
## manual, por lo tanto la idea es obtener posibles aproximaciones para el nugget, rango y partial sill.
## Cortesía de GPT
init_variog_params_unbinned <- function(emp_tv) {
  o <- order(emp_tv$u)
  u <- emp_tv$u[o]; v <- emp_tv$v[o]
  k <- max(5, floor(0.05*length(v)))
  nug <- median(v[1:k])
  sill <- max(quantile(v, 0.95) - nug, .Machine$double.eps)
  target <- nug + 0.95*sill
  phi <- u[which.min(abs(v - target))]
  list(nugget = nug, sigma2 = sill, phi = max(phi, 1e-6))
}

## Ajuste variograma traza.
fit_trace_variog <- function(emp_tv, init_pars, model = "exponential", fix_nugget = FALSE) {
  fit.tracevariog(
    emp.trace.vari = emp_tv,
    sigma2.0 = init_pars$sigma2,
    phi.0    = init_pars$phi,
    nugget   = init_pars$nugget,
    models   = model,
    fix.nugget = fix_nugget)
}

## Obtener curva ajustada, utilizando los parámetros de fit_trace_variog
calculate_fitted_line <- function(model_name, plot_distances, fit_nugget, fit_psill, fit_range, best_model) {
  fit_kappa <- if (!is.null(best_model$kappa)) best_model$kappa else 0.5
  C_h <- cov.spatial(
    obj       = plot_distances,
    cov.model = model_name,
    cov.pars  = c(fit_psill, fit_range),
    kappa     = fit_kappa)
  gamma_h <- fit_nugget + fit_psill - C_h
  gamma_h[plot_distances == 0] <- fit_nugget
  return(gamma_h)
}

## Graficar variograma traza empírico (binned) junto a la curva ajustada. Lamentablemente
## trabajar con la versión binned es complicado así que se hace manualmente
plot_manual_variogram <- function(emp_variogram, fitted_variogram, scale_factor = 1) {
  unscale_div <- scale_factor^2
  
  emp_dist <- emp_variogram$u
  emp_semivar <- emp_variogram$v / unscale_div
  
  best_model <- fitted_variogram$best
  model_name <- best_model$cov.model
  fit_nugget <- best_model$nugget / unscale_div
  fit_psill  <- best_model$cov.pars[1] / unscale_div
  fit_range  <- best_model$cov.pars[2]
  
  max_dist <- max(emp_dist)
  plot_distances <- seq(0, max_dist, len = 100)
  
  fit_line_y <- calculate_fitted_line(model_name, plot_distances, fit_nugget, fit_psill, fit_range, best_model)
  
  all_y_values <- c(emp_semivar, fit_line_y)
  plot_ylim <- c(min(0, min(all_y_values)), max(all_y_values) * 1.05)
  
  plot(emp_dist, emp_semivar, pch = 1, col = "black", ylim = plot_ylim)
  lines(plot_distances, fit_line_y, col = "red", lwd = 2)
  abline(h = fit_nugget, col = "grey", lty = 2)
  abline(h = fit_nugget + fit_psill, col = "blue", lty = 2)
  
  legend("bottomright",
         legend = c("Empirical bins",
                    paste("Fitted", model_name, "model"),
                    "Total sill",
                    "Nugget"),
         col = c("black", "red", "blue", "grey"),
         pch = c(1, NA, NA, NA),
         lty = c(NA, 1, 2, 2),
         lwd = c(NA, 2, 1, 1),
         bty = "n")
}

## Cálculo de matriz trazavariograma estimada y ESS
compute_variance_term <- function(fit_obj, coords_df, scale_factor = 1) {
  
  n <- nrow(coords_df)
  
  unscale_div <- scale_factor^2
  best_model  <- fit_obj$best
  model_name  <- best_model$cov.model
  
  fit_nugget <- best_model$nugget / unscale_div
  fit_psill  <- best_model$cov.pars[1] / unscale_div
  fit_range  <- best_model$cov.pars[2]
  
  sigma_tr_0 <- fit_nugget + fit_psill
  
  coords  <- as.matrix(coords_df[, c("lon", "lat")])
  dist_h  <- as.matrix(dist(coords))
  
  gamma_matrix <- calculate_fitted_line(
    model_name      = model_name,
    plot_distances  = dist_h,
    fit_nugget      = fit_nugget,
    fit_psill       = fit_psill,
    fit_range       = fit_range,
    best_model      = best_model
  )
  
  sigma_matrix <- matrix(sigma_tr_0, nrow = n, ncol = n) - gamma_matrix
  
  denom <- sum(sigma_matrix)
  
  # Guard: if denominator is zero/negative/NaN, return NA instead of Inf
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  
  result <- (n^2 * sigma_tr_0) / denom
  
  return(result)
}

###

generate_square_grid <- function(xmin, xmax, ymin, ymax, dx, dy = dx) {
  xs <- seq(xmin, xmax - dx, by = dx)
  ys <- seq(ymin, ymax - dy, by = dy)
  
  bboxes <- vector("list", length(xs) * length(ys))
  k <- 1L
  for (x in xs) {
    for (y in ys) {
      bboxes[[k]] <- list(
        xmin = x,
        xmax = x + dx,
        ymin = y,
        ymax = y + dy
      )
      k <- k + 1L
    }
  }
  bboxes
}

run_square_pipeline <- function(anom,
                                bbox,
                                level = 1,
                                scale_factor = 1000) {
  dfr <- anom$df
  
  dfr_bbox <- area_filter(bbox$xmin, bbox$xmax,
                          bbox$ymin, bbox$ymax, dfr)
  
  filter_result   <- filter_na_depths(dfr_bbox, anom$levels)
  dfr_bbox_clean  <- filter_result$df_filtered
  levels_clean    <- filter_result$levels_filtered
  
  # Skip if there are no levels or no salinity columns left
  if (length(levels_clean) == 0L) return(NULL)
  if (nrow(dfr_bbox_clean) == 0L) return(NULL)
  
  dfr_proj <- sinusoidal_projection(dfr_bbox_clean)
  v        <- depths(levels_clean, dfr_proj)
  obj      <- make_fdobj(dfr_proj, v)
  sal_spline <- smooth_cols(obj$sal_mat, v)
  
  coord <- as.matrix(cbind(dfr_proj$lon, dfr_proj$lat))
  emp   <- empirical_tracevariog(coord, obj$fd_obj, obj$base_spline,
                                 scale_factor = scale_factor)
  
  init <- init_variog_params_unbinned(emp$emp_tv)
  fit  <- fit_trace_variog(emp$emp_tv, init,
                           model = "exponential",
                           fix_nugget = FALSE)
  
  variance_term <- compute_variance_term(fit_obj = fit,
                                         coords_df = dfr_proj,
                                         scale_factor = scale_factor)
  
  g <- area_graphs(bbox$xmin, bbox$xmax,
                   bbox$ymin, bbox$ymax,
                   dfr, dfr_bbox_clean, level)
  
  list(
    anomaly_df      = dfr,
    bbox_df         = dfr_bbox_clean,
    projected_df    = dfr_proj,
    depth_axis_v    = v,
    fda             = obj,
    sal_spline      = sal_spline,
    emp_variog      = emp$emp_tv,
    emp_variog_bin  = emp$emp_tv_bin,
    fit             = fit,
    variance_term   = variance_term,
    plots           = g,
    scale_factor    = scale_factor,
    bbox            = bbox
  )
}

search_squares_over_months <- function(path_2024 = "salt.2024.nc",
                                       path_ltm  = "salt.mon.ltm.1991-2020.nc",
                                       months    = 1:12,
                                       global_bbox,
                                       level         = 1,
                                       scale_factor  = 1000,
                                       min_locs      = 150,
                                       max_locs      = 500,
                                       target_var    = 1,
                                       dx,
                                       dy = dx) {
  results <- list()
  
  # All candidate boxes in the big region
  candidate_bboxes <- generate_square_grid(global_bbox$xmin,
                                           global_bbox$xmax,
                                           global_bbox$ymin,
                                           global_bbox$ymax,
                                           dx = dx, dy = dy)
  
  for (m in months) {
    message("Processing month ", m, " ...")
    
    # Anomalies once per month
    anom <- anomaly_df(path_2024, path_ltm, month = m)
    month_result <- NULL
    
    for (bb in candidate_bboxes) {
      # Quick location count with NA filtering
      dfr_bbox <- area_filter(bb$xmin, bb$xmax, bb$ymin, bb$ymax, anom$df)
      filter_result  <- filter_na_depths(dfr_bbox, anom$levels)
      dfr_bbox_clean <- filter_result$df_filtered
      levels_clean   <- filter_result$levels_filtered
      
      # Skip if no salinity columns/levels remain
      if (length(levels_clean) == 0L) next
      
      n_loc <- nrow(dfr_bbox_clean)
      if (n_loc < min_locs || n_loc > max_locs) next

      
      # Full pipeline for this square
      r <- run_square_pipeline(anom = anom,
                               bbox = bb,
                               level = level,
                               scale_factor = scale_factor)
      if (is.null(r)) next
      
      if (is.finite(r$variance_term) && r$variance_term >= target_var) {
        month_result <- c(
          list(month = m,
               n_locations = n_loc),
          r
        )
        message("  month ", m,
                ": found square with variance_term = ",
                round(r$variance_term, 2),
                " and n_locations = ", n_loc)
        break
      }
    }
    
    if (is.null(month_result)) {
      message("  month ", m,
              ": no square reached the target variance term.")
    }
    
    results[[paste0("month_", m)]] <- month_result
  }
  
  results
}


