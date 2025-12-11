## Obtención data frame con las anomalías (año a evaluar - media histórica)
anomaly_df <- function(path_a, path_b, month, xmin, xmax, ymin, ymax) {
  A <- read_sal_nc(path_a, month, xmin, xmax, ymin, ymax)
  B <- read_sal_nc(path_b, month, xmin, xmax, ymin, ymax)
  diff_array <- A$array - B$array
  dfr <- remove_na_rows(data.frame(lon = A$lon, lat = A$lat, diff_array), diff_array)
  list(df = dfr, levels = A$levels)
}

## Leer archivo .nc en el mes especificado. Creación data array, vectores lon/lat y level
read_sal_nc <- function(filename, month, xmin, xmax, ymin, ymax) {
  nc <- nc_open(filename)
  on.exit(nc_close(nc))
  lon0 <- ncvar_get(nc, "lon")
  lat0 <- ncvar_get(nc, "lat")
  lev <- ncvar_get(nc, "level")
  ix <- which(lon0 >= xmin & lon0 <= xmax)
  iy <- which(lat0 >= ymin & lat0 <= ymax)
  nlon <- length(ix)
  nlat <- length(iy)
  nlev <- length(lev)
  dat_sub <- ncvar_get(nc, "salt", start = c(min(ix), min(iy), 1, month), count = c(nlon, nlat, nlev, 1))
  lon_sub <- lon0[ix]
  lat_sub <- lat0[iy]
  coord <- expand.grid(lon_sub, lat_sub)
  out <- matrix(NA_real_, nrow = nlon * nlat, ncol = nlev)
  
  for (k in seq_len(nlev)) {
    out[, k] <- c(dat_sub[, , k])
  }
  
  list(array = out, lon = coord[, 1], lat = coord[, 2], levels = lev)
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
  df[sal_col_names] <- scale(df[sal_col_names], center = TRUE, scale = FALSE)
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
  sal_df <- df[, sal_col_names]
  na_col_indexes <- which(sapply(sal_df, function(x) any(is.na(x))))
  if (length(na_col_indexes) > 0) {
    cols_to_remove <- sal_col_names[na_col_indexes]
    df <- df[, !(names(df) %in% cols_to_remove)]
    levels <- levels[-na_col_indexes]
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

## Vector de profundidad. Se pasa a logaritmo
depths <- function(levels, df_like) {
  k <- ncol(df_like) - 2
  v <- log(as.vector(levels)[1:k])
  return(v)
}


## Cálculo de matriz trazavariograma estimada y ESS
compute_ESS <- function(nugget, psill, range, coords_df) {
  n <- nrow(coords_df)
  model_name <- "exponential"
  sigma_tr_0 <- nugget + psill
  coords <- as.matrix(coords_df[, c("lon", "lat")])
  dist_h <- as.matrix(dist(coords))
  
  gamma_matrix <- calculate_fitted_line(
    model_name     = model_name,
    plot_distances = dist_h,
    fit_nugget     = nugget,
    fit_psill      = psill,
    fit_range      = range)
  
  sigma_matrix <- matrix(sigma_tr_0, nrow = n, ncol = n) - gamma_matrix
  denom <- sum(sigma_matrix)
  result <- ifelse(abs(denom) < 1e-10, NA, (n^2 * sigma_tr_0) / denom)
  return(result)
}

## Obtener curva ajustada, utilizando los parámetros de fit_trace_variog
calculate_fitted_line <- function(model_name, plot_distances, fit_nugget, fit_psill, fit_range) {
  fit_kappa <- 0.5
  C_h <- cov.spatial(obj = plot_distances, cov.model = model_name, cov.pars = c(fit_psill, fit_range), kappa = fit_kappa)
  gamma_h <- fit_nugget + fit_psill - C_h
  gamma_h[plot_distances == 0] <- fit_nugget
  return(gamma_h)
}

## rectángulo aleatorio dentro de los límites del mapa y con lados no mayores a 30

random_region <- function(lon, lat, min_side = 5, max_side = 30) {
  valid_lon <- lon[lon >= 0 & lon <= 360]
  valid_lat <- lat[lat >= -74.5 & lat <= 64.5]
  
  get_bounds <- function(vec, min_s, max_s) {
    for (attempt in 1:100) {
      val1 <- sample(vec, 1)
      candidates <- vec[vec >= (val1 + min_s) & vec <= (val1 + max_s)]
      if (length(candidates) > 0) {
        val2 <- sample(candidates, 1)
        return(c(val1, val2))
      }
    }
    stop("error")
  }
  
  x_coords <- get_bounds(valid_lon, min_side, max_side)
  y_coords <- get_bounds(valid_lat, min_side, max_side)
  
  list(xmin = x_coords[1], xmax = x_coords[2], 
       ymin = y_coords[1], ymax = y_coords[2])
}

##

count_cells <- function(region, lon, lat) {
  ix <- which(lon >= region$xmin & lon <= region$xmax)
  iy <- which(lat >= region$ymin & lat <= region$ymax)
  length(ix) * length(iy)
}