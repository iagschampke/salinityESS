library(geofd)
library(geoR)
library(spatial)
library(lattice)
library(viridis)
library(RNetCDF)
library(latticeExtra)
library(dplyr)
library(sf)
library(ggplot2)
library(fda)

source("salinity_functions.R")

############### INICIALIZADOR ############### 

pipeline_salinidad <- function(path_2024 = "salt.2024.nc", 
                               path_ltm = "salt.mon.ltm.1991-2020.nc", 
                               month = 1, 
                               bbox = list(xmin, xmax, ymin, ymax), 
                               level = 1,
                               scale_factor = 1000) {
  
  anom <- anomaly_df(path_2024, path_ltm, month)
  dfr  <- anom$df
  
  dfr_bbox <- area_filter(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax, dfr)
  g <- area_graphs(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax, dfr, dfr_bbox, level)
  filter_result <- filter_na_depths(dfr_bbox, anom$levels)
  
  dfr_proj <- sinusoidal_projection(filter_result$df_filtered)
  levels <- filter_result$levels_filtered
  
  v <- depths(levels, dfr_proj)
  obj <- make_fdobj(dfr_proj, v)
  sal_spline <- smooth_cols(obj$sal_mat, v)
  
  par(mfrow = c(1, 3))
  matplot(v, obj$sal_mat, type = "p", pch = 16, cex = 0.5)
  matplot(v, sal_spline, type = "l", lty = 1)
  plot(obj$fd_obj)
  par(mfrow = c(1, 1))
  
  coord <- as.matrix(cbind(dfr_proj$lon, dfr_proj$lat))
  emp <- empirical_tracevariog(coord, obj$fd_obj, obj$base_spline, scale_factor = scale_factor)
  
  init <- init_variog_params_unbinned(emp$emp_tv)
  fit <- fit_trace_variog(emp$emp_tv, init, model = "exponential", fix_nugget = FALSE)
  
  list(anomaly_df      = dfr,
       bbox_df         = dfr_bbox,
       projected_df    = dfr_proj,
       depth_axis_v    = v,
       fda             = obj,
       sal_spline      = sal_spline,
       emp_variog      = emp$emp_tv,
       emp_variog_bin  = emp$emp_tv_bin,
       fit             = fit,
       plots           = g,
       scale_factor    = scale_factor)
}

############### INICIALIZACIÓN ###############

sf <- 1000
res <- pipeline_salinidad(month = 2, bbox = list(xmin = -135, xmax = -125, ymin = 5, ymax = 15), level = 10, scale_factor = sf)

#res$plots$map
plot_manual_variogram(res$emp_variog_bin, res$fit, scale_factor = res$scale_factor)
res$plots$zoom
res$fit$best
length(res$bbox_df$lon)

variance_term <- compute_ESS(fit_obj = res$fit, coords_df = res$projected_df, 
                                       scale_factor = res$scale_factor)
print(variance_term)



