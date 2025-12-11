library(geofd)
library(spatial)
library(lattice)
library(viridis)
library(latticeExtra)
library(dplyr)
library(ggplot2)
library(fda)
library(geoR)
library(geoFourierFDA)
library(skimr)
library(rainbow)
library(fdagstat)
library(profvis)
library(ncdf4)
library(parallel)
library(pbapply)

source("salinity_functions.R")

############### PIPELINE (OBSOLETO POR AHORA) ############### 

pipeline_salinidad <- function(path_2024 = "salt.2024.nc", 
                               path_ltm = "salt.mon.ltm.1991-2020.nc", 
                               month = 1, 
                               bbox = list(xmin, xmax, ymin, ymax), 
                               level = 1) {
  
  anom <- anomaly_df(path_2024, path_ltm, month)
  dfr  <- anom$df
  
  dfr_bbox <- area_filter(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax, dfr)
  g <- area_graphs(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax, dfr, dfr_bbox, level)
  filter_result <- filter_na_depths(dfr_bbox, anom$levels)
  
  dfr_proj <- sinusoidal_projection(filter_result$df_filtered)
  levels <- filter_result$levels_filtered
  
  v <- depths(levels, dfr_proj)
  obj <- make_fdobj(dfr_proj, v)
  coords <- as.matrix(cbind(dfr_proj$lon, dfr_proj$lat))
  
  list(anomaly_df      = dfr,
       bbox_df         = dfr_bbox,
       projected_df    = dfr_proj,
       depth_axis_v    = v,
       fda             = obj,
       coords          = coords,
       plots           = g)
  }

############### INICIALIZACIÓN ###############
path_2024 = "salt.2024.nc"
path_ltm = "salt.mon.ltm.1991-2020.nc"
iter = 100 # cuidado con la RAM
min_side = 5
max_side = 30
min_ub = 200
max_ub = 300
min_depths = 35

nc <- nc_open("salt.2024.nc")
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

# uso de extra cores, PUEDE SER QUE NO FUNCIONE EN MAC!!!!!!!! voy a ver si puedo
# implementarlo

n_cores <- detectCores() - 1L
cl <- makeCluster(n_cores)
pkgs <- grep("^package:", search(), value = TRUE)
pkgs <- sub("^package:", "", pkgs)
clusterExport(cl, ls())
clusterEvalQ(cl, {
  lapply(pkgs, require, character.only = TRUE)
  NULL
})

res <- pblapply(1:iter, cl = cl, function(i) {
  repeat {
    region <- random_region(lon, lat, min_side, max_side)
    n_obs  <- count_cells(region, lon, lat)
    
    if (!(min_ub < n_obs && n_obs < max_ub)) next
    
    anom <- anomaly_df(path_2024, path_ltm, month, region$xmin, region$xmax, region$ymin, region$ymax)
    
    if (nrow(anom$df) == 0) next
    
    filt_df  <- filter_na_depths(anom$df, anom$levels)
    df_f     <- filt_df$df_filtered
    levels_f <- filt_df$levels_filtered

    if (nrow(df_f) == 0 || length(levels_f) <= min_depths) next
    
    df_f$n_obs <- rep(n_obs, nrow(df_f))

    return(list(i = i, region = region, n_obs = n_obs, df = df_f, levels = levels_f))
  }
})


stopCluster(cl)



#########################  BOXPLOTS PRESENTACIÓN

??boxplot
boxplot(res$fda$fd_obj, method = "MBD", prob = c(0.75,0.5,0.25), color = c(1,2,3))

fd_smooth <- res$fda$fd_obj
fd <- fd_smooth$fd

grid <- seq(min(res$depth_axis_v), max(res$depth_axis_v), length.out = 101)

Y <- eval.fd(grid, fd)

fds_obj <- fds(x = grid, y = Y, xname = "log-profundidad", yname = "Anomalías de salinidad")

?fboxplot

fboxplot(data = fds_obj, plot.type = "functional", type = "bag", projmethod = "PCAproj", plotlegend = FALSE, factor = 2)

################ PRUEBA GEOFOURIERFDA
matriz <- res$projected_df[,3:ncol(res$projected_df)]
coords <- res$projected_df[,1:2]
coords <- coords %>% relocate(lat, .before = lon)

??skimr
?geo_model

data(canada)


skim(matriz)
modelo <- geo_model(as.vector(matriz), as.matrix(coords))
################### PRUEBA FDAGSTAT

??fdagstat ### ??????????
data("3ParamFWPR")

Coordinates <- data.frame(lat = coords$lat, lon = coords$lon)
Functions <- data.frame(t(matriz))

ArgStep <- mean(diff(res$depth_axis_v))

D <- as.matrix(dist(Coordinates[, 1:2]))
max_dist <- max(D) / 2

g <- fstat(NULL, vName = "Salinity", Coordinates = Coordinates, Functions = Functions, scalar = FALSE)
g <- estimateDrift("~.", g, Intercept = TRUE)
g <- fvariogram("~.", g, Nlags = 50, LagMax = max_dist, ArgStep = ArgStep, useResidual = FALSE, comments=TRUE)
plotVariogram(g)
head(g$variogram)

g$variogram

emp <- g$variogram

nugget0 <- min(emp$gamma)
sill0 <- mean(emp$gamma[emp$dist > 500]) - nugget0
range0 <- max(emp$dist) / 3

m0 <- vgm(psill = sill0, model = "Exp", range = range0, nugget = nugget0)

m0

g <- fitVariograms(g, model = m0, fitSills = TRUE, fitRanges = TRUE)
plotVariogram(g)
g$model

g <- addCovariance(g, type = 'omni')
str(g, 1)

model <- g$model$omni$Salinity

nugget <- model$psill
psill  <- model$psill[2]
sill   <- sum(model$psill)
range  <- model$range[2]

nugget

