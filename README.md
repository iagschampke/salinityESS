Code utilized for calculating the Effective Sample Size of salinity data.

Data source: https://psl.noaa.gov/data/gridded/data.godas.html

Note: the geofd R package is NOT available to download from current CRAN as it was archived in 2020. To install it, run the following commands:

```r
old_repo <- c(CRAN = "https://packagemanager.posit.co/cran/2020-02-05")
install.packages("geofd", repos = old_repo)
