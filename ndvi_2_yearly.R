library(terra)


ms <- vect("data/landscape_data/mississippi_ACEA.shp") 

ndvi_l <- list.files(path="data/ndvi/", pattern=".tif$", full.names=TRUE)
ndvi_w <- rast(ndvi_l[1:5]) 
ndvi_s <- rast(ndvi_l[6:10])

ndvi_w <- project(ndvi_w, crs(ms))
ndvi_s <- project(ndvi_s, crs(ms))

ndvi_w <- mask(ndvi_w, ms)
ndvi_s <- mask(ndvi_s, ms)

writeRaster(ndvi_w, "data/landscape_data/ndvi5year_jan1.tif")
writeRaster(ndvi_s, "data/landscape_data/ndvi5year_jun10.tif")
