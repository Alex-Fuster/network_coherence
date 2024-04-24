# Pakcages ----------------------------------------------------------------
source("./a_code/functions/ipak.R")
ipak(c("sf", "terra"))

# Download shapefile ------------------------------------------------------

# dir.create("b_data/shapefile")
# download.file("http://terrabrasilis.dpi.inpe.br/download/dataset/mata-atlantica-aux/vector/biome_border.zip", destfile = "b_data/shapefile/ma_shapefile.zip")
# unzip("b_data/shapefile/ma_shapefile.zip", exdir = "b_data/shapefile")


# Read Atlantic Rainforest shapefile --------------------------------------

ma_shapefile <- terra::simplifyGeom(terra::union(terra::vect("b_data/shapefile/biome_border.shp")))
ma_shapefile <- terra::forceCCW(ma_shapefile)


# Get polygons WKT --------------------------------------------------------

# ma_shapefile_wkt <- terra::geom(ma_shapefile, wkt = TRUE)



