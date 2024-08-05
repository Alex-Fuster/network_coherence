# Pakcages ----------------------------------------------------------------
source("./a_code/functions/ipak.R")
ipak(c("sf", "terra"))

# Download shapefile ------------------------------------------------------

# dir.create("b_data/shapefile") # Create folder to save shapefile
# download.file("http://terrabrasilis.dpi.inpe.br/download/dataset/mata-atlantica-aux/vector/biome_border.zip", destfile = "b_data/shapefile/ma_shapefile.zip") # Download shapefile
# unzip("b_data/shapefile/ma_shapefile.zip", exdir = "b_data/shapefile") # Unzip shapefile


# Read Atlantic Rainforest shapefile --------------------------------------

# Read in shapefile, union layers (this shp contains multiple polygons) and simplify geometries
ma_shapefile <- terra::simplifyGeom(terra::union(terra::vect("b_data/shapefile/biome_border.shp")))

# Force counter-clockwise polygons
ma_shapefile <- terra::forceCCW(ma_shapefile)

# Get polygons WKT --------------------------------------------------------

# Get WKT of the polygons to use on the GBIF query
# ma_shapefile_wkt <- terra::geom(ma_shapefile, wkt = TRUE)



