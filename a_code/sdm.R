# Loading packages --------------------------------------------------------

# Use the code bellow in case the Rocc package is not installed

# if (!"devtools" %in% installed.packages()){install.packages("devtools")}  
# devtools::install_github("liibre/Rocc")

# remotes::install_github("sjevelazco/flexsdm")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("readxl","Rocc","terra","sdmpredictors","rgdal","sf","flexsdm", 
              "tidyverse","ENMTML")

ipak(packages); rm(packages,ipak)

# Species occurrences -----------------------------------------------------

ssp <- readxl::read_excel(path = "./SDMs/data/raw/Heike_Traits.xlsx", sheet = 1)

no_coordinates = matrix(nrow = length(ssp), ncol = 1)

for (i in 1:nrow(ssp)) {
  data = rgbif2(species = ssp[i,1], save = FALSE, limit = 99500)
  if(nrow(data) == 0){
    no_coordinates[i,1] = ssp[i]
    no_coordinates = na.omit(as.data.frame(no_coordinates))
    write.table(x = no_coordinates, file = "./no_coordinates.txt", col.names = "species", row.names = FALSE)
  }
  if(nrow(data) !=  0){
    data = data[,c("species_search","decimalLongitude","decimalLatitude")]
    if(i == 1){
      gbif = data
    } else{
      gbif = rbind(gbif, data)
    }
  }
  print(i)
}; rm(data, no_coordinates)


# Cleaning occurrence data ------------------------------------------------

# Removing duplicates

gbif <- unique(gbif)

# Subsetting the records based on the extent of the marine biogeographical realm of the study area
# https://www.nature.com/articles/s41467-017-01121-2

gbif <- subset(gbif, decimalLongitude >= -72.161264 & decimalLongitude <= -1.199208 & decimalLatitude >= 29.549038 & decimalLatitude <= 70.034463)

write.table(gbif,"./SDMs/data/clean/occ_clean.txt", row.names = FALSE)

# Getting environmental data ----------------------------------------------

# Listing data sets available on sdmpredictors package

list_env_layers <- list_layers(datasets=c("Freshwater", "Bio-ORACLE", "MARSPEC"), 
                               terrestrial = FALSE, 
                               marine = TRUE, 
                               freshwater =TRUE)$layer_code

subset_env_layers <- str_subset(list_env_layers, 
                                str_c(c(".temp.", ".sal.", ".oxy."), 
                                      collapse = "|"))

# environ_layers <- load_layers(subset_env_layers, rasterstack = TRUE, datadir = "b_data/environmental")

# Random sample of environmental variables

environ_layers <- load_layers(sample(subset_env_layers, 10, replace = FALSE), 
                              rasterstack = TRUE, 
                              datadir = "b_data/environmental")

# Define a boundary

boundary = extent(c(xmin = -72.161264, xmax = -1.199208, ymin = 29.549038, ymax = 70.034463))

# Cropping environmental variables based on the extent of the marine biogeographical realm of the study area

environ_layers_cropped <- crop(environ_layers, boundary)

# Performing a PCA using flexsdm package

environ_layers_pca <- correct_colinvar(environ_layers_cropped, method = c('pca'), proj = NULL, maxcell = NULL)

# Saving PCs derived from PCA

writeRaster(environ_layers_pca$env_layer$PC1, filename="./SDMs/PCA/PC1.tif")
writeRaster(environ_layers_pca$env_layer$PC2, filename="./SDMs/PCA/PC2.tif")
writeRaster(environ_layers_pca$env_layer$PC3, filename="./SDMs/PCA/PC3.tif")
writeRaster(environ_layers_pca$env_layer$PC4, filename="./SDMs/PCA/PC4.tif")
writeRaster(environ_layers_pca$env_layer$PC5, filename="./SDMs/PCA/PC5.tif")

# Ecological Niche Modeling

ENMTML (pred_dir = "./PCA", 
        occ_file = "./occ_clean.txt",
        sp = "species_search",
        x = "decimalLongitude",
        y = "decimalLatitude",
        min_occ = 10,
        thin_occ = NULL,
        part = c(method = "BLOCK"),
        colin_var = NULL,
        sp_accessible_area = NULL,
        pseudoabs_method = c(method= 'ENV_CONST'),
        pres_abs_ratio = 1,
        algorithm = c("RDF", "MLK", "MXS", "GAU", "SVM"),
        thr = c(type = "SORENSEN"),
        msdm = NULL,#c(method = 'PRES'),
        ensemble = c(method= 'SUP', metric= 'Sorensen'),
        cores = 4)




