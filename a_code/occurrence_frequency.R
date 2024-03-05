
# Packages loading --------------------------------------------------------

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("tidyverse", "terra", "sdmpredictors")

ipak(packages); rm(packages,ipak)




# Load species data and splitting -----------------------------------------

ssp_cleaned <- read.table("b_data/Occurrences_Cleaned.txt", header = TRUE)

species_names <- unique(ssp_cleaned$sp)

ssp_01 <- filter(ssp_cleaned, sp == species_names[1])

ssp_02 <- filter(ssp_cleaned, sp == species_names[2])

# Getting environmental data ----------------------------------------------

# Listing data sets available on sdmpredictors package

list_env_layers <- list_layers(datasets=c("Freshwater", "Bio-ORACLE", "MARSPEC"), 
                               terrestrial = FALSE, 
                               marine = TRUE, 
                               freshwater =TRUE)$layer_code

subset_env_layers <- str_subset(list_env_layers, 
                                str_c(c(".temp.", ".sal.", ".oxy."), 
                                      collapse = "|"))

# Random sample of environmental variables

environ_layers <- load_layers(sample(subset_env_layers, 1, replace = FALSE), 
                              rasterstack = TRUE, 
                              datadir = "b_data/environmental")

# Define a boundary

boundary = c(xmin = -72.161264, xmax = -1.199208, ymin = 29.549038, ymax = 70.034463)

# Cropping environmental variables based on the extent of the marine biogeographical realm of the study area

environ_var1 <- crop(environ_layers, boundary)

rm(boundary, environ_layers)

# Get environ values for occurrence points --------------------------------

environ_values_sp1 <- extract(environ_var1, ssp_01[2:3])

environ_values_sp2 <- extract(environ_var1, ssp_02[2:3])

combined_df <- rbind(cbind(ssp_01, environ_values_sp1), 
                     cbind(ssp_02, environ_values_sp2))


# Get environmental curve -------------------------------------------------

ggplot(combined_df, aes(x=BO21_salinityltmax_ss, color = sp))+
  geom_density()


# Trying to get derivatives -----------------------------------------------


