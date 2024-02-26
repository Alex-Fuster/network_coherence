# Script to get a diet matrix from GloBI for the Newfoundland shelf groundfish 
# community

# groundfish-data-analysis contains the data folder that is copied from a clone of:
# https://github.com/eric-pedersen/groundfish-data-analysis
load("groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

# get species names from the abundance matrix
taxon_names <- colnames(Year_Geom_Means_all)
taxon_names <- gsub("_", " ", taxon_names)
taxon_names <- stringr::str_to_sentence(taxon_names)

taxon_names_fishbase = taxon_names
# misspellings
taxon_names_fishbase = gsub("Lycodes esmarki", "Lycodes esmarkii", taxon_names_fishbase)
taxon_names_fishbase = gsub("Nezumia bairdi", "Nezumia bairdii", taxon_names_fishbase)
# synonyms
# https://fishbase.mnhn.fr/Nomenclature/SynonymSummary.php?ID=167362&GSID=57841&GenusName=Notocanthus&SpeciesName=nasus&SpecCode=2661&SynonymsRef=92135
taxon_names_fishbase = gsub("Notacanthus nasus", "Notacanthus chemnitzii", taxon_names_fishbase)
taxon_names_fishbase = gsub("Raja spinicauda", "Bathyraja spinicauda", taxon_names_fishbase)
taxon_names_fishbase = gsub("Sebastes marinus", "Sebastes norvegicus", taxon_names_fishbase)
taxon_names_fishbase = gsub("Urophycis chesteri", "Phycis chesteri", taxon_names_fishbase)
taxon_names_fishbase = gsub("Raja radiata", "Amblyraja radiata", taxon_names_fishbase)
# taxon_names_fishbase = gsub("Raja senta", "Malacoraja senta", taxon_names_fishbase)


# get interaction matrix from GloBI
interactions <- rglobi::get_interaction_matrix(source.taxon.names = taxon_names_fishbase,
                                               target.taxon.names = taxon_names_fishbase,
                                               interaction.type = "eats")




######################################


# Further search effort for those species missing interactions

# Raja senta <- Enchelyopus cimbrius
#(reference: https://repository.library.noaa.gov/view/noaa/3755)

# Limanda ferruginea  <- Clupea harengus
#(reference: https://repository.library.noaa.gov/view/noaa/3755)

# Lycodes -> Gadus morhua, Anarhichas lupus
#(reference: GLOBI)

# Raja spinicauda <- sebastes metella
#(reference: https://www.nafo.int/Portals/0/PDFs/sc/2002/scr02-093.pdf)


######################################


# format into a matrix for saving
diet <- lapply(interactions, unlist) 
diet <- do.call(cbind, diet)
rownames(diet) = diet[,1]
diet <- diet[,-1]

write.csv(diet, "b_data/clean/globi_diet_matrix.csv", row.names = TRUE)

# prepare the adjacency matrix
X <- diet
X <- as.matrix(X)
diag(X) <- 0

# plot the network
library(DiagrammeR)
g <- from_adj_matrix(X)
render_graph(g, layout = "nicely", width = 600, height = 600)

# export as a figure
# this looks pretty ugly, so nevermind!! I just took a screenshot called globi_foodweb.png
# export_graph(g, file_name = "figures/globi_foodweb.svg", 
#              width = 1000, height = 1000)

# save the adjacency matrix
write.csv(X, "b_data/clean/globi_adjacencymatrix.csv", row.names = TRUE)
