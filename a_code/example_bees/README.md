# Code rationale for the pollinators example

This series of scripts reproduces in full our example with plants and pollinators in our paper about ecological network coherence. In this example, we investigate how similarly plants and pollinators from a network in the Atlantic Rainforest respond to specific environmental conditions. The overall workflow is:

1. Get network data as raw data, from which we will get species names and get biome shapefile to restrict analysis
2. Get occurrence data (in our case, GBIF) to map the environmental niche of the species
3. Get the environmental data from the occurrence points
4. Clean the data to remove biases
5. Calculate the derivatives from each species' niche at each spatial point within the biome
6. Calculate the correlations between plants' and pollinators' derivatives

## Network and biome data

These are obtained with scripts 01 and 01.1. The network is a weighted interaction matrix with plants (bromeliads) and pollinators (bees, butterflies and hummingbirds). [Line 29](https://github.com/Alex-Fuster/network_coherence/blob/0d7fa6f4f10db8d9680dab9a27b4ebf050971e14/a_code/example_bees/01-get_interaction_data.R#L29) removes species identified only up to Genus, and we did that because of the immense diversity of these groups. Keeping them really messed up the rest of the analyses.

We chose to delimit the analyses to the Atlantic Rainforest because this was the biome where the interactions were sampled. It makes ecological sense to keep the extrapolation of the interactions within a biome because this reduces the variation that could influence our ability to assume that an interaction would occur at each co-occurrence point since a particular combination of environmental conditions characterizes a biome.

## Occurrence data

The next step was to get occurrence data for each species within the Atlantic Rainforest. We used GBIF as a first exploration and we like how this raises important concerns about data gaps. GBIF is one of the most used biodiversity databases, but still, it lacks a lot of local occurrences especially in "Global South" countries like Brazil. We could fill these gaps by getting occurrences from experts, from the literature, and from Brazilian biodiversity data sources, but we chose to keep it simple for the sake of demonstration.

We detail how we cleaned this data below.

## Environmental data

The script 03 gets the environmental data for each cell. We use the `sdmpredictors` package and select only the terrestrial variables (because our species are terrestrial) related to mean temperature, for no specific reason. Then we randomly retrieve one of them (because we don't need any specific environmental variable to demonstrate that the method works) and crop the raster using the biome shapefile.

Then we get environmental values for each cell where there is an occurrence of plants, and another dataset is created with environmental values for each cell where there is an occurrence of pollinators, keeping the cell ID, which we'll use later.

## Data cleaning

With script 04 we get rid of duplicated occurrences. We decided to keep only one record per raster cell to avoid bias inflation. Also, we filtered out species with only one occurrence because they wouldn't give us information on the slope or curve of environmental tolerances, which are needed to calculate the network coherence. Finally, we excluded from the interaction matrix the species with no occurrences left.

The choice to keep unique occurrences comes from the rationale that the abundance of GBIF occurrences is rarely a reflection of the species' preferences. By keeping unique occurrences only we still get the minimum temperature, the maximum temperature, and the temperature at which we find the species more frequently. Of course, this is masked by the distribution of temperature of our raster: if our species has one occurrence in each cell of our raster, its temperature curve will be identical to the temperature curve of the raster. We don't think this is an ecological infraction as both represent the truth for that sampling area, and this correlation could be easily dissolved with the addition of more variables.

Moreover, the removal of species with only one occurrence led to results where the network coherence is null because the interacting species don't co-occur.

## Derivatives

The derivatives were then calculated on script 05. First, we calculate the function of the curve using the maximum, minimum, mean and standard deviation of temperatures for each species in the whole Atlantic Rainforest. Then, we create a function to calculate the derivative of this function at each cell of our raster and its environmental value. We then loop it for all species and the cells where they occur.

## Correlations

Finally, correlations between derivatives are calculated on script 06. The loop includes lines to calculate the network coherence for interacting and non-interacting species. We can uncomment lines 15 and 22-24 to calculate it only for interacting species.

## Figures

Figures are all created on scripts starting with "99-".
