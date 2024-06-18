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
