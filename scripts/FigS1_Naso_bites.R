## Script to produce Figure S1 that shows percentage of bites Naso lituratus took on two categories of algae: 1) mature macroalgae and 2) immature macroalgae and turf

# Data was collected during 20 minute follows of foraging Naso lituratus (n=19 individuals) on the north shore lagoon of Moorea. Data describe the number of bites a fish took on 5 categories of algae: turf and immature macroalgae, Sargassum, Turbinaria, Dictyota, and Amansia.


# Packages -----
library(tidyverse)
library(ggplot2)
library(readr)

# Data -----
# Fish feeding data
bite.data <- read_csv("data/Naso_lituratus_bite_data_Moorea_2017.csv") %>% 
  select(fish_ID, species, turf_bites, sargassum_bites, turbinaria_bites, dictyota_bites, amansia_bites)

# Metadata to group algal taxa 
algae_group <- read_csv("data/Naso_lituratus_bite_data_Moorea_2017_algae_groups.csv")

# Data wrangling -----
# Get df into format to calculate total bites for each algal taxa (pool fish)
naso.data <- bite.data %>% 
  pivot_longer(cols = ends_with("bites"),
               names_to = "taxa",
               values_to = "bites") %>% 
  mutate_at("taxa", str_replace, "_bites", "")

# Add algae group
naso.data <- merge(naso.data, algae_group, by = "taxa")

# Calculate total bites per algal taxa
naso.data.summary <- naso.data %>% 
  group_by(species, algae_group) %>% 
  summarize(bites = sum(bites)) %>% 
  mutate(total_bites = sum(bites)) %>% 
  mutate(percent_bites = bites/total_bites*100)

# Visualization: percent of total bites taken by Naso lituratus -----
ggplot(naso.data.summary, 
       aes(x = species, 
           y = percent_bites,
           fill = algae_group)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("") +
  ylab("Percent of total bites") +
  scale_fill_discrete(name = "Algal category") +
  theme_classic() +
  theme(text = element_text(size = 12))

