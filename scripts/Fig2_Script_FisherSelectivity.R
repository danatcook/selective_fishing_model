### Script and data to produce figure 2 that shows theoretical selectivity relationships and actual fisher selectivity of herbivores in Moorea ###
# Data used for fisher selectivity of herbivores in Moorea is from Figure 5 of
# Rassweiler, A., Lauer, M., Lester, S.E. et al. Perceptions and responses of Pacific Island fishers to changing coral reefs. Ambio 49, 130–143 (2020). https://doi.org/10.1007/s13280-019-01154-5


# Packages -----
library(tidyverse)
library(ggplot2)
library(readr)

# Assign colors -----
Pcol <- '#fd8d3c' # Parrotfish (Chlorurus/Scarus) or non-browsing herbivores color
Ucol <- '#99d8c9' # Unicornfish (Naso) or browsing herbivores color

# Visualization: theoretical selectivity curves -----

## Build catch-abundance relationship
f <- 0.2 # Set total fishing effort
prefset <- c(.1,.3,.5,.7,.9) # Set of fishing preferences for parrotfish

# Make a plot to hold the data we'll create in the next step
par(pty="s") # sets aspect ratio of plot = 1 (i.e., makes plot a square)
plot(c(0,1),c(0,1), type='n', las=1, xaxt = "n", yaxt = "n", xlab='Relative fishable biomass on the reef', ylab='Relative biomass in the catch') 
axis(side=1, at=c(0, 0.25, 0.5, 0.75, 1))
axis(side=2, at=c(0, 0.25, 0.5, 0.75, 1))

TFP <- 100 # Set total fish population
Uset <- seq(from = 0, to = TFP, by = 1) # Full theoretical range of unicornfish population
Pset <- TFP - Uset # Make parrotfish population 1-Uset
propU <- Uset/(Uset + Pset) # Proportion of unicornfish in the population
for(j in 1:length(prefset)){ # Loop over the different choices of preferences
  p <- prefset[j] # Choose the preference for parrotfish
  u <- 1 - p # Compute the corresponding preference for unicornfish
  propUcatch <-(1-p)*f*Uset/((1-p)*f*Uset+p*f*Pset) # Compute the proportion of catch comprised of unicornfish
  lines(propU,propUcatch,lwd=5*u) # Add a line to the plot showing this
}


# Visualization: Fisher selectivity in Moorea -----
# Here we're replotting reef abundance and catch data from Rassweiler et al. 2020 Figure 5

## Load data
# Data: abundance of fish on reef and in the catch in Moorea, taken from Rassweiler et al. 2020, Ambio, Figure 5
fish <- read_csv("data/Rassweiler_2020_Fig5_data.csv")

# Set order of taxa groups for plot
fish$Taxa <- factor(fish$Taxa, levels = c("Naso", "Scaridae", "Acanthurus/Ctenochaetus"))

# Replot Rassweiler et al. 2020 reef abundance-catch data 
ggplot(fish, aes(x=Reef_Biomass, y=Catch_Biomass, color=Taxa)) +
  geom_point(size=5) +
  xlab("Relative biomass on the reef") +
  ylab("Relative biomass in the catch") + 
  scale_color_manual(labels=c("Naso", "Scaridae", "Acanthurus/Ctenochaetus"), values = c(Ucol, Pcol, 'grey0')) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size=15, colour="black"),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black")) +
  theme(axis.text = element_text(size = 15)) +
  #geom_abline(intercept = 0, slope = 1) + # adds 1:1 line
  xlim(0,1) +
  ylim(0,1) +
  theme(aspect.ratio=1)  # sets aspect ratio = 1 (i.e., makes plot a square)


