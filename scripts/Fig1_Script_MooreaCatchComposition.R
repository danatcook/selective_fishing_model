### Script and data to produce figure 1 that shows composition of the catch in Moorea ###
# Dataset: SEES Roadside Catch in 2014 and 2015
# Data describes the identity and size of fish obtained from photographs of roadside stands around the island

# Packages -----
library(tidyverse)
library(readr)
library(ggplot2)

# Assign colors -----
Hercol <- 'green' # herbivore color
Carcol <- 'red' # carnivore color
Placol <- 'blue' # planktivore color
Pcol <- '#fd8d3c' # non-browsing herbivore color (parrotfish)
Ucol <- '#99d8c9' # browsing herbivore color (unicornfish)

# Data -----
SEES <- read_csv("data/SEES_roadside_data_2014_2015_bcodmo_dataset_709963_712b_5843_9069.csv") 
SEES <- SEES[-1,]

# List of functional groups for fish in the catch
func_SEES <- read_csv("data/trophic_and_functional_groups_SEES_roadeside_data_2014_2015_bcodmo_dataset_709963_712b_5843_9069.csv")



# Data wrangling: proportional abundance of trophic and functional groups in catch -----
# Extract year from date column, currently in format 'mm/dd/yyyy' and we want 'yyyy'
SEES$date <- as.POSIXct(SEES$date, format = "%Y/%m/%d")
SEES$date <- format(SEES$date, format="%Y")

## Summarize SEES dataset: calculate total and proportional abundance of 'browsers', 'scrapers/excavators', and 'all other fish' per year
# Subset dataset to only keep taxonomy and length columns
SEES <- SEES %>% 
  select(genus,length, date) 
# Merge data with functional group metadata
SEES <- merge(SEES,func_SEES, by = "genus")
# Sum the number of fish in each group
SEES <- SEES %>% 
  mutate(Number = 1)

# Drop fish that we can't assign to a functional group (only 3 observations)
SEES_totalYear <- SEES %>% 
  group_by(Trophic_Level, Functional_Group) %>% 
  summarise(Number = sum(Number)) %>% 
  filter(Functional_Group != "Labridae") %>%  # Drop 'Labridae' 
  filter(Functional_Group != "Unknown") %>%  # Drop 'Unknown'
  filter(Functional_Group != "Other") %>%   # Drop 'Other'
  ungroup()

# Calculate total number of fish caught and proportion of each group in the catch for both years
SEES_totalYear <- SEES_totalYear %>%
  mutate(Total_Number = sum(Number)) %>% 
  mutate(Prop_Number = Number/Total_Number)



# Visualization: proportional abundance of trophic and functional groups in catch -----

# Plot the total proportional abundances of each group in the SEES catch (total catch for 2014 and 2015)
# Grouped barplot
ggplot(SEES_totalYear,
       aes(x = factor(Trophic_Level, level = c('Herbivore','Carnivore','Planktivore')), y=Prop_Number, fill=Functional_Group))+
  geom_bar(position="stack", stat="identity") + 
  xlab("") +
  ylab("Proportional abundance in catch\n") +
  scale_fill_manual(values=c(Ucol, Pcol, 'grey0', 'grey0'), 
                    name="Functional Group",
                    breaks=c("Browsing herbivore", "Non-browsing herbivore", "", ""),
                    labels=c("Browsing herbivores", "Non-browsing herbivores", "", "")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size=15, colour="black"),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black")) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0))

