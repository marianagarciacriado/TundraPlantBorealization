## Arctic plant borealization
## Mariana Garcia Criado
## Script 2. Distance to treeline & permafrost
## February 2024


## PACKAGES ----
library(tidyverse)
library(brms)


## DATA ----

# previous file with treeline data
subs.data <- read.csv("data/subsites_treeline_completed_dec2024_final.csv") %>%
  select(-c(LAT, LONG, Column1))

# correct lat-long, there is some issue with conversion
# subs.latlong is the final subset of subsites
subs.latlong <- read.csv("data/subsites_treeline_oct2024.csv")

# merge with the correct subsites and coordinates
subs.correct <- subs.latlong %>% left_join(., subs.data, by = "SiteSubsite")


# Distance to treeline
subs.distance <- subs.correct %>% 
  mutate(ElevationDifference = ElevationSubsite - ElevationTreeline) %>%
  mutate(DistanceTreelineCorrected = case_when(biome == 6 ~ sqrt((ElevationDifference)^2 + (DistanceToTreeline)^2),
                                               biome == 11 ~ DistanceToTreeline)) %>%
  mutate(InterruptedCorrected = case_when(Interrupted %in% c("small lakes and glacial river", "many little lakes",
                                                             "little lakes", "small sea and lakes") ~ "Small water bodies",
                                          Interrupted %in% c("mountain ridge", "mountain range") ~ "Mountains",
                                          Interrupted == "ocean" ~ "Large water bodies",
                                          Interrupted == "no" ~ "Uninterrupted"))

write.csv(subs.distance, "data/subs_distance_treeline_corrected_dec2024.csv")


## PERMAFROST ----

# Permafrost data from Obu et al. (2019) can be accessed at 
# https://www.sciencedirect.com/science/article/pii/S0012825218305907

## Data previously extracted by Joe Everest
permafrost <- read.csv("data/output_boreal_itex_permafrost_fromjoe.csv")

# join with ITEX mastersheet
load("data/boreal_master_dec2024.RData")

boreal_master_perm <- itex.full16 %>% left_join(., permafrost, by = c("LAT", "LONG"))

# fill in the data for the 3 extra sites that were not extracted 
# double-checked against the Obu map and the surrounding sites
boreal_master_perm2 <- boreal_master_perm %>% 
  mutate(permafrost = case_when(SITE == "INCLINE_LAV" ~ "none",
                                SITE %in% c("BILLEFJORDEN", "PYRAMIDEN") ~ "continuous",
                                TRUE ~ permafrost))

save(boreal_master_perm2, file = "data/boreal_master_perm.RData")

