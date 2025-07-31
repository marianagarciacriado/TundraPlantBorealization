## Arctic plant borealization
## Mariana Garcia Criado 
## Script 5. Trait data
## May 2024

# The raw data files with traits are too large to upload to Zenodo 
# TRY data (v6.0) can be downloaded at https://www.try-db.org/TryWeb/Home.php
# The cleaned up version of the data ready for analysis is produced at the end of this script ('species_traits_dataset_full.csv')


## LIBRARIES ----
library(arrow)
library(tidyverse)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## TRAIT DATA ----

# db with traits and coordinates (too large and not uploaded to Zenodo/GitHub)
try6 <- read_parquet("M:/# Manuscripts in progress/# Borealization/trait data/TRY_6_coords.parquet")

unique(try6$TraitName)
unique(sort(try6$Dataset))

# species list for filtering - use the one with synonyms
# to ensure that we get trait data for as many species as possible

# species list for filtering - this db has only named species (no morphos)
sp.list.ab <- read.csv("data/species_ab_class_jan2025.csv") #287
sp.list.vector <- sp.list.ab$SPECIES_CLEAN

# this db has some synonyms from script 05-caff-boreal-categories
sp.list.syno <- read.csv("data/species_list_with_synos.csv")
sp.list.syno.vector <- sp.list.syno$Synonym

# Standardise trait names
try6.0 <- try6 %>% 
  mutate(TraitNameNew = case_when(TraitName %in% c("Plant height", "Plant height vegetative", "Plant height generative") ~ "PlantHeight",
                                  TraitName %in% c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included", 
                                                   "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded") ~ "SLA",
                                  TraitName == "Seed dry mass" ~ "SeedMass",
                                  TraitName == "Leaf nitrogen (N) content per leaf dry mass" ~ "LeafN",
                                  TraitName == "Leaf carbon/nitrogen (C/N) ratio" ~ "LeafCN")) %>%
  filter(LAT_site > 50 | is.na(LAT_site)) %>%
  filter(AccSpeciesName %in% c(sp.list.vector, sp.list.syno.vector))

# Can we figure out if trait > 50deg latitude from database name for those with no coordinate data?
na.db <- try6.0 %>% filter(is.na(LAT_site)) %>% distinct(Dataset)

## INSIDE RANGE
# Ecological Flora of the British Isles
# PLANTATT - Attributes of British and Irish Plants
# Sheffield Database
# The Netherlands Plant Height Database

## OUTSIDE RANGE
# BASECO: a floristic and ecological database of Mediterranean 
# BROT 2.0 (Med)
# Caucasus Plant Traits Database
# Cedar Creek Savanna SLA, C, N Database
# DISEQU-ALP (Alps)
# Freschet et al. 2015 - Mount Hutt (New Zealand)
# Global A, N, P, SLA Database
# Global Respiration Database
# Global Seed Mass, Plant Height Database
# Italian Alps Plant Traits Database
# Leaf trait records of rare and endangered plant species in the Pannonian flora
# LECA - Traits of the European Alpine Flora (Europ Alps)
# Plant Traits from Romania
# Species and trait shifts in Apennine grasslands
# SW Michigan restored prairies
# Trait Data from Niwot Ridge LTER (2016)
# Traits for Herbaceous Species from Andorra

## UNCLEAR
# 1000Seedweight
# Baccara - Plant Traits of European Forests
# BiolFlor Database - Germany only but parts of Germany are >50
# BIOME-BGC Parameterization Database
# BIOPOP: Functional Traits for Nature Conservation - Germany
# Chinese Leaf Traits Database - from boreal to 
# Cold Tolerance, Seed Size and Height of North American Forest Tree Species - big range from N to Florida
# Dispersal Traits Database
# Functional Flowering Plant Traits
# Functional traits explaining variation in plant life history strategies
# KEW Seed Information Database (SID)
# Leaf N-Retention Database
# Leaf, trap, and seed traits of carnivorous plants
# Leaf Physiology Database
# Linking hard and soft traits
# Maximum Height of Chinese Tree Species (From Silva Sinica)
# Midwestern and Southern US Herbaceous Species Trait Database
# Nutrient Resorption Efficiency Database
# Onoda 2017 leaf dataset
# Pladias: Life forms and heights of the Czech flora (czechia ><50)
# Plant Physiology Database
# PLANTSdata USDA
# Reich-Oleksyn Global Leaf N, P Database
# Roots Of the World (ROW) Database
# Seed Characteristics of Ericaceae - from Spain to Sweden
# Seed Information Database (SID) Seed Mass 2010
# Seed Information Database, Royal Botanic Gardens, Kew
# Seed Mass from Literature
# Sheffield & Spain Woody Database
# The Americas N&P database
# The DIRECT Plant Trait Database
# The Functional Ecology of Trees (FET) Database  - Jena
# The LEDA Traitbase (some regions <50)
# The Tansley Review LMA Database
# The Xylem/Phloem Database
# Traits from the Wildfire Project
# Traits of 59 grassland species
# Tundra Plant Traits Database
# Tundra Trait Team
# UV-B Radiation Sensitivity of Hieracium Pilosella
# Wetland Dunes Database

# datasets that are clearly inside range
keep.db <- c("Ecological Flora of the British Isles", 
             "PLANTATT - Attributes of British and Irish Plants",
             "Sheffield Database",
             "The Netherlands Plant Height Database")

# extract these records
try.keep <- try6.0 %>% filter(Dataset %in% keep.db)

# remove NAs from main database
try6.1 <- try6.0 %>% filter(!is.na(LAT_site))
  
# bind together
try.final <- rbind(try6.1, try.keep)

# extract list of unique ObservationID
obs.id.filter <- try.final %>% distinct(ObservationID)

obs.id.vector <- obs.id.filter$ObservationID
  
## info on control vs treatments comes from the covariates database
try6covs <- read_parquet("M:/# Manuscripts in progress/# Borealization/trait data/TRY_covariates_info.parquet")

## covariate table -- filter for everything else (not control) 
## remove records from observation table (based on ObsID) 
## which leaves control ones and those with no treatment 
## thanks Camila for this code!

try.final #contains the right species from the right area and traits

# only treatment 
treatment <- try6covs %>% filter(str_detect(DataName, "Treatment"))

keywords <- c("natural", "ambient", "none", "control", "not applicable")

# filtering  OrigValueSt containing keywords
treatment |> distinct(DataName, OrigValueStr, StdValueStr) |>
mutate(OrigValueStr = iconv(OrigValueStr, to = "UTF-8")) |> # Convert to proper encoding
  filter(str_detect(tolower(OrigValueStr), paste(keywords, collapse = "|"))) |> 
  write_csv("treatment.csv")

# standardize StdValueStr ambient manually updated
treatment_std <- read_csv("treatment.csv")

# Filtering only StdValueStr to have the data with information
treatment_words <- c("Ambient", "Control plot", "Natural environment", "Not applicable")

treatment_std <- treatment_std |> filter(!is.na(StdValueStr)) |> 
  filter(StdValueStr %in% treatment_words) |> 
  select(-StdValueStr)

treatment <- treatment |> anti_join(treatment_std) |> 
  distinct(ObservationID)

try.final.control <- try.final |> 
  anti_join(treatment, by = "ObservationID")
#62165 obs

save(try.final.control, file = "data/try_final_control.RData")

# checking species names that we have trait data for
sp.names.test <- try.final.control %>% distinct(AccSpeciesName)


## SUMMARY PER SPECIES ----
load("data/try_final_control.RData")

# Check if all units are the same - correct
units <- try.final.control %>% group_by(TraitNameNew) %>% 
  summarise(unique(UnitName)) %>% ungroup()

# How many records per species and trait?
records <- try.final.control %>% 
  group_by(AccSpeciesName, TraitNameNew) %>%
  summarise(n = n()) %>% ungroup()

# trait records vary widely from 1-4237
# establish a minimum of 5 records per species and trait

# Check outlier values
try.out <- try.final.control %>% 
  group_by(TraitNameNew, AccSpeciesName) %>% 
  mutate(NumRecords = length(AccSpeciesName)) %>% 
  ungroup() %>%
  filter(NumRecords > 4) %>%
  group_by(TraitNameNew, AccSpeciesName) %>%
  mutate(MeanTraitValue = mean(StdValue), 
         SD = sd(StdValue)) %>% 
    ungroup() %>%
    mutate(SD5 = SD*5,
         Mean_SD5 = MeanTraitValue + SD5,
         Outlier = case_when(StdValue > Mean_SD5 ~ "Yes", TRUE ~ "No"))

# examine potential outliers: 118 obs
pos.out <- try.out %>% filter(Outlier == "Yes")


# Final trait dataset (without outliers)
trait.db <- try.final.control %>% 
  group_by(AccSpeciesName, TraitNameNew) %>%
  mutate(NumRecords = length(AccSpeciesName)) %>% 
  ungroup() %>%
  filter(NumRecords > 4) %>%
  group_by(TraitNameNew, AccSpeciesName) %>%
  mutate(MeanTraitValue = mean(StdValue),
         SD = sd(StdValue)) %>% 
  ungroup() %>%
  mutate(SD5 = SD*5,
         Mean_SD5 = MeanTraitValue + SD5,
         Outlier = case_when(StdValue > Mean_SD5 ~ "Yes", TRUE ~ "No")) %>% 
  ungroup() %>%
  filter(Outlier == "No") %>%
  distinct(AccSpeciesName, TraitNameNew, .keep_all = TRUE) %>%
  select(AccSpeciesName, TraitNameNew, MeanTraitValue) %>%
  pivot_wider(values_from = MeanTraitValue, names_from = TraitNameNew) %>%
  rename(SPECIES_CLEAN = AccSpeciesName)

# NAs are introduced for those with no trait values



## SPECIES DATASET ----

# Join with species-level dataset
# Use species_ab_class.csv as default (as it doesn't contain the morphospecies)
species_ab_class <- read.csv("data/species_ab_class_jan2025.csv") %>% select(SPECIES_CLEAN, ClassNew)

# join this with the synonym list
species_ab_class_syno <- species_ab_class %>% 
  left_join(., sp.list.syno, by = "SPECIES_CLEAN") %>% select(-X)

# multiple matches but all the same synonym, leave just one row per species - back to 287 species
species_ab_class_syno_dist <- species_ab_class_syno %>%
  distinct(SPECIES_CLEAN, ClassNew, Synonym)


# join this with the categorical trait list
cat.trait.check <- read.csv("data/cat_trait_check.csv")

species_ab_class_syno_cat <- species_ab_class_syno_dist %>% 
  left_join(., cat.trait.check, by = "SPECIES_CLEAN") %>% select(-X)

# join this with the numerical trait dataset: by SPECIES_CLEAN
species_dataset <- species_ab_class_syno_cat %>% 
  left_join(., trait.db, by = "SPECIES_CLEAN")

# join with numerical trait dataset: by synonym
trait.db.syn <- trait.db %>% rename(Synonym = SPECIES_CLEAN)

species.dataset.syno <- species_ab_class_syno_cat %>% 
  left_join(., trait.db.syn, by = "Synonym")

# now we need to remove the species with data that came from the synonyms from the 1st dataset
# get the SPECIES_CLEAN name for the species with data from species.dataset.syno

# list of species with trait that came from the synonym database
species.dataset.syno.list <- c("Lysimachia europaea", "Cardamine digitalis", "Stellaria crassipes",
                              "Lycopodium complanatum", "Persicaria vivipara", "Oreojuncus trifidus",
                              "Oreomecon radicata", "Lagotis minor", "Carex × turfosa", "Carex microcarpa", 
                              "Carex fimbriata", "Festuca richardsonii", "Rumex alpestris",
                              "Salix daphnoides", "Potentilla fruticosa", "Cerastium glabratum")

species.dataset.syno.keep <- species.dataset.syno %>% 
  filter(SPECIES_CLEAN %in% species.dataset.syno.list)

# check that these species with synonym data are not being removed from the main db with other traits
species_dataset_check <- species_dataset %>%
  filter(SPECIES_CLEAN %in% species.dataset.syno.list)
# perfect - no trait data being removed here

# remove these species then from the original dataset
species_dataset_shorter <- species_dataset %>%
  filter(SPECIES_CLEAN %notin% species.dataset.syno.list)

# join with the synonym dataset: 287 species (without morphos)
species_traits_full <- rbind(species_dataset_shorter, species.dataset.syno.keep)

# For how many species do we have trait data?
# Total unique: 348, total without morphos: 287
n.height <- species_traits_full %>% filter(!is.na(PlantHeight)) # 191
n.sla <- species_traits_full %>% filter(!is.na(SLA))  # 166
n.seed <- species_traits_full %>% filter(!is.na(SeedMass)) # 83
n.leafn <- species_traits_full %>% filter(!is.na(LeafN)) # 120
n.leafcn <- species_traits_full %>% filter(!is.na(LeafCN)) # 54

write.csv(species_traits_full, "data/species_traits_dataset_full.csv")



