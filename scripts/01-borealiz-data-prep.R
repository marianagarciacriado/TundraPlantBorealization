## Arctic plant borealization
## Mariana Garcia Criado
## Script 1. ITEX data cleaning
## April 2024

# This script is provided to show the process of generating the main input file of plant composition 
# (boreal_master_dec2024, loaded in R as itex_full16.RData), but the original raw files to produce this file 
# are not provided in the Zenodo or GitHub respositories. 
# The original raw files from the ITEX+ database will be openly available when the data paper 
# led by Dr Anne Bjorkman is published.


## PACKAGES ----
library(tidyverse)
library(brms)
library(WorldFlora)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## DATA ----
load("data/input_itex/pfxy_all.RData")
load("data/input_itex/pfplot_all.RData")

# updated version from Oct2024 including Petr's extra years of collected data
# at Billefjorden and Pyramiden
load("data/input_itex/perccov_all_oct2024.RData") # 40718

# remove these data from the main as Petr sent an updated version
perccov_all_rem <- perccov_all %>% filter(SITE %notin% c("BILLEFJORDEN", "PYRAMIDEN"))

str(perccov_all_rem)

# load Petr's new data
petr <- read.csv("data/input_itex/PYR_AWS_fixed_data.csv")

# standardise into main format
petr2 <- petr %>% rename(SUBSITE = SITE) %>% 
  mutate(SITE = case_when(SUBSITE == "AWS" ~ "BILLEFJORDEN",
                          SUBSITE == "PYR" ~ "PYRAMIDEN"),
         ValueType = "percent_cover",
         TISSUE = NA) %>%
  select(., -GFNARROWarft)

petr2$PLOT <- as.character(petr2$PLOT)


# final cover database with updated data from Petr
perccov_all_final <- bind_rows(perccov_all_rem, petr2)  



# test if QHI data is correct in this ITEX version
qhi.test <- pfxy_all %>% filter(SITE == "QHI")
unique(qhi.test$YEAR) # this is missing years so we use the corrected version

qhi <- read.csv("data/input_itex/qhi-1999-2022-clean-nov22.csv") #43064

unique(perccov_all$SITE)
unique(pfxy_all$SITE)
unique(perccov_all$SITE)

# Consistent QHI variable structure
qhi.good <- qhi %>% select(-1)
qhi.good$PLOT <- as.character(qhi.good$PLOT)
qhi.good$X <- as.character(qhi.good$X)
qhi.good$Y <- as.character(qhi.good$Y)
qhi.good$HIT <- as.character(qhi.good$HIT)
qhi.good$ABUNDANCE <- as.numeric(qhi.good$ABUNDANCE)

# Removing QHI from ITEX dataset and include the most recent cleaned one
pfxy_all.qhi <- pfxy_all %>% filter(SITE != "QHI") %>% bind_rows(., qhi.good) #1560837


# Raw figures of plots and species number
length(unique(pfxy_all.qhi$PLOT)) #1019
length(unique(pfplot_all$PLOT)) #925
length(unique(perccov_all_final$PLOT)) #6810
#total = 8754

length(unique(pfxy_all.qhi$SPECIES_NAME)) #997
length(unique(pfplot_all$SPECIES_NAME)) #608
length(unique(perccov_all_final$SPECIES_NAME)) #613

sort(unique(pfxy_all.qhi$SITE)) # "ZACKENBERG", "LATNJA"
sort(unique(pfplot_all$SITE)) # KANGER, "JAMESONLAND", "JOATKA", "NAKKALA", "KILPISJARVI", "ABISKO", "LATNJA", "ADVENT"
sort(unique(perccov_all_final$SITE)) #"RIRI", "LOGH", "LORI", "FURI"


# Compare species list from Petr's site with input species name into WorldFlora
sort(unique(petr2$SPECIES_NAME))
# Equisetums are removed automatically because they are not shrubs/grams/forbs
# XXXBARE which doesn't count as species, "XXXDRABA:PYR_AWS", "XXXMOSS+LICHEN", "XXXPOA:PYR_AWS"
# "Minuartia rossii" is the only new species, but its status is accepted on WFO
# check that no synonyms are in the final databse (Arenaria rossii, Alsinanthe rossii, Alsinopsis rossii)


## STRUCTURE FIXES ----

# First of all fix Latnjajaure's name to Latnja to standardise
pfplot_all0 <- pfplot_all %>% mutate(SITE = case_when(SITE == "LATNJAJAURE" ~ "LATNJA", TRUE ~ SITE))


# Now fix plot names that include year - important because later on we will group by plot, not by plotXyear
pfplot_allX <- pfplot_all0 %>% 
  mutate(YearInPlot = ifelse(str_detect(PLOT, fixed(as.character(YEAR))), "YES", "NO")) %>%
  mutate(PLOT = ifelse(YearInPlot == "YES", str_remove(PLOT, pattern = as.character(YEAR)), PLOT)) %>%
  mutate(PLOT = str_remove(PLOT, pattern = "^_+"), # Removes leading "_"
         PLOT = str_remove(PLOT, pattern = "_$")) %>% # Removes trailing "_"  
  select(., -YearInPlot) 

# This is not the case in xy dataframe

perccov_allX <- perccov_all_final %>% 
  mutate(YearInPlot = ifelse(str_detect(PLOT, fixed(as.character(YEAR))), "YES", "NO")) %>%
  mutate(PLOT = ifelse(YearInPlot == "YES", str_remove(PLOT, pattern = as.character(YEAR)), PLOT)) %>%
  mutate(PLOT = str_remove(PLOT, pattern = "^_+"), # Removes leading "_"
         PLOT = str_remove(PLOT, pattern = "_$")) %>% # Removes trailing "_"  
  select(., -YearInPlot) 

unique(perccov_allX$SITE)

# Not permanently marked plots are identified in the manner "651a, 651b and 651c" (as plot names) 
# where the letter indicates a different sampling year (also identified in the year column).
# I am removing these as they will certainly include artificial extinctions/colonisations.

# These plots are all permanently marked
pfplot_allXX <- pfplot_allX %>% 
  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO")) 

unique(pfplot_allXX$PLOT)

# These are named with some letters but it's always the same plot name, and the year is specified correctly
pfxy_all.qhiXX <- pfxy_all.qhi %>% 
  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO"))

unique(pfxy_all.qhiXX$PLOT)

# These are all plots where a = start year and b, c = subsequent years
perccov_allXX <- perccov_allX %>% 
  mutate(LetterInPlot = ifelse(str_detect(PLOT, "[abc]$"), "YES", "NO"))
  
# confirm that this is the case for all
letter.check <- perccov_allXX %>% filter(LetterInPlot == "YES")

sort(unique(letter.check$SITE))

# Database without the letter plots
perccov_noletter <- perccov_allXX %>% filter(LetterInPlot == "NO")





## FILTERING ----

# Keep control plots only
unique(perccov_noletter$TREATMENT) # "CTL"     "CTLNG"   "OTC"     "CONTROL" "CTL_EXCLOSED" "WARMING"
unique(pfplot_allX$TREATMENT) # "OTC"     "CTL"     "OTCW"    "OTCWS"   "CONTROL" "WARMING" "DAMAGE""CTL_snowfence" "OTC_snowfence"
unique(pfxy_all.qhi$TREATMENT) # "OTC"      "CTL"      "CTLNG"    "CONTROL"  "WARMING"  NA         "CLT_GRA"  "CTL_NGRA"

# "CTLNG", "CTL_NGRA" are exclusion experiments so ambient conditions are with grazing included
#ctlng <- perccov_noletter %>% filter(TREATMENT == "CTLNG") #FAROE
#ctlxy <- pfxy_all.qhi %>% filter(TREATMENT %in% c("CTLNG",  "CTL_NGRA")) # "TIBET", "AUDKULUHEIDI"
#unique(ctlxy$SITE)

# control treatments only
trtmt.vector <- c("CTL", "CONTROL", "CLT_GRA", "CTL_snowfence", NA)

perccov_all2 <- perccov_noletter %>% filter(TREATMENT %in% trtmt.vector)
pfxy_all2 <- pfxy_all.qhi %>% filter(TREATMENT %in% trtmt.vector)
pfplot_all2 <- pfplot_allX  %>% filter(TREATMENT %in% trtmt.vector)

sort(unique(perccov_all2$SITE))



# only keeping alive hits and NA
unique(perccov_all2$STATUS) # "LIVE"         "LITTER"       "OTHER"        "DEAD"         "UNKNOWN"  NA
unique(pfplot_all2$STATUS) # "LIVE"   "LITTER"  "OTHER"  "DEAD"    "UNKNOWN"        
unique(pfxy_all2$STATUS) # "LIVE" "OTHER" "LITTER" "DEAD" "UNKNOWN" "STANDINGDEAD" NA "Standing dead" 
# "Live" "Litter" "" "Status"   "N/A" "Dead"

# Check the empty space - it's 2 records from QHI so we assume they are the same as NA
# emptyspace <- pfxy_all2 %>% filter(STATUS == "")


# Keep unknown as the majority of plots that have this category refer to all hits
# abiotic hits/litter will get filtered out below
status.vector <- c("LIVE", NA, "Alive", "Live", "N/A", "", "UNKNOWN")

perccov_all3 <- perccov_all2 %>% filter(STATUS %in% status.vector)
pfxy_all3 <- pfxy_all2 %>% filter(STATUS %in% status.vector)
pfplot_all3 <- pfplot_all2 %>% filter(STATUS %in% status.vector)

# quick site check
sort(unique(perccov_all3$SITE))
sort(unique(pfxy_all3$SITE))
sort(unique(pfplot_all3$SITE))

# Functional groups
unique(perccov_all3$GFNARROWwalker)
unique(pfplot_all3$GFNARROWwalker)
unique(pfxy_all3$GFNARROWwalker)

# Check FG = unknown/NA
unk.cov <- filter(perccov_all3, GFNARROWwalker %in% c("OTHER", "UNKNOWN", NA)) # "Bupleurum americanum"
unique(unk.cov$SPECIES_NAME)

unk.pfxy <- filter(pfxy_all3, GFNARROWwalker %in% c("XXXNA", "OTHER", NA))
unique(unk.pfxy$SPECIES_NAME)
  
unk.pfall <- filter(pfplot_all3, GFNARROWwalker %in% c("UNK", "OTHER")) 
unique(unk.pfall$SPECIES_NAME)



# Add in functional group to the clear vasculars
perccov_all4 <- perccov_all3 %>% 
  mutate(GFNARROWwalker = case_when(is.na(GFNARROWwalker) & SPECIES_NAME == "Bupleurum americanum" ~ "FORB", 
                                    TRUE ~ GFNARROWwalker))

# Retain vascular plants only
fg.vector <- c("FORB", "SEVER", "SDECI", "SHRUBU", "SHRUB", "GRAMINOIDU", "GRASS", "SEDGE", "RUSH", 
               "GRAMU", "WOODYU")

# Note: did NOT retain NA (as they're fixed below), and "XXXNA" as status = other and will be removed regardless because of dead/experimental.
perccov_all5 <- perccov_all4 %>% filter(GFNARROWwalker %in% fg.vector)
pfplot_all4 <- pfplot_all3 %>% filter(GFNARROWwalker %in% fg.vector)
pfxy_all4 <- pfxy_all3 %>% filter(GFNARROWwalker %in% fg.vector)

# Check abundance = NA
abun.na <- perccov_all5 %>% filter(is.na(ABUNDANCE))
abun.na.pf <- pfplot_all4 %>% filter(is.na(ABUNDANCE)) # all ok
abun.na.xy <- pfxy_all4 %>% filter(is.na(ABUNDANCE))
#1 NA abundance value for QHI:HE only, should be 1 instead as no 0s are recorded at this site

# Create unique plotXyear ID
perccov_all6 <- perccov_all5 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  #tidyr::replace_na(list(ABUNDANCE = 1)) %>% # "INCLINE_ULV", "INCLINE_LAV" which are NA should be 1 (see below)
  filter(SiteSubsite != "SADVENT:WET_PHOTO") # incomplete subsite, removing it to avoid issues in the cover calc below

pfplot_all5 <- pfplot_all4 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>% 
  select(., -COVER_UNDERSTORY)

pfxy_all5 <- pfxy_all4 %>% 
  unite(SiteSubsitePlotYear, c("SITE", "SUBSITE", "PLOT", "YEAR"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsitePlot, c("SITE", "SUBSITE", "PLOT"), sep = ":", remove = FALSE) %>%
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>%
  tidyr::replace_na(list(ABUNDANCE = 1)) #1 NA abundance value for QHI:HE only, should be 1 instead as no 0s are recorded


# "sum" doesn't belong there
unique(pfxy_all5$HIT)

# "SUM is Gavia, which will get filtered out below
# "sum" is Zackenberg so moving it to pfplot_all
summed <- pfxy_all5 %>% filter(HIT == "sum") %>% select(., -c(X, Y, HIT)) #66 obs

# Remove from XY database
pfxy_all6 <- pfxy_all5 %>% filter(HIT != "sum" | is.na(HIT))

# add in to sum database
pfplot_all6 <- rbind(pfplot_all5, summed)

sort(unique(pfplot_all6$SITE))




## COVER CONVERSION ----

# Do the cover values add up to 100?
cov <- perccov_all6 %>% filter(ValueType == "percent_cover") %>% 
  group_by(SiteSubsitePlotYear) %>% summarise(sum = sum(ABUNDANCE))
# Quite a lot of values over 100 so they need to be made proportional too so all values are comparable


#### Cover-equivalent ####

# Confirm that 1 row = 1 species
cov.test <- perccov_all6 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()

# Check if this is because of the missing species names or actually there are repeated species names
dif <- cov.test %>% filter(SameOrNot == "Different")

# Add up values per species so we end up with only one row per species
dif2 <- dif %>% 
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% 
  mutate(AbundanceFixed = sum(ABUNDANCE)) %>% 
  ungroup() %>% group_by(SiteSubsitePlotYear) %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()

# Confirm that this has worked and 1 row = 1 species
dif3 <- dif2 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()


# 1) Remove inconsistent plots from original dataset
perccov_all7 <- cov.test %>% filter(SameOrNot != "Different")

# 2) Dataframe with summed up values
sum.fixed <- dif2 %>% mutate(ABUNDANCE = AbundanceFixed) %>% select(., -AbundanceFixed)

# Bind all three into one fixed cover dataset
perccov_fixed0 <- rbind(perccov_all7, sum.fixed)

# Keep relevant columns only
perccov_fixed <- perccov_fixed0 %>% filter(ABUNDANCE > 0) %>% 
  select(., -c(NumberRows, NumberSpecies, SameOrNot))



# Convert all values to relative cover
itex.cov <- perccov_fixed %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% ungroup() 
# 13516 with Petr's sites

# Confirm that total cover values add up to 100 in every plotXyear
cov.check <- itex.cov %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)




#### Point-framing (summed) ####

# Confirm that 1 row = 1 species
pfsum.test <- pfplot_all6 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()

# Check if this is because of the missing species names or actually there are repeated species names
dif.sum <- pfsum.test %>% filter(SameOrNot == "Different")

# Duplicates: exactly the same species and abundance values on repeat
zac.vector <- c("ZACKENBERG:SALIX_SITE:S3_CONT:2009", "ZACKENBERG:SALIX_SITE:S1_CONT:2009",
               "ZACKENBERG:SALIX_SITE:S2_CONT:2009", "ZACKENBERG:SALIX_SITE:S4_CONT:2009",
               "ZACKENBERG:SALIX_SITE:S5_CONT:2009", "ZACKENBERG:CASSIOPE_SITE:C1_CONT:2009",
               "ZACKENBERG:CASSIOPE_SITE:C2_CONT:2009", "ZACKENBERG:CASSIOPE_SITE:C3_CONT:2009",
               "ZACKENBERG:CASSIOPE_SITE:C4_CONT:2009", "ZACKENBERG:CASSIOPE_SITE:C5_CONT:2009")

# Multiple plots have duplicate values: all records have exactly the same values twice
zac.dup <- dif.sum %>% filter(SiteSubsitePlotYear %in% zac.vector) %>% 
  group_by(SiteSubsitePlotYear) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% ungroup()

# Add up values per species so we end up with only one row per species
dif.sum2 <- dif.sum %>% filter(SiteSubsitePlotYear %notin% zac.vector) %>% 
  group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% 
  mutate(AbundanceFixed = sum(ABUNDANCE)) %>% ungroup() %>%
  group_by(SiteSubsitePlotYear) %>% distinct(SPECIES_NAME, .keep_all = TRUE) %>% 
  ungroup()

# Confirm that this has worked and 1 row = 1 species
dif.sum.test <- dif.sum2 %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(NumberRows = n()) %>% mutate(NumberSpecies = length(unique(SPECIES_NAME))) %>%
  mutate(SameOrNot = ifelse(NumberRows == NumberSpecies, "Same", "Different")) %>% ungroup()


# Dataframes to merge:

# 1) Remove inconsistent plots from original dataset
pfplot_all7 <- pfsum.test %>% filter(SameOrNot != "Different")

# 2) Dataframe with added values
dif.sum3 <- dif.sum2 %>% mutate(ABUNDANCE = AbundanceFixed) %>% select(., -AbundanceFixed)

# 3) Zackenberg with no duplicates (zac.dup)

# Bind all three into one fixed cover dataset
pfsum_fixed0 <- rbind(pfplot_all7, dif.sum3, zac.dup)

# Keep relevant columns only
pfsum_fixed <- pfsum_fixed0 %>% filter(ABUNDANCE > 0) %>% 
  select(., -c(NumberRows, NumberSpecies, SameOrNot))

# Fix Kilpisjarvi records (0.25 value for species in the plot but not touched by pin)
pfsum_fixed_kil <- pfsum_fixed %>% filter(SITE != "KILPSIJARVI" & ABUNDANCE != 0.25)


# Convert all values to relative cover
itex.pfsum <- pfsum_fixed_kil %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalAbundance = sum(ABUNDANCE)) %>%
  mutate(RelCover = (ABUNDANCE/TotalAbundance)*100) %>% ungroup() # 15532 obs

# Confirm that total cover values add up to 100 in every plotXyear
pfsum.check <- itex.pfsum %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)





#### Point-framing (XY) ####

# There are XY data that only have one coordinate - I'm filling in the NA cell so we avoid problems with cover calculation
na.x <- pfxy_all6 %>% filter(is.na(X)) # all filled in
na.y <- pfxy_all6 %>% filter(is.na(Y)) # 139062 obs

# Replace Y coords that are NA by 0s, create a unique coordinate
pfxy_all7 <- pfxy_all6 %>% tidyr::replace_na(list(Y = "0")) %>% 
  unite(XY, c("X", "Y"), sep = "_", remove = FALSE)

# No point on checking that 1 row = 1 species because multiple entries of the same species per plotXyear


# STEP 1: Convert species abundance to presence/absence 
# (2D, not considering multiple hits of the same species at each xy coord, just 1)
pfxy_all_pa <- pfxy_all7 %>% mutate(ABUNDANCE = ifelse(ABUNDANCE > 1, 1, ABUNDANCE)) %>%
  filter(ABUNDANCE > 0) %>%
  group_by(SiteSubsitePlotYear, XY) %>% 
  distinct(SPECIES_NAME, .keep_all = TRUE) %>% 
  ungroup()
  

# STEP 2: Calculate unique species hits per plot and total unique species hits per plot
pfxy_all_pa2 <- pfxy_all_pa %>% group_by(SiteSubsitePlotYear, SPECIES_NAME) %>% 
  mutate(UniqueSpHitsPlot = n()) %>% 
  distinct(SiteSubsitePlotYear, SPECIES_NAME, .keep_all = TRUE) %>% 
  ungroup() %>% select(., -c(X, Y, XY, HIT)) %>% group_by(SiteSubsitePlotYear) %>%
  mutate(TotalUniqueSpHitsPlot = sum(UniqueSpHitsPlot)) %>% ungroup()

# STEP 3: Calculate cover per species
pfxy_all_cov <- pfxy_all_pa2 %>% 
  mutate(RelCover = (UniqueSpHitsPlot/TotalUniqueSpHitsPlot)*100) #58639


# Confirm that total cover values add up to 100 in every plotXyear
pfxy.check <- pfxy_all_cov %>% group_by(SiteSubsitePlotYear) %>% 
  mutate(TotalCover = sum(RelCover)) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)



## BINDING ----

# Keep same number of relevant columns
itex.cov.f <- itex.cov %>% select(., -c(TotalAbundance, ABUNDANCE, LetterInPlot))
itex.pfsum.f <- itex.pfsum %>% select(., -c(TotalAbundance, ABUNDANCE)) 
pfxy_all_cov.f <- pfxy_all_cov %>% select(., -c(ABUNDANCE, UniqueSpHitsPlot, TotalUniqueSpHitsPlot))

# Bind all methods in one database
itex.all <- rbind(itex.cov.f, itex.pfsum.f, pfxy_all_cov.f) 
#87687 

# Quick site check
sort(unique(itex.all$SiteSubsite))
sort(unique(itex.all$SITE))



## SITE CHECKS ----

# Plots that had inconsistent surveyed areas/methods, identify species in the whole subsite
# didn't have permanent plots according to metadata, or had a different number of hits over the years
out.subsites <- c("SADVENT:WET_PHOTO", "SADVENT:MES_PHOTO",
                  "ALEXFIORD:LEVDOLOMITE", "ALEXFIORD:LEVGRANITE", 
                  "LATNJA:MESIC_MEADOW", "LATNJA:CAREX", "NAKKALA:8", "NAKKALA:9", "NAKKALA:10",
                  "CAPEBOUNTY:ARCTIC_WATERSHED_OBSERVATORY", "ADVENT:ADVENT_1",
                  "ADVENT:ADVENT_2", "ADVENT:ADVENT_3", "ADVENT:ADVENT_4", "ADVENT:ADVENT_5", 
                  "ADVENT:ADVENT_6", "ADVENT:ADVENT_7", 
                  "BYLOT:MESPOLYGON", "BYLOT:MESPRAIRIE", 
                  "SVERDRUP:SVERDRUP",
                  "STEPSTONES:TUNDRA2", "STEPSTONES:TUNDRA1",
                  "ZACKENBERG:SALIX_SITE", "ZACKENBERG:CASSIOPE_SITE")

# Because of diff number of hits: ADVENT, VARIABLE - 12 in 2016, 25 in 2017
# BYLOT: 100 hits per plot in first year; 81 hits per plot in 2008
# SVERDRUP: 0.25 m2 plots in 1992; 0.49 m2 in 2009
# Zackenberg salix_site: in 2009 only first hit recorded
# Zackenberg cassiope_site: in 2009 only first hit recorded
# SADVENT	WET_COVER and MES_PHOTO -- estimates from photos
# STEPSTONES -- estimates from photos
# Alex levgranite: 0.25 m2 plots in 1992, 0.49 m2 in 2009; 
# Alex levdolomite: 0.25 m2 plots in 1992, 0.49 m2 in 2009

# Incline: ULV (Ulvehaugen and all); 2018 data is pre-OTC set-up, i.e. all plots are not warmed then





# Remove inconsistent sites
itex.all2 <- itex.all %>% filter(SiteSubsite %notin% out.subsites) # 84859

# Remove inconsistent years from otherwise consistent plots
# ATQASUK only did top and bottom hits in first two survey periods (1995-97 and 2000) but done all hits since (as per email chain)
# BARROW: Note from Joe: only did top and bottom hits in first two survey periods (1995-97 and 2000) but done all hits since (as per email chain)
itex.all20 <- itex.all2 %>% filter(SITE %notin% c("ATQASUK", "BARROW"))

# keep only sites > year 2000 as consistent methods from there on
bobs.sites <- itex.all2 %>% filter(SITE %in% c("ATQASUK", "BARROW")) %>% filter(YEAR > 2000)

# bind together 
itex.all200 <- rbind(itex.all20, bobs.sites) #82095




## METADATA ----

# Updated metadata file November 2024
metadata0 <- read.csv("data/input_itex/itex_metadata_nov2024.csv")

# Keep only relevant columns
metadata <- metadata0 %>% 
  unite(SiteSubsite, c("SITE", "SUBSITE"), sep = ":", remove = FALSE) %>% 
  select(SiteSubsite, Moisture.category, LAT, LONG, ELEV, AZONE, SurveyedArea, 
         Data.Contributor, DomGrazer, GrazerIntensity)

# Merge with composition data
itex.full <- left_join(itex.all200, metadata, by = "SiteSubsite")

# Keep N hemisphere, boreal and Arctic plots only (all sizes)
itex.full2 <- itex.full %>% filter(LAT > 40| is.na(LAT))


# Check for missing NA
unique(itex.full2$Moisture.category)
unique(itex.full2$SurveyedArea)
unique(itex.full2$LAT)
unique(itex.full2$LONG)
unique(itex.full2$ELEV)
unique(itex.full2$DomGrazer)
unique(itex.full2$GrazerIntensity)

# Which subsites didn't get metadata?
nomoist <- itex.full2 %>% filter(Moisture.category == "" | is.na(Moisture.category))
sort(unique(nomoist$SiteSubsite)) #14 subsites

nolat <- itex.full2 %>% filter(is.na(LAT))
sort(unique(nolat$SiteSubsite)) #4 subsites without latitude data

nolong <- itex.full2 %>% filter(is.na(LONG))
sort(unique(nolat$SiteSubsite)) #4 subsites without longitude data

nosize <- itex.full2 %>% filter(SurveyedArea == "" |is.na(SurveyedArea))
sort(unique(nosize$SiteSubsite)) # 1 subsite

# Fill in missing surveyed area - thanks Joe Everest for Jameson code!
itex.full3 <- itex.full2 %>% 
  mutate(SurveyedArea = ifelse(SiteSubsite == "AUYUITTUQ:OWL RIVER", 1, SurveyedArea), 
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT"), 0.2178, SurveyedArea),
         SurveyedArea = ifelse(SiteSubsite %in% c("JAMESONLAND:TYSKIT") & 
                                 str_detect(PLOT, paste(c("13", "CASEMP", "BBN"), collapse = '|')), 0.1089, SurveyedArea))


# Identify sites known to have LAT and LONG that fall in the ocean (identified during climate extraction)
water.LAT.LONG <- itex.full3 %>% 
  dplyr::select(SiteSubsite, YEAR, LAT, LONG) %>% 
  filter(LAT %in% c(74.28, 74.29)) %>% 
  distinct(SiteSubsite, .keep_all = TRUE)

water.LAT.LONG <- unique(water.LAT.LONG$SiteSubsite)

# Incorrect LAT and LONG: 8 x ZACKENGERG subsites

# Manually input new LAT and LONG info for those missing or in water (all above 60 so not removed for being 'non-Arctic')
itex.full4 <- itex.full3 %>% 
  mutate(LAT = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", 69.27, LAT),
         LAT = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", 69.43, LAT),
         LAT = ifelse(SiteSubsite %in% water.LAT.LONG, 74.47427, LAT),
         LONG = ifelse(SiteSubsite == "DISKO:DRYHEATH_FLUX", -53.45, LONG),
         LONG = ifelse(SiteSubsite == "DISKO:WETFEN_FLUX", -53.78, LONG),
         LONG = ifelse(SiteSubsite %in% water.LAT.LONG,-20.52895, LONG),
         LAT = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", 61.16, LAT),
         LONG = ifelse(SiteSubsite == "NARSARSUAQ:HIGH_ELEVATION", -45.40, LONG),
         LAT = ifelse(SiteSubsite == "KYTALYK:LAKEBED", 70.81, LAT))
# Info taken from paper
# https://cdnsciencepub.com/doi/10.1139/AS-2020-0041 and
# https://iopscience.iop.org/article/10.1088/1748-9326/6/3/035502/meta 

# These have all been inputted now
nolat2 <- itex.full4 %>% filter(is.na(LAT))
nolong2 <- itex.full4 %>% filter(is.na(LONG))

# Spot checks
hist(itex.full4$LAT) # makes sense
hist(itex.full4$LONG) # ok
unique(itex.full4$SiteSubsitePlot)

unique(itex.full4$RelCover)
sort(unique(itex.full4$SITE))



## SITE SELECTION ----

# Extract subsite-level coordinates
subs.coord <- itex.full4 %>% distinct(SiteSubsite, .keep_all = TRUE) %>% 
  dplyr::select(SiteSubsite, LAT, LONG)

# save subsites and coordinates
write.csv(subs.coord, "data/subs_coords.csv")

# Site selection occurs on QGIS, retaining subsites in zones 6 and 11 of WWF ecoregions
# (Arctic and boreal)
# Load subsite_coords.csv on QGIS and make a list of sites that don't overlap zones

# remove sites outside zones 6 and 11
low.sites <- c("FAROE", "NIWOT", "TAISETSU", "GAVIA", "VALBERCLA")

itex.full5 <- itex.full4 %>% filter(SITE %notin% low.sites)

# all makes sense
sort(unique(itex.full5$SITE))


# How many plots were surveyed more than once and for a minimum of 5 years?
duration.filter <- itex.full5 %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% ungroup() %>%
  filter(Duration > 4) %>% distinct(SiteSubsitePlot, .keep_all = TRUE) 
# 1293

duration.vector <- duration.filter$SiteSubsitePlot

# filter only those plots with min 5 years monitoring
itex.full6 <- itex.full5 %>% filter(SiteSubsitePlot %in% duration.vector) 
# 41247

# down to 34 sites
sort(unique(itex.full6$SITE))




## SPECIES NAMES ----
spp <- unique(itex.full6$SPECIES_NAME)

# Identify empty species names
empty <- itex.full6 %>% filter(SPECIES_NAME == " ") # no empty cells

# Remove trailing white space so the same species are comparable
itex.full6 <- itex.full6 %>% mutate(SPECIES_NAME = str_trim(SPECIES_NAME))

# Check full species name again
unique(itex.full6$SPECIES_NAME)

# Standardise subspecies to species level
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Ledum palustre subsp. groenlandicum"] <- "Ledum palustre"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Ledum palustre subsp. decumbens"] <- "Ledum palustre"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Tephroseris integrifolia subsp. atropurpurea"] <- "Tephroseris integrifolia"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Silene uralensis subsp. apetala"] <- "Silene uralensis"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Cardamine bellidifolia subsp. alpina"] <- "Cardamine bellidifolia"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Carex aquatilis var. minor"] <- "Carex aquatilis"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Eriophorum angustifolium subsp. triste"] <- "Eriophorum angustifolium"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Empetrum nigrum subsp. hermaphroditum"] <- "Empetrum nigrum"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Anthoxanthum odoratum subsp. nipponicum"] <- "Anthoxanthum odoratum"
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Silene uralensis subsp. apetala"] <- "Silene uralensis"  
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Salix lanata subsp. richardsonii"] <- "Salix lanata"  
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Rumex alpestris subsp. lapponicus"] <- "Rumex alpestris"   
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Empetrum nigrum/Empetrum nigrum subsp. hermaphroditum"] <- "Empetrum nigrum"   
itex.full6$SPECIES_NAME[itex.full6$SPECIES_NAME == "Salix arctica/Salix arctica"] <- "Salix arctica"



# Check that species name are correct
library(WorldFlora)
WFO.remember('C:/TundraDivHub/WFO_Backbone/classification.csv')
WFO.remember('WFO_Backbone/classification.csv')

# dataframe with unique species names only
sp.names.only <- itex.full6 %>% distinct(SPECIES_NAME)
write.csv(sp.names.only, "data/species_names_list_dec2024.csv")

#sp.names.list <- read.csv("data/species_names_list_dec2024.csv")

# run the taxonomic checker
taxon.check <- WFO.match(spec.data = sp.names.list, 
                         spec.name = "SPECIES_NAME",
                         WFO.file = 'WFO_Backbone/classification.csv',
                         no.dates = TRUE)

write.csv(taxon.check, "data/worldflora_taxon_check_oct2024.csv")

#taxon.check <- read.csv("data/worldflora_taxon_check_oct2024.csv")

# check those that don't match
tofix <- taxon.check %>% 
  mutate(SPECIES_CLEAN = case_when(SPECIES_NAME != scientificName ~ scientificName,
                                                          TRUE ~ SPECIES_NAME)) %>%
  select(SPECIES_NAME, scientificName, family, SPECIES_CLEAN) %>%
  distinct(SPECIES_NAME, .keep_all = TRUE)

# convert SPECIES_CLEAN genus to morphospecies
genus.morp.vector <- scan(text = "Anemone
Antennaria Arnica Astragalus Calamagrostis Cardamine
Carex Deschampsia Draba Eriophorum Festuca Gentiana
Hedysarum Dryas Tofieldia Hepatica Juncus Luzula
Minuartia Oxytropis Pedicularis Petasites Poa
Polemonium Polygonum Sagina Salix Saxifraga
Senecio Stellaria", what="")

# Transform genus only into morphospecies
tofix2 <- tofix %>% 
  mutate(SPECIES_CLEAN = ifelse(SPECIES_CLEAN %in% genus.morp.vector, 
                                paste0("XXX", SPECIES_CLEAN), SPECIES_CLEAN)) %>%
  mutate(SPECIES_CLEAN = ifelse(SPECIES_NAME %in% c("Carex microchaeta&rupestris",
                                                    "Deschampsia flexuosa/Juncus trifidus",
                                                    "Empetrum nigrum/Phyllodoce caerulea",
                                                    "Eriophorum scheuchzeri/chamissonis",
                                                    "Luzula spicata/confusa",
                                                    "Minuartia NA",
                                                    "Salix arctica/arctophila",
                                                    "Graminoid unknown"),
                                paste0("XXX", SPECIES_NAME), SPECIES_CLEAN)) %>%
  mutate(SPECIES_CLEAN = ifelse(SPECIES_NAME == "Cardamine digitalis", 
                                "Cardamine digitalis", SPECIES_CLEAN)) %>%
  mutate(SPECIES_CLEAN = ifelse(SPECIES_NAME == "Deschampsia caespitosa", 
                                "Deschampsia cespitosa", SPECIES_CLEAN)) 


# remove ssp introduced by the taxo checker
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Carex marina  subsp. marina"] <- "Carex marina"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Deschampsia cespitosa  subsp. cespitosa"] <- "Deschampsia cespitosa"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Doronicum grandiflorum  subsp. grandiflorum"] <- "Doronicum grandiflorum"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Empetrum nigrum  subsp. hermaphroditum"] <- "Empetrum nigrum"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Eriophorum angustifolium  subsp. triste"] <- "Eriophorum angustifolium"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Hedysarum boreale  subsp. mackenziei"] <- "Hedysarum boreale"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Tephroseris integrifolia  subsp. atropurpurea"] <- "Tephroseris integrifolia"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Arabidopsis lyrata  subsp. petraea"] <- "Arabidopsis lyrata"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Arnica griscomii  subsp. frigida"] <- "Arnica griscomii"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Artemisia norvegica  subsp. saxatilis"] <- "Artemisia norvegica"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Saussurea angustifolia  var. viscida"] <- "Saussurea angustifolia"
tofix2$SPECIES_CLEAN[tofix2$SPECIES_CLEAN == "Rumex aquaticus  subsp. arcticus" ] <- "Rumex aquaticus" 



# Double-check species names
spppp <- sort(unique(tofix2$SPECIES_CLEAN))

# Check empty families
na.fam <- tofix2 %>% filter(is.na(family))

sort(unique(na.fam$SPECIES_CLEAN))

# fill in empty families that we have enough resolution to know
# extract families from tofix2
tofix2$family[tofix2$SPECIES_CLEAN == "XXXALCHEMILLA"] <- "Rosaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXasteraceae"] <- "Asteraceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXCAMPANULA"] <- "Campanulaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXcarex"] <- "Cyperaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXCAREX"] <- "Cyperaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXCERASTIUM"] <- "Caryophyllaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXDRABA"] <- "Brassicaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXDUPONTIA"] <- "Poaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXEPILOBIUM"] <- "Onagraceae"                          
tofix2$family[tofix2$SPECIES_CLEAN == "XXXJuncus"] <- "Juncaceae"  
tofix2$family[tofix2$SPECIES_CLEAN == "XXXkobresia"] <- "Cyperaceae"  
tofix2$family[tofix2$SPECIES_CLEAN == "XXXluzula"] <- "Juncaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXLUZULA"] <- "Juncaceae"             
tofix2$family[tofix2$SPECIES_CLEAN == "XXXMinuartia NA"] <- "Caryophyllaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXotheraster"] <- "Asteraceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXoxytropis"] <- "Fabaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXPEDICULARIS"] <- "Orobanchaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXPOA"] <- "Poaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXPOLYGONUM"] <- "Polygonaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXPyrola"] <- "Ericaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXPYROLA"] <- "Ericaceae"
tofix2$family[tofix2$SPECIES_CLEAN == "XXXTARAXACUM"] <- "Asteraceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXTofieldia"] <- "Tofieldiaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXunkOxytropis"] <- "Fabaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXVIOLA"] <- "Violaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXWOODYDRYAS"] <- "Rosaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXWOODYSALIX"] <- "Salicaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXTaraxacum"] <- "Asteraceae"  
tofix2$family[tofix2$SPECIES_CLEAN == "XXXSTELLARIA"] <- "Caryophyllaceae"  
tofix2$family[tofix2$SPECIES_CLEAN == "XXXSalix"] <- "Salicaceae" 
tofix2$family[tofix2$SPECIES_CLEAN == "XXXSALIX"] <- "Salicaceae"


# check again 
na.fam2 <- tofix2 %>% filter(is.na(family))

# fill in NA and empty families that are not possible to figure out
tofix3 <- tofix2 %>% mutate(family = case_when(is.na(family) ~ "Unknown",
                                               family == "" ~ "Unknown",
                            TRUE ~ family)) %>%
  select(., -scientificName)

write.csv(tofix3, "data/clean_sp_names_families_dec2024")




## FUNCTIONAL GROUPS ----

# join new species names with full database
itex.full7 <- itex.full6 %>% left_join(., tofix3, by = "SPECIES_NAME") %>%
  mutate(SPECIES_CLEAN = case_when(SPECIES_NAME == "Salix lanata x reticulata" ~ "Salix lanata x reticulata",
                                   SPECIES_NAME == "XXXDRABA:PYR_AWS" ~ "XXXDRABA:PYR_AWS",
                                   SPECIES_NAME == "XXXPOA:PYR_AWS" ~ "XXXPOA:PYR_AWS",
                                   TRUE ~ SPECIES_CLEAN)) %>%
  mutate(family = case_when(SPECIES_NAME == "Salix lanata x reticulata" ~ "Salicaceae",
                            SPECIES_NAME == "XXXDRABA:PYR_AWS" ~ "Brassicaceae",
                            SPECIES_NAME == "XXXPOA:PYR_AWS" ~ "Poaceae",
                                   TRUE ~ family))

# species names looking good
sort(unique(itex.full7$SPECIES_CLEAN))

# check functional groups
unique(itex.full7$GFNARROWwalker)


# Check that 1 func group per species and per genus

# How many func groups per genus?
gf.check <- itex.full7 %>% 
  group_by(GENUS) %>% 
  mutate(count_fg = length(unique(GFNARROWwalker))) %>% 
  ungroup() %>% 
  filter(count_fg > 1, !is.na(GENUS)) %>% 
  dplyr::select(SPECIES_CLEAN, GENUS, GFNARROWwalker, count_fg) %>% 
  unique() %>% 
  arrange(GENUS)


# Manually correct incorrect records (Potentillas are ok)
itex.full8 <- itex.full7 %>% 
  mutate(GFNARROWwalker = case_when(SPECIES_CLEAN == "XXXCarex microchaeta&rupestris" ~ "SEDGE",
                                    SPECIES_CLEAN == "Diapensia lapponica" ~ "SEVER",
                                    SPECIES_CLEAN == "XXXkobresia" ~ "SEDGE",
                                    SPECIES_CLEAN == "Rhododendron tomentosum" ~ "SEVER",
                                    SPECIES_CLEAN == "Linnaea borealis" ~ "SEVER",
                                    SPECIES_CLEAN == "Dryas octopetala" ~ "SEVER",
                                    TRUE ~ GFNARROWwalker))

# All good now - remaining genera with multiple FGs are correct (manually checked)
gf.check2 <- itex.full8 %>% 
  group_by(GENUS) %>% 
  mutate(count_fg = length(unique(GFNARROWwalker))) %>% 
  ungroup() %>% 
  filter(count_fg > 1, !is.na(GENUS)) %>% 
  dplyr::select(SPECIES_CLEAN, GENUS, GFNARROWwalker, count_fg) %>% 
  unique() %>% 
  arrange(GENUS)

# Check that all species have a func group
na.fg.check <- itex.full8 %>% filter(is.na(GFNARROWwalker)) # all good

sort(unique(itex.full8$GFNARROWwalker))

# Standardise and add woodiness
itex.full9 <- itex.full8 %>% 
  mutate(FunctionalGroup = case_when(GFNARROWwalker == "FORB" ~ "Forb",
                                     GFNARROWwalker %in% c("GRAMINOIDU", "GRAMU", "GRASS", "RUSH", "SEDGE") ~ "Graminoid",
                                     GFNARROWwalker %in% c("SDECI", "SEVER") ~ "Shrub")) %>%
  mutate(Woodiness = case_when(FunctionalGroup == "Shrub" ~ 1, 
                               TRUE ~ 0)) %>%
  mutate(Deciduousness = case_when(GFNARROWwalker == "SEVER" ~ "Evergreen",
                                   GFNARROWwalker == "SDECI" ~ "Deciduous",
                                   TRUE ~ "Not applicable"))


# Berry categories
berry <- read.csv("data/berry-species-list.csv")

itex.full10 <- itex.full9 %>% left_join(., berry, by = "SPECIES_CLEAN") %>%
  mutate(BerryCat = case_when(is.na(BerryCat) ~ "Not berry", TRUE ~ BerryCat))


# N fixers
nfix <- read.csv("data/n-fixers-species-list.csv")

itex.full11 <- itex.full10 %>% left_join(., nfix, by = "GENUS") %>%
  mutate(N_fixer = case_when(is.na(N_fixer) ~ 0, TRUE ~ N_fixer))


         
## GEO DATA ----

# Check sites
sort(unique(itex.full11$SITE)) # 34 sites

# group sites by biogeographical regions
eurasia <- c("ENDALEN", "SADVENT", "KYTALYK", "ABISKO", "FURI", "JOATKA", "KILPISJARVI", 
             "LATNJA", "LOGH", "LORI", "RIRI", "INCLINE_GUD", "INCLINE_SKJ", "INCLINE_ULV", 
             "INCLINE_LAV", "PYRAMIDEN", "BILLEFJORDEN")
greenice <- c("KANGER", "ZACKENBERG", "DISKO", "AUDKULUHEIDI", "THINGVELLIR", "THUFUVER")
na.east <- c("ALEXFIORD", "BYLOT", "AUYUITTUQ")
na.west <- c("ANWR", "ATQASUK", "BARROW", "BROOKS", "QHI", 
             "TOOLIK", "KLUANE", "WOLFCREEK", "TORNGATS")

# add in regions to the database
itex.full12 <- itex.full11 %>%
  mutate(lat_grid = plyr::round_any(LAT, 0.5, f = floor)) %>%
  mutate(lon_grid = ifelse(LONG >0, plyr::round_any(LONG, 0.5, f = floor),
                           plyr::round_any(LONG, 0.5, f = ceiling))) %>%
  mutate(gridcell = paste0("_", lat_grid, "_", lon_grid)) %>%
  mutate(Region = ifelse(SITE %in% eurasia, "Eurasia",
                         ifelse(SITE %in% greenice, "GreenIceLand",
                                ifelse(SITE %in% na.east, "North America-East",
                                       ifelse(SITE %in% na.west, "North America-West", NA)))))



## MORPHOSPECIES ----
length(unique(itex.full12$SiteSubsitePlotYear)) # 5586 plotXyear in total

# Check the unidentified species names
xxxtest <- itex.full12 %>% 
  filter(str_detect(SPECIES_CLEAN, 'XXX|xxx')) %>% 
  group_by(SiteSubsitePlotYear) %>% 
  mutate(MorphoCover = sum(RelCover)) %>%
  ungroup() %>% select(SITE, SiteSubsite, SiteSubsitePlot, SiteSubsitePlotYear, MorphoCover) %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE)

# 1962 unique plotXyear that contain morphospecies
length(unique(xxxtest$SiteSubsitePlotYear))

# Check how many plots would be removed with different cut-offs
morpho10 <- xxxtest %>% filter(MorphoCover > 10) # 10% morphospecies cutoff would remove 889 plotXyear (16% of database)
morpho15 <- xxxtest %>% filter(MorphoCover > 15) # 15% morphospecies cutoff would remove 666 plotXyear (12% of database)
morpho20 <- xxxtest %>% filter(MorphoCover > 20) # 20% morphospecies cutoff would remove 529 plotXyear (9.5% of database)
morpho25 <- xxxtest %>% filter(MorphoCover > 25) # 25% morphospecies cutoff would remove 406 plotXyear (7.3% of database)

# Mean morphospecies cover is 16.2%
mean(xxxtest$MorphoCover)

# Visualize morphospecies cover across plots
(xxxmorph <- ggplot(xxxtest, aes(x=MorphoCover)) + 
    geom_histogram(binwidth = 0.5) + 
    geom_vline(aes(xintercept = mean(MorphoCover)), colour = "red", linetype = "dashed", size = 1) +
    xlab("Cover of morphospecies in plots"))

# Extract as vector from the main dataset those plotXyear with >10% morphospecies cover
morpho10.vec <- unique(morpho10$SiteSubsitePlotYear)

# save this vector
write.csv(morpho10.vec, "data/morpho10_vector.csv")

# Remove those plots
itex.full13 <- itex.full12 %>% filter(SiteSubsitePlotYear %notin% morpho10.vec) #35398

sort(unique(itex.full13$SiteSubsite)) #116 subSites

# Apply 5-year duration filter again because after removing plotXyears due to morphospecies
# that leaves some subsites with one year per plot only
duration2 <- itex.full13 %>% 
  distinct(SiteSubsitePlotYear, .keep_all = TRUE) %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% 
  ungroup() %>%
  filter(Duration > 4) %>% 
  distinct(SiteSubsitePlot, .keep_all = TRUE) # 1131 plots

# remove the shorter plots
duration.vector2 <- duration2$SiteSubsitePlot

itex.full14 <- itex.full13 %>% filter(SiteSubsitePlot %in% duration.vector2) # 34550

sort(unique(itex.full14$SiteSubsite)) #113 subsites, removed 3 by re-applying the filter


# let's double check that all sites have a minimum of 5-year duration
itex.full15 <- itex.full14 %>% group_by(SiteSubsitePlot) %>% 
  mutate(Duration = max(YEAR) - min(YEAR)) %>% ungroup()
  


# coordinates for subsites to extract distance to treeline (this was done before final site selection/filter)
subs.treeline <- itex.full15 %>% 
  distinct(SiteSubsite, .keep_all = T) %>% 
  select(SiteSubsite, LAT, LONG)

write.csv(subs.treeline, "data/subsites_treeline_oct2024.csv")

sites <- unique(itex.full15$SITE) # 32 sites retained
subsite <- unique(itex.full15$SiteSubsite) # 113 subsites
plots <- unique(itex.full15$SiteSubsitePlot) # 1137 plots





## METADATA CLEANING ----
unique(itex.full15$Moisture.category) #all ok
unique(itex.full15$SurveyedArea) # all ok
unique(itex.full15$ELEV) # all ok
unique(itex.full15$DomGrazer) # much to standardise
unique(itex.full15$GrazerIntensity) # needs standardising

# standardise categorical variables
itex.full16 <- itex.full15 %>% 
  rename(Moisture = Moisture.category) %>%
  mutate(Moisture = case_when(Moisture == "MOIST" ~ "Moist",
                              Moisture == "DRY" ~ "Dry",
                              Moisture == "WET" ~ "Wet",
                              Moisture == "MIXED" ~ "Mixed")) %>%
  mutate(GrazerIntensity = case_when(GrazerIntensity %in% c("LOW", "low", "very low") ~ "Low",
                                     GrazerIntensity == "high" ~ "High",
                                     GrazerIntensity == "medium" ~ "Medium",
                                     GrazerIntensity == "" ~ NA,
                                     TRUE ~ GrazerIntensity)) %>%
  #mutate(ELEV = case_when(ELEV == "" ~ NA, TRUE ~ ELEV)) %>%
  mutate(DomGrazer = case_when(DomGrazer %in% c("insects small", "insects") ~ "Insects",
                               DomGrazer %in% c("large", "LARGE", "large (reindeer)") ~ "Large",
                               DomGrazer %in% c("small", "Small") ~ "Small",
                               DomGrazer %in% c("none", "none ", "no vertebrate grazers") ~ "None",
                               DomGrazer %in% c("large, small", "large,small", "large.small", "insects, large", "large, birds") ~ "Mixed",
                               DomGrazer == "birds" ~ "Birds",
                               DomGrazer == "" ~ NA,
                               TRUE ~ DomGrazer))

# all ok now - missing NA values are in Brooks (never heard back from PI)
unique(itex.full16$DomGrazer) 
unique(itex.full16$GrazerIntensity)
unique(itex.full16$Moisture)

# check variable type
str(itex.full16)

itex.full16$ELEV <- as.numeric(itex.full16$ELEV)
itex.full16$SurveyedArea <- as.numeric(itex.full16$SurveyedArea)

# save mastersheet
save(itex.full16, file = "data/boreal_master_dec2024.RData")

