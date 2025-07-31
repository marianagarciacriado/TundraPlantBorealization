## Arctic plant borealization
## Mariana Garcia Criado 
## Null model
## May 2025



## PACKAGES ----
library(tidyverse)
library(brms)
library(ggeffects)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


# Set seed for reproducibility
set.seed(123)


## DATA ----

# join relevant datasets into one
species.class0 <- read.csv("data/species_ab_class_jan2025.csv") 

# ITEX data
load("data/boreal_master_perm.RData")

# join class to main database and removes morphospecies
master.class.com <- boreal_master_perm2 %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  tidyr::drop_na(Class) # all NA re morphospecies

# save smaller reference dataset
plots.only <- master.class.com %>% distinct(SiteSubsitePlot)

subsite.only <- master.class.com %>% 
  distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SiteSubsite)



## OBSERVED BIC ----

# Define trends per species and plot
master.class2.com <- master.class.com %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR)) %>% mutate(MaxYear = max(YEAR)) %>%
  filter(MinYear != MaxYear) %>%
  filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>%
  ungroup() %>%
  group_by(SiteSubsitePlot, SPECIES_CLEAN) %>% 
  mutate(SpeciesTimePoint = n()) %>% 
  mutate(trend = case_when(SpeciesTimePoint > 1 ~ "Persisting",
                           SpeciesTimePoint = 1 & YEAR == MinYear ~ "Extinct",
                           SpeciesTimePoint = 1 & YEAR == MaxYear ~ "Coloniser",
                           TRUE ~ "none")) %>% ungroup()

# calculate the Borealization Colonization Index 
master.class3.com <- master.class2.com %>% 
  filter(trend == "Coloniser") %>%
  group_by(SiteSubsitePlot, Class) %>% 
  summarise(n = n()) %>%
  #left_join(., duration.only, by = "SiteSubsitePlot") %>%
  ungroup() %>%
  group_by(SiteSubsitePlot) %>%
  mutate(TotalColonisers = sum(n)) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  mutate(BIC_B_LAB = (((Boreal + Low_Arctic_Boreal)/TotalColonisers)))

length(unique(master.class3.com$SiteSubsitePlot)) #891

# fill up with 0 the plots without colonisers (boreal or otherwise)
master.class4.com <- plots.only %>% 
  left_join(master.class3.com, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0) %>%
  left_join(., subsite.only, by = "SiteSubsitePlot")

# keep only observed values for comparison later on
bic.observed <- master.class4.com %>% select(SiteSubsitePlot, BIC_B_LAB)



## SIMULATED BCI ----


# for each site, extract the total site-level species list
master.list.site <- master.class2.com %>% 
  select(SITE, SPECIES_CLEAN, Zone_A:ClassNew) %>%
  distinct()

# this is the input file for randomization for now

# define the function
fun_null_sb_xyz <- function(data, plots, sp.list.site){
  
  sp_list_xyz <- sp.list.site %>%
    group_by(SITE) %>%
    # randomize species classes (here easiest to do by randomizing the species names)
    mutate(SPECIES_CLEAN_random = sample(SPECIES_CLEAN, n(), replace = FALSE)) %>%
    ungroup() %>%
    select(-SPECIES_CLEAN)
  
  result_sb_xyz <- data %>% 
    # remove the "real" classes
    select(-c(Zone_A,Zone_B,Zone_C,Zone_D,Zone_E,BorealSpecies,ArcticSpecies,ArcticEndemic,BorderlineSpecies,Acronym,Class,ClassNew)) %>%
    left_join(., sp_list_xyz, by=c("SPECIES_CLEAN" = "SPECIES_CLEAN_random", "SITE"="SITE")) %>%
    # calculate BIC
    filter(trend == "Coloniser") %>%
    group_by(SiteSubsitePlot, ClassNew) %>% 
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(SiteSubsitePlot) %>%
    mutate(TotalColonisers = sum(n)) %>%
    pivot_wider(names_from = ClassNew, values_from = n, values_fill = 0) %>%
    ungroup() %>%
    mutate(Null_BIC = (((`Boreal-tundra boundary` + `Boreal specialist`)/TotalColonisers))) %>%
    # fill up with 0s the plots without colonisers (boreal or otherwise)
    left_join(plots, ., by = "SiteSubsitePlot") %>%
    replace(is.na(.), 0) %>%
    dplyr::select(SiteSubsitePlot, Null_BIC) 
  
  return(result_sb_xyz)
}


# run 999 simulations
bic.null.sb.xyz <- do.call(rbind, replicate(999, fun_null_sb_xyz(data = master.class2.com, 
                                                                 plots = plots.only,
                                                                 sp.list.site = master.list.site), simplify = FALSE))


# calculate mean and 2SD range per plot
bic.null.sb.xyz2 <- bic.null.sb.xyz %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Null_Mean_BIC = mean(Null_BIC),
         Null_Sample_Size = length(Null_BIC),
         Null_2SD = sd(Null_BIC)*1.96) %>%
  select(-Null_BIC) %>%
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  mutate(Null_Lower_2SD = Null_Mean_BIC - Null_2SD,
         Null_Upper_2SD = Null_Mean_BIC + Null_2SD)


# merge the two: observed + null
bic.both.xyz <- left_join(bic.null.sb.xyz2, bic.observed, by = "SiteSubsitePlot")

# compare both values
bic.compare.xyz <- bic.both.xyz %>%
  mutate(Comparison_2SD = case_when(BIC_B_LAB > Null_Upper_2SD ~ "Greater",
                                   BIC_B_LAB < Null_Lower_2SD ~ "Lower",
                                   TRUE ~ "Similar"))

# how many times is it actually similar? 
bic.compare.summary.xyz <- bic.compare.xyz %>% group_by(Comparison_2SD) %>% summarise(n = n())

# with 2SD, 28 plots out of 1137 are different (2.5%)




## OBSERVED BAI ----

# already calculated original BAI values
orig.bia <- read.csv("data/bia-abund-lab-b-startend.csv")


## SIMULATED BAI ----

# define the function
fun_null_bia_xyz <- function(data, sp.list.site){
  
  sp_list_bia_xyz <- sp.list.site %>%
    group_by(SITE) %>%
    # randomize species classes (here easiest to do by randomizing the species names)
    mutate(SPECIES_CLEAN_random = sample(SPECIES_CLEAN, n(), replace = FALSE)) %>%
    ungroup() %>%
    select(-SPECIES_CLEAN)
  
  result_bia_xyz <- data %>% 
    # remove the "real" classes
    select(-c(Zone_A,Zone_B,Zone_C,Zone_D,Zone_E,BorealSpecies,ArcticSpecies,ArcticEndemic,BorderlineSpecies,Acronym,Class,ClassNew)) %>%
    left_join(., sp_list_bia_xyz, by=c("SPECIES_CLEAN" = "SPECIES_CLEAN_random", "SITE"="SITE")) %>%
    # calculate BIA per plot
    group_by(SiteSubsitePlot) %>% 
    mutate(MinYear = min(YEAR),
           MaxYear = max(YEAR)) %>%
    filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>%
    ungroup() %>%
    group_by(SiteSubsitePlotYear, ClassNew) %>% 
    mutate(CoverSumPerClass = sum(RelCover)) %>%
    ungroup() %>%
    distinct(SiteSubsitePlotYear, ClassNew, .keep_all = T) %>%
    select(-RelCover) %>%
    mutate(Timepoint = case_when(YEAR == MinYear ~ "Start",
                                 YEAR == MaxYear ~ "End")) %>%
    pivot_wider(names_from = c(Timepoint, ClassNew), values_from = CoverSumPerClass) %>%
    group_by(SiteSubsitePlot, .keep_all = T) %>% 
    summarise_each(funs(first(.[!is.na(.)]))) %>%
    ungroup() %>%
    # keep only relevant columns
    select(c(1, 37, 54:61)) %>%
    # replacing NAs by 0
    replace(is.na(.), 0) %>%
    mutate(Boreal_Start_Cover = `Start_Boreal specialist` + `Start_Boreal-tundra boundary`,
           Boreal_End_Cover = `End_Boreal specialist` + `End_Boreal-tundra boundary`,
           Boreal_Cover_Change = (Boreal_End_Cover - Boreal_Start_Cover)/Duration)
  
  return(result_bia_xyz)
  
}


# replicate 999 times
bia.null.xyz <- do.call(rbind, replicate(999, fun_null_bia_xyz(data = master.class.com, 
                                                                 sp.list.site = master.list.site), 
                                         simplify = FALSE))


# calculate mean and 2SD for each plot
bia.null.xyz2 <- bia.null.xyz %>% 
  group_by(SiteSubsitePlot) %>% 
  mutate(Null_Mean_BIA = mean(Boreal_Cover_Change),
         Null_Sample_Size = length(Boreal_Cover_Change),
         Null_2SD = sd(Boreal_Cover_Change)*1.96) %>%
  ungroup() %>%
  distinct(SiteSubsitePlot, .keep_all = TRUE) %>%
  mutate(Null_Lower_2SD = Null_Mean_BIA - Null_2SD,
         Null_Upper_2SD = Null_Mean_BIA + Null_2SD)


# bia observed - clean
bia.observed.clean <- orig.bia %>% select(SiteSubsitePlot, CoverChange)

# merge the two: observed + null
bia.both.xyz <- left_join(bia.null.xyz2, bia.observed.clean, by = "SiteSubsitePlot")

# compare both values
bia.compare.xyz <- bia.both.xyz %>%
  mutate(Comparison_2SD = case_when(CoverChange > Null_Upper_2SD ~ "Greater",
                                    CoverChange < Null_Lower_2SD ~ "Lower",
                                    TRUE ~ "Similar"))

# how many times is it actually similar? 
bia.compare.summary.xyz <- bia.compare.xyz %>% group_by(Comparison_2SD) %>% summarise(n = n())

# 84 are different (greater/lower) if we use the 2SD - 7.4%
