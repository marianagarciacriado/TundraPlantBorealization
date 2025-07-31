## Arctic plant borealization
## Mariana Garcia Criado
## Script 6. Community-level analyses
## May 2024


## PACKAGES ----
library(tidyverse)
library(broom)
library(brms)
library(corrplot)
library(viridis)
library(ggOceanMaps)
library(ggrepel)
library(ggeffects)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## DATA ----

# join relevant datasets into one
species.class0 <- read.csv("data/species_ab_class_jan2025.csv") 

# ITEX data
load("data/boreal_master_perm.RData")

# join class to main database and removes morphospecies
master.class.com <- boreal_master_perm2 %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  tidyr::drop_na(Class) # all NA re morphospecies

## TOTAL RECORDS: 33388 (final figure of records from this file as it's the one without morphospecies)

# database summary stats
length(unique(master.class.com$SITE)) #32
length(unique(master.class.com$SiteSubsite)) #113
length(unique(master.class.com$SiteSubsitePlot)) #1137
min(master.class.com$YEAR) #1981
max(master.class.com$YEAR) #2023
length(unique(master.class.com$SPECIES_CLEAN)) #287

# plots per subsite
plotsXsubs <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SiteSubsite) %>% group_by(SiteSubsite) %>%
  summarise(n = n())

mean(plotsXsubs$n) #10 plots per subsite (1-83)


# subsites per study area
subsXsite <- master.class.com %>% distinct(SiteSubsite, .keep_all = T) %>% 
  select(SiteSubsite, SITE) %>% group_by(SITE) %>%
  summarise(n = n())

mean(subsXsite$n) #3.5 subsites per study area (1-31)


# duration per plot
durXplot <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, Duration)

mean(durXplot$Duration)
# 5-28, mean 15

# save smaller reference datasets
duration.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, Duration)
write.csv(duration.only, "data/duration_only.csv")

site.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SITE)
write.csv(site.only, "data/site_plot_only.csv")

subsite.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SiteSubsite)
write.csv(site.only, "data/subsite_plot_only.csv")

plots.only <- master.class.com %>% distinct(SiteSubsitePlot)
write.csv(plots.only, "data/plots_only.csv")

site.long <- master.class.com %>% distinct(SITE, .keep_all = T) %>%
  select(SITE, LONG)
write.csv(site.long, "data/site_long.csv")



## 1) BCI CALCULATION ----
master.class2.com <- master.class.com %>% group_by(SiteSubsitePlot) %>% 
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
  ungroup() %>%
  group_by(SiteSubsitePlot) %>%
  mutate(TotalColonisers = sum(n)) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  left_join(., duration.only, by = "SiteSubsitePlot") %>%
  mutate(BIC_B_LAB = (((Boreal + Low_Arctic_Boreal)/TotalColonisers)),
         BIC_B = (Boreal/TotalColonisers))

length(unique(master.class3.com$SiteSubsitePlot)) #891

# plots without colonisers (boreal or otherwise) are missing from this dataset
# eed to fill them up with 0 as they still count as zero borealization (even if they had no colonisations at all)
master.class4.com <- plots.only %>% 
  left_join(master.class3.com, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0) %>%
  left_join(., subsite.only, by = "SiteSubsitePlot") %>%
  left_join(., site.only, by = "SiteSubsitePlot") %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))

length(unique(master.class4.com$SiteSubsitePlot)) #1137

# Are shorter/longer timeseries picking up greater rates of borealization?
(dur.bic <- ggplot(master.class4.com, aes(x = Duration, y = BIC_B_LAB)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))

# model
dur.bic.mod.re <- brm(bf(BIC_B_LAB ~ Duration + (1|SiteSubsite), zoi ~ 1, coi ~ 0),
                      data = master.class4.com, 
                      family = 'zero_one_inflated_beta',
                      prior = prior(gamma(0.01, 0.01), class = phi), inits = 0,
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/dur_bic_mod_re")

summary(dur.bic.mod.re) #ns
print(summary(dur.bic.mod.re), digits = 4) 

# BIC (all values)
hist(master.class4.com$BIC_B_LAB, breaks = 40)
mean(master.class4.com$BIC_B_LAB) # 0.4
max(master.class4.com$BIC_B_LAB) #1
min(master.class4.com$BIC_B_LAB) #0

# BIC (only positive values)
bic_lab <- master.class4.com %>% filter(BIC_B_LAB > 0)
hist(bic_lab$BIC_B_LAB)
mean(bic_lab$BIC_B_LAB) # 0.77
max(bic_lab$BIC_B_LAB) # 1
min(bic_lab$BIC_B_LAB) # 0.16



## figure per site
points.bic <- master.class4.com %>% 
  select(SiteSubsitePlot, BIC_B_LAB) %>%
  filter(BIC_B_LAB > 0)

points.bic.site <- left_join(points.bic, site.only, by = "SiteSubsitePlot") %>% 
  left_join(., site.long, by = "SITE")

length(unique(points.bic.site$SITE)) #28
# 4 sites out because BIC_LAB = 0

length(unique(master.class4.com$SITE))
# ENDALEN, SADVENT, KANGER, WOLFCREEK
# also not in the map


# we need coords for each site
site.coords <- boreal_master_perm2 %>% 
  distinct(SITE, .keep_all = T) %>%
  select(SITE, LAT, LONG) %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))

# map of average cover change for each site
bic.blab.site.mean <- points.bic.site %>% group_by(SITE) %>% 
  summarise(MeanBIC = mean(BIC_B_LAB)) %>%
  left_join(., site.coords, by = "SITE")


# FIGURE 1b
# plot per site and plots, order per longitude
(bic.points <- ggplot() +
    geom_point(data = points.bic.site, aes(x = reorder(SITE, LONG), y = BIC_B_LAB), 
               alpha = 0.3, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_point(data = bic.blab.site.mean, aes(x = SITE, y = MeanBIC), 
               size = 5, pch = 21, stroke = 1.5, colour = "black") +
    xlab("Study area") + ylab("Borealization (BCI)\n") +
    theme(axis.text.x  = element_text(angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 20),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(face="bold", size = 20),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          legend.background = element_blank()))

ggplot2::ggsave(bic.points, 
                filename = "figures/points_blab_bic_new.png", 
                width = 30, height = 20, units = "cm")



# FIGURE 1a
# Convert coords to utm
site.coords.utm2 <- transform_coord(x = bic.blab.site.mean, lon = "LONG", lat = "LAT", 
                                   new.names = c("lon.utm", "lat.utm"), 
                                   verbose = FALSE, bind = TRUE)

# Plot map at site level - unique species per site (with labels)
(site.bic.blab.map <- basemap(limits = 50, land.col = "#c9c7c1",
                              grid.size = 1, grid.col = "#8f8f8f") + 
    geom_point(data = site.coords.utm2, 
               aes(x = lon.utm, y = lat.utm, colour = MeanBIC), 
               size = 15, alpha = 0.7) + 
    scale_colour_gradient(low = "#f2e81d", high = "#2c6611") +
    labs(colour='Mean BCI') +
    geom_label_repel(data = site.coords.utm2, aes(lon.utm, lat.utm, label = SITE), 
                     color = "black", box.padding = 2, 
                     segment.color = "black", segment.size = 1.2, 
                     fill = "white", label.size = 0.7,  
                     size = 9, max.iter = 150000, max.overlaps = 50) +
    theme(legend.title = element_text(size = 22), 
          legend.text = element_text(size = 22), 
          legend.position ="right", 
          legend.margin = margin(1,1,1,1),
          legend.box.margin = margin(10,10,10,10), 
          plot.margin = grid::unit(c(0,0,0,0), "mm"),
          panel.background = element_blank(),
          panel.ontop = FALSE))

ggsave(site.bic.blab.map, filename = "figures/map_blab_bic_new.png", 
       width = 30, height = 30, units = "cm")






## 2) BAI CALCULATION ----
bia.all <- master.class.com %>% 
  select(SiteSubsitePlotYear, SiteSubsitePlot, YEAR, Class, RelCover, Duration) %>%
  group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR),
         MaxYear = max(YEAR)) %>%
  filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>%
  ungroup() %>%
  group_by(SiteSubsitePlotYear, Class) %>% 
  mutate(CoverSum = sum(RelCover)) %>%
  ungroup() %>%
  distinct(SiteSubsitePlotYear, Class, .keep_all = T) %>%
  select(-RelCover) %>%
  pivot_wider(names_from = Class, values_from = CoverSum, values_fill = 0) %>%
  pivot_longer(cols = 7:11, names_to = "Class", values_to = "CoverSum") %>%
  mutate(Timepoint = case_when(YEAR == MinYear ~ "Start",
                               YEAR == MaxYear ~ "End")) %>%
  pivot_wider(names_from = Timepoint, values_from = CoverSum) %>%
  group_by(SiteSubsitePlot, Class, .keep_all = T) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>%
  ungroup() %>%
  mutate(CoverChange = (End - Start)/Duration) %>%
  select(SiteSubsitePlot, Class, CoverChange, Duration)

write.csv(bia.all, "data/bia-abund-startend.csv")

length(unique(bia.all$SiteSubsitePlot))


## ABUNDANCES / LAB+B ONLY
bia.labb <- master.class.com %>% 
  select(SiteSubsitePlotYear, SiteSubsitePlot, YEAR, Class, RelCover, Duration) %>%
  group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR),
         MaxYear = max(YEAR)) %>%
  filter(YEAR == min(YEAR) | YEAR == max(YEAR)) %>%
  ungroup() %>%
  filter(Class %in% c("Boreal", "Low_Arctic_Boreal")) %>%
  group_by(SiteSubsitePlotYear) %>% 
  mutate(B_LAB_CoverSum = sum(RelCover)) %>%
  ungroup() %>%
  distinct(SiteSubsitePlotYear, .keep_all = T) %>%
  select(-RelCover) %>%
  mutate(Timepoint = case_when(YEAR == MinYear ~ "Start",
                               YEAR == MaxYear ~ "End")) %>%
  pivot_wider(names_from = Timepoint, values_from = B_LAB_CoverSum) %>%
  group_by(SiteSubsitePlot, .keep_all = T) %>% 
  summarise_each(funs(first(.[!is.na(.)]))) %>%
  ungroup() %>%
# NAs introduced when the species was not present at the start or at the end
# there are only NAs in the Start and End columns so we are fine to do (.) below, we won't replace anything else
# not a problem above because an extra pivot_wider that includes values_fill = 0
  replace(is.na(.), 0) %>%
  mutate(CoverChange = (End - Start)/Duration) %>%
  select(SiteSubsitePlot, CoverChange, Duration)

write.csv(bia.labb, "data/bia-abund-lab-b-startend.csv")

length(unique(bia.labb$SiteSubsitePlot))

# this leaves 1000 plots, there are a total of 1137 unique plots
# the ones missing either a) never had B/LAB species at the start or end
# or b) they had B+LAB species in between timepoints, but not at the start and/or end 
missing.check <- plots.only %>% left_join(., bia.labb, by = "SiteSubsitePlot") %>% filter(is.na(CoverChange))
missing.vector <- missing.check$SiteSubsitePlot
missing.plots <- master.class.com %>% filter(SiteSubsitePlot %in% missing.vector)

# Fill in the remaining plots with 0
# Since this is a total of 0% cover change (even if they never had boreal species to begin with)
bia.labb2 <- plots.only %>%
  left_join(bia.labb, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0)

write.csv(bia.labb2, "data/bia-abund-lab-b-startend-full.csv")



# BIA (all values)
hist(bia.labb2$CoverChange, breaks = 40)
mean(bia.labb2$CoverChange) # -0.09
max(bia.labb2$CoverChange) # 5.79
min(bia.labb2$CoverChange) # -5.82

# average BIA (full range)
bia.labb.sb <- bia.labb2 %>% left_join(., subsite.only, by = "SiteSubsitePlot")

bia.int.only.mod <- brm(CoverChange ~ 1 + (1|SiteSubsite),
                                    data = bia.labb.sb, family = 'gaussian',
                                    iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/bia_int_only")
summary(bia.int.only.mod)


# BIA (only positive values)
bia.labb.pos <- bia.labb2 %>% filter(CoverChange > 0) #488 plots -> 42.9% plots experienced increases in LAB spps

mean(bia.labb.pos$CoverChange) # 0.93
max(bia.labb.pos$CoverChange) # 5.79
min(bia.labb.pos$CoverChange) # 0.007

# BIA (only negative values)
bia.labb.neg <- bia.labb2 %>% filter(CoverChange < 0) #501 plots -> 44.06% plots experienced decreases in LAB spps

# BIA (just zeroes to double-check)
bia.zero <- bia.labb2 %>% filter(CoverChange == 0) #148 plots -> 13% plots no change in LAB spps



# data prep for figures
bia.labb2

bia.blab.site <- bia.labb2 %>% 
  left_join(., site.only, by = "SiteSubsitePlot") %>% 
  filter(CoverChange > 0) %>%
  left_join(., site.long, by = "SITE") %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))

# map of average cover change for each site
bia.blab.site.mean <- bia.blab.site %>% group_by(SITE) %>% 
  summarise(MeanCoverChange = mean(CoverChange)) %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))


# FIGURE 1d
(bic.abun.blab.plot <- ggplot() +
    geom_point(data = bia.blab.site, aes(x = reorder(SITE, LONG), y = CoverChange),
               alpha = 0.3, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_point(data = bia.blab.site.mean, aes(x = SITE, y = MeanCoverChange), 
               size = 5, pch = 21, stroke = 1.5, colour = "black") +
    xlab("Study area") + ylab("Borealization (BAI)\n") +
    theme(axis.text.x  = element_text(angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          axis.text.y  = element_text(size = 15, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 20),
          axis.title.y = element_text(face="bold", size = 20),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          #legend.position = c(0.95, 0.88), legend.key = element_blank(),
          legend.background = element_blank()))

ggsave(bic.abun.blab.plot, filename = "figures/points_blab_bia_new.png", 
       width = 30, height = 20, units = "cm")



# join back in
bia.blab.site.mean.coord <- bia.blab.site.mean %>% left_join(site.coords, by = "SITE")

# Convert coords to utm
site.coords.utm <- transform_coord(x = bia.blab.site.mean.coord, lon = "LONG", lat = "LAT", 
                                     new.names = c("lon.utm", "lat.utm"), 
                                     verbose = FALSE, bind = TRUE)

# FIGURE 1c
(site.bia.blab.map <- basemap(limits = 50, land.col = "#c9c7c1",
                              grid.size = 1, grid.col = "#8f8f8f") + 
    geom_point(data = site.coords.utm, 
               aes(x = lon.utm, y = lat.utm, colour = MeanCoverChange), 
               size = 15, alpha = 0.7) + 
    scale_colour_gradient(low = "#f2e81d", high = "#2c6611") +
    labs(colour='Mean BAI') +
    geom_label_repel(data = site.coords.utm, aes(lon.utm, lat.utm, label = SITE), 
                    color = "black", box.padding = 2, 
                    segment.color = "black", segment.size = 1.2, 
                    fill = "white", label.size = 0.7,  
                    size = 9, max.iter = 150000, max.overlaps = 50) +
    theme(legend.title = element_text(size = 22), 
          legend.text=element_text(size = 22), 
          legend.position="right", 
          legend.margin=margin(1,1,1,1),
          legend.box.margin=margin(10,10,10,10), 
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.background = element_blank(),
          panel.ontop = FALSE))

ggsave(site.bia.blab.map, filename = "figures/map_blab_bia_new.png", 
       width = 30, height = 30, units = "cm")








## MODEL PREP ----
bia.labb2

# prepare treeline data
treeline <- read.csv("data/subs_distance_treeline_corrected_dec2024.csv") %>%
  select(SiteSubsite, DistanceTreelineCorrected, InterruptedCorrected, biome, ElevationSubsite) %>%
  mutate(DistanceTreelineCorrectedKm = DistanceTreelineCorrected/1000,
         DistanceTreelineCorrectedKmCentred = DistanceTreelineCorrectedKm - mean(DistanceTreelineCorrectedKm),
         DistanceTreelineLog = (log(DistanceTreelineCorrected+0.0001/1000)))

hist(treeline$DistanceTreelineCorrected, breaks = 100)
hist(treeline$DistanceTreelineCorrectedKmCentred, breaks = 100)
hist(treeline$DistanceTreelineLog, breaks = 100)

# join treeline with other relevant variables
load("data/boreal_master_perm.RData")

master_perm_site <- boreal_master_perm2 %>% 
  select(SiteSubsite, SiteSubsitePlot, Moisture, LAT, LONG, SurveyedArea, DomGrazer,
         GrazerIntensity, Region, Duration, permafrost) %>% 
  distinct(SiteSubsitePlot, .keep_all = T) %>% 
  left_join(., treeline, by = "SiteSubsite")


# join all variables with BAI data
master_lab_b <- master_perm_site %>% 
  left_join(., bia.labb2, by = "SiteSubsitePlot") %>%
  select(-Duration.y)

# this includes values for both B+LAB together, perfect Gaussian distribution
hist(master_lab_b$CoverChange)

# check structure of different variables
glimpse(master_lab_b)
str(master_lab_b)

# join with climate data
clima <- read.csv("data/climate_data.csv")

abun.data.clim <- left_join(master_lab_b, clima, by = c("SiteSubsitePlot"))



## Calculate start borealization

# calculate number of B+LAB species at t0
start.bor.spps <- master.class2.com %>% 
  group_by(SiteSubsitePlot) %>% 
  filter(YEAR == min(YEAR)) %>%
  filter(Class %in% c("Boreal", "Low_Arctic_Boreal")) %>%
  mutate(StartBorSpps = length(Class)) %>% ungroup() %>%
  distinct(SiteSubsitePlot, .keep_all = T) %>%
  select(SiteSubsitePlot, StartBorSpps)

# calculate starting cover of B+LAB species at t0
start.bor.cov <- master.class.com %>% 
  select(SiteSubsitePlotYear, SiteSubsitePlot, YEAR, Class, RelCover, Duration) %>%
  group_by(SiteSubsitePlot) %>% 
  filter(YEAR == min(YEAR)) %>%
  ungroup() %>%
  group_by(SiteSubsitePlotYear, Class) %>% 
  mutate(CoverSum = sum(RelCover)) %>%
  ungroup() %>%
  distinct(SiteSubsitePlotYear, Class, .keep_all = T) %>%
  select(-RelCover) %>%
  filter(Class %in% c("Boreal", "Low_Arctic_Boreal")) %>%
  group_by(SiteSubsitePlot) %>%
  mutate(StartBorAbun = sum(CoverSum)) %>%
  distinct(SiteSubsitePlot, .keep_all = T) %>%
  select(SiteSubsitePlot, StartBorAbun)

# join with main and replace NAs by 0
abun.data.clim2 <- left_join(abun.data.clim, start.bor.cov, by = "SiteSubsitePlot") %>%
  left_join(., start.bor.spps, by = "SiteSubsitePlot") %>%
  tidyr::replace_na(list(StartBorSpps = 0, StartBorAbun = 0))



## PAIRWISE CORRELATIONS ----

# Check correlations between all different drivers

# leave only relevant variables (all types)
pair.cor <- abun.data.clim2 %>% select(-c(1, 2, 5, 10, 12, 16, 18:23)) %>%
  rename(Latitude = LAT.x,
         DistanceTreeline = DistanceTreelineCorrectedKmCentred)

# convert categorical variables to numerical
pair.cor.all <- pair.cor %>% mutate(Moisture = case_when(Moisture == "Dry" ~ 0,
                                                         Moisture == "Moist" ~ 1,
                                                         Moisture == "Mixed" ~ 2,
                                                         Moisture == "Wet" ~ 3),
                                    DomGrazer = case_when(DomGrazer == "None" ~ 0,
                                                          DomGrazer == "Insects" ~ 1,
                                                          DomGrazer == "Small" ~ 2,
                                                          DomGrazer == "Birds" ~ 3,
                                                          DomGrazer == "Mixed" ~ 4,
                                                          DomGrazer == "Large" ~ 5),
                                    GrazerIntensity = case_when(GrazerIntensity == "Low" ~ 1,
                                                                GrazerIntensity == "Medium" ~ 2,
                                                                GrazerIntensity == "High" ~ 3),
                                  Permafrost = case_when(permafrost == "none" ~ 0,
                                                           permafrost == "sporadic" ~ 1,
                                                           permafrost == "continuous" ~ 2),
                                    Interrupted = case_when(InterruptedCorrected == "Uninterrupted" ~ 0,
                                                            InterruptedCorrected == "Small water bodies" ~ 1,
                                                            InterruptedCorrected == "Mountains" ~ 2,
                                                            InterruptedCorrected == "Large water bodies" ~ 3),
                                    Biome = case_when(biome == 6 ~ 0,
                                                      biome == 11 ~ 1),
                                    Region = case_when(Region == "Eurasia" ~ 1,
                                                       Region == "North America-West" ~ 2,
                                                       Region == "North America-East" ~ 3,
                                                       Region == "GreenIceLand" ~ 4)) %>%
  select(-c(InterruptedCorrected, biome, permafrost))

# use only pairwise-complete observations
pair.cor.all.matrix.complete = cor(pair.cor.all, method = "spearman", use = "pairwise.complete.obs")

# visualize correlogram (with numbers)
(pairwise.all.numbers <- corrplot(pair.cor.all.matrix.complete, method = "number", 
                                  type = 'lower', diag = FALSE, tl.col = 'black', tl.srt = 45))




## ABUNDANCE MODS ----

## B+LAB / INCREASES ONLY ##

# keep increases in cover only
abun.data.inc <- abun.data.clim2 %>% filter(CoverChange > 0)

mean(abun.data.inc$CoverChange)
max(abun.data.inc$CoverChange)
min(abun.data.inc$CoverChange)


#### 01) Biogeographic model (inc only) ####

# check correlations between variables
abun.data.biog <- abun.data.inc %>% 
  select(LAT.x, Region, biome, DistanceTreelineCorrectedKmCentred, InterruptedCorrected) %>%
  mutate(Biome = case_when(biome == 6 ~ 0, biome == 11 ~ 1),
         Interrupted = case_when(InterruptedCorrected == "Uninterrupted" ~ 0,
                                 InterruptedCorrected == "Small water bodies" ~ 1,
                                 InterruptedCorrected == "Mountains" ~ 2,
                                 InterruptedCorrected == "Large water bodies" ~ 3),
         Region = case_when(Region == "Eurasia" ~ 1,
                            Region == "North America-West" ~ 2,
                            Region == "North America-East" ~ 3,
                            Region == "GreenIceLand" ~ 4)) %>%
  rename(Latitude = LAT.x, DistTreeline = DistanceTreelineCorrectedKmCentred) %>%
  select(-c(InterruptedCorrected, biome))


# use only pairwise-complete observations
pair.cor.abun.mod1 = cor(abun.data.biog, method = "spearman")

# FIGURE S2b
# visualize correlogram with coloured square - save with 1000 width
(pairwise.abun.mod1.num.v2 <- corrplot(pair.cor.abun.mod1, type = 'lower', method = "color", 
                                    diag = FALSE, tl.col = 'black', tl.srt = 45, 
                                    addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.8,
                                    cl.cex = 1))

# we remove latitude (correlated with most), interrupted (correlated with most) and biome (same)
# we leave region and distance to treeline 


# final model (increases only)
abun_blab_biogeo_mod_dec2024 <- brm(bf(CoverChange ~ Region + 
                                     DistanceTreelineCorrectedKmCentred + 
                                     (1|SiteSubsite)),
                                data = abun.data.inc, family = 'gaussian',
                                iter = 2000, chains = 4, warmup = 400, 
                                file = "models/dec2024_abun_blab_biogeo_mod")

print(summary(abun_blab_biogeo_mod_dec2024), digits = 4) #dist negative
conditional_effects(abun_blab_biogeo_mod_dec2024) # Eurasia greater than and NAMW


# dist treeline predictions
abun_blab_treeline_pred <- ggpredict(abun_blab_biogeo_mod_dec2024, 
                                    terms = "DistanceTreelineCorrectedKmCentred [sample = 50]")
colnames(abun_blab_treeline_pred) = c('DistanceTreelineCorrectedKmCentred', 'fit', 'lwr', 'upr', 'dunno')

# Mean value used for centering 
mean(treeline$DistanceTreelineCorrectedKm) #342.2577

# Mean value of the predicted values for Y (fit)
mean(abun_blab_treeline_pred$fit) # 0.613618

# Back-centre 
abun_blab_treeline_pred_more <- abun_blab_treeline_pred %>% 
  mutate(treeline_raw = DistanceTreelineCorrectedKmCentred + 342.2577)



# region predictions
abun_blab_region_pred <- ggpredict(abun_blab_biogeo_mod_dec2024, 
                                     terms = "Region") %>% select(-group) 
colnames(abun_blab_region_pred) = c('Region', 'fit', 'lwr', 'upr')

abun_blab_region_pred2 <- abun_blab_region_pred %>% 
  mutate(Region = case_when(Region == "North America-West" ~ "WNA", 
                            Region == "North America-East" ~ "ENA", 
                            Region == "GreenIceLand" ~ "GI",
                            Region == "Eurasia" ~ "EA")) %>%
  mutate(Order = case_when(Region == "WNA" ~ 1, 
                            Region == "ENA" ~ 2, 
                            Region == "GI" ~ 3,
                           Region == "EA" ~ 4)) %>%
  mutate(ModelType = "Borealization only")


#### 01) Biogeographic model (full range) ####
abun_blab_biogeo_mod_all_dec2024 <- brm(bf(CoverChange ~ Region + 
                                     DistanceTreelineCorrectedKmCentred + 
                                     (1|SiteSubsite)),
                                data = abun.data.clim2, family = 'gaussian',
                                iter = 2000, chains = 4, warmup = 400, 
                                file = "models/dec2024_abun_blab_biogeo_mod_all_complete")

summary(abun_blab_biogeo_mod_all_dec2024) #distance ns
conditional_effects(abun_blab_biogeo_mod_all_dec2024) # region ns 


# region predictions (full range)
abun_blab_region_fullrange_pred <- ggpredict(abun_blab_biogeo_mod_all_dec2024, 
                                   terms = "Region") %>% select(-group) 

colnames(abun_blab_region_fullrange_pred) = c('Region', 'fit', 'lwr', 'upr')

abun_blab_region_fullrange_pred2 <- abun_blab_region_fullrange_pred %>% 
  mutate(Region = case_when(Region == "North America-West" ~ "WNA", 
                            Region == "North America-East" ~ "ENA", 
                            Region == "GreenIceLand" ~ "GI",
                            Region == "Eurasia" ~ "EA")) %>%
  mutate(Order = case_when(Region == "WNA" ~ 1, 
                           Region == "ENA" ~ 2, 
                           Region == "GI" ~ 3,
                           Region == "EA" ~ 4)) %>%
  mutate(ModelType = "Full range")


# merge the two datasets
abun_blab_two_mods <- bind_rows(abun_blab_region_pred2, abun_blab_region_fullrange_pred2)


# FIGURE 2h
(abun_blab_region_plot <- ggplot() +
    geom_errorbar(data = abun_blab_two_mods, aes(x = reorder(Region, Order), ymin = lwr, ymax = upr, fill = ModelType), 
                  size = 1, width = 0.3, position = position_dodge(0.5)) + 
    geom_point(data = abun_blab_two_mods, aes(x = reorder(Region, Order), y = fit, fill = ModelType), 
               size = 9, pch = 21, colour = "black", position = position_dodge(0.5)) + 
    scale_fill_manual(values = c("Borealization only" = "#2c6611", "Full range" = "#a19f9f")) +
    xlab("Biogeographic region") + 
    ylab("Borealization (BAI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position = "none",
          legend.key = element_blank(),
          legend.title = element_blank(), 
          legend.text = element_blank(),
          legend.background = element_blank()))

ggsave(abun_blab_region_plot, filename = "figures/dec2024_fig2_bia_region.png", 
       width = 20, height = 15, units = "cm")





# dist treeline predictions
abun_blab_treeline_all_pred <- ggpredict(abun_blab_biogeo_mod_all_dec2024, 
                                     terms = "DistanceTreelineCorrectedKmCentred [sample = 30]")
colnames(abun_blab_treeline_all_pred) = c('DistanceTreelineCorrectedKmCentred', 'fit', 'lwr', 'upr', 'dunno')


# Mean value used for centering 
mean(treeline$DistanceTreelineCorrectedKm) #342.2577

# Mean value of the predicted values for Y (fit)
mean(abun_blab_treeline_all_pred$fit) # -0.01941174

# Back-centre dist treeline values
abun_blab_treeline_all_pred_more <- abun_blab_treeline_all_pred %>% 
  mutate(treeline_raw = DistanceTreelineCorrectedKmCentred + 342.2577)
#mutate(y_constant = 341.5876+fit)

abun.data.neg.xx <- abun.data.clim2 %>% filter(CoverChange <= 0)

# FIGURE 2i
(abun_blab_treeline_plot_backt <- ggplot() +
    geom_point(data = abun.data.clim2, aes(x = DistanceTreelineCorrectedKm, y = CoverChange), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = abun.data.neg.xx, aes(x = DistanceTreelineCorrectedKm, y = CoverChange), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = abun_blab_treeline_pred_more, aes(x = treeline_raw, y = fit), colour = "#2c6611", linetype = "solid", linewidth = 2) +
    geom_ribbon(data = abun_blab_treeline_pred_more, aes(x = treeline_raw, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = abun_blab_treeline_all_pred_more, aes(x = treeline_raw, y = fit), colour = "black", linetype = "dashed", linewidth = 2) +
    geom_ribbon(data = abun_blab_treeline_all_pred_more, aes(x = treeline_raw, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Distance to treeline (km)") + 
    ylab("Borealization (BAI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(abun_blab_treeline_plot_backt, filename = "figures/dec2024_fig2_bia_treeline.png", 
       width = 20, height = 15, units = "cm")








#### 02) Climatic model (inc only) ####
abun.data.clim.mod.data <- abun.data.inc %>% 
  select(WarmQSlope, PrecSlope, MinTempSlope, tmin_clim, warmq_clim, prec_clim)

# use only pairwise-complete observations
pair.cor.abun.mod2 = cor(abun.data.clim.mod.data, method = "spearman", use = "pairwise.complete.obs")

# FIGURE S2d
# visualize correlogram with coloured square - save with 1000 width
(pairwise.abun.mod2.num.v2 <- corrplot(pair.cor.abun.mod2, type = 'lower', method = "color", 
                                       diag = FALSE, tl.col = 'black', tl.srt = 45, 
                                       addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.8,
                                       cl.cex = 1))

# tmin_clim related to prec_clim, we remove tmin_clim as prec_clim is the only othe prec variable

# model
abun_blab_clim_mod_dec2024 <- brm(bf(CoverChange ~ WarmQSlope + PrecSlope + MinTempSlope +
                            warmq_clim + prec_clim + (1|SiteSubsite)),
                       data = abun.data.inc, family = 'gaussian',
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/dec2024_abun_blab_clim_mod")

summary(abun_blab_clim_mod_dec2024) # warm q slope ns, prec slope negative, all others ns
conditional_effects(abun_blab_clim_mod_dec2024)
# where there has been greater increases in prec over time, there has been less increases in
# abundance of B+LAB species

# warmq predictions
abun_blab_clim_pred <- ggpredict(abun_blab_clim_mod_dec2024, 
                                 terms = "WarmQSlope [sample = 30]")
colnames(abun_blab_clim_pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'dunno')




#### 02) Climatic model (full range) ####

# there is an outlier value for INCLINE_SKJ_7_1. After consulting with the site PIs, it was decided
# to remove this value as it does not align with site-specific knowledge and it's very possibly caused by
# climatic data having issues in mountain/fjord areas. The rest of climatic values for the site are fine.
abun.data.clim.fixed <- abun.data.clim2 %>% 
  mutate(PrecSlope = case_when(SiteSubsitePlot == "INCLINE_SKJ:SKJ:SKJ_7_1" ~ NA,
                               TRUE ~ PrecSlope))

# model with all
abun_blab_clim_mod_all_dec2024 <- brm(bf(CoverChange ~ WarmQSlope + PrecSlope + MinTempSlope +
                                   warmq_clim + prec_clim + (1|SiteSubsite)),
                              data = abun.data.clim.fixed, family = 'gaussian',
                              iter = 2000, chains = 4, warmup = 400, 
                              file = "models/dec2024_abun_blab_clim_mod_all_complete")

summary(abun_blab_clim_mod_all_dec2024) # warm q slope ns, prec slope ns, min temp slope ns, temp clim ns, prec clim ns
conditional_effects(abun_blab_clim_mod_all_dec2024)



# prec change predictions (inc only)
abun_blab_prec_pred <- ggpredict(abun_blab_clim_mod_dec2024, 
                                 terms = "PrecSlope [sample = 50]")
colnames(abun_blab_prec_pred) = c('PrecSlope', 'fit', 'lwr', 'upr', 'dunno')


# prec change predictions (all)
abun_blab_prec_fullrange_pred <- ggpredict(abun_blab_clim_mod_all_dec2024, 
                                 terms = "PrecSlope [sample = 60]")
colnames(abun_blab_prec_fullrange_pred) = c('PrecSlope', 'fit', 'lwr', 'upr', 'dunno')

# FIGURE 2j
(abun_blab_precslope_plot <- ggplot() +
    geom_point(data = abun.data.inc, aes(x = PrecSlope, y = CoverChange), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = abun.data.neg, aes(x = PrecSlope, y = CoverChange), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = abun_blab_prec_pred, aes(x = PrecSlope, y = fit), colour = "#2c6611", linewidth = 2) +
    geom_ribbon(data = abun_blab_prec_pred, aes(x = PrecSlope, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = abun_blab_prec_fullrange_pred, aes(x = PrecSlope, y = fit), colour = "black", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = abun_blab_prec_fullrange_pred, aes(x = PrecSlope, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Precipitation change (mm per year)") + ylab("Borealization (BAI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(abun_blab_precslope_plot, filename = "figures/dec2024_fig2_bia_precch_new.png", 
       width = 20, height = 15, units = "cm")







#### 03) Local model (inc only) ####
abun.data.local <- abun.data.inc %>% 
  select(Moisture, DomGrazer, GrazerIntensity, permafrost, 
         SurveyedArea, ElevationSubsite, StartBorAbun)

# convert categorical variables to numerical for correlogram
abun.data.local.num <- abun.data.local %>% mutate(Moisture = case_when(Moisture == "Dry" ~ 0,
                                                         Moisture == "Moist" ~ 1,
                                                         Moisture == "Mixed" ~ 2,
                                                         Moisture == "Wet" ~ 3),
                                    DomGrazer = case_when(DomGrazer == "None" ~ 0,
                                                          DomGrazer == "Insects" ~ 1,
                                                          DomGrazer == "Small" ~ 2,
                                                          DomGrazer == "Birds" ~ 3,
                                                          DomGrazer == "Mixed" ~ 4,
                                                          DomGrazer == "Large" ~ 5),
                                    GrazerIntensity = case_when(GrazerIntensity == "Low" ~ 1,
                                                                GrazerIntensity == "Medium" ~ 2,
                                                                GrazerIntensity == "High" ~ 3),
                                    permafrost = case_when(permafrost == "none" ~ 0,
                                                           permafrost == "sporadic" ~ 1,
                                                           permafrost == "continuous" ~ 2))

abun.data.local.num2 <- abun.data.local.num %>% rename(Elevation = ElevationSubsite,
                                                       GrazerInt = GrazerIntensity,
                                                       Permafrost = permafrost,
                                                       PlotSize = SurveyedArea,
                                                       InitialBor = StartBorAbun)

# use only pairwise-complete observations
pair.cor.abun.mod3 = cor(abun.data.local.num2, method = "spearman", use="pairwise.complete.obs")

# FIGURE S2f
# visualize correlogram
(pairwise.abun.mod3.num <- corrplot(pair.cor.abun.mod3, type = 'lower', method = "color", 
                                    diag = FALSE, tl.col = 'black', tl.srt = 45,
                                    addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.5,
                                    cl.cex = 1))
# nothing majorly correlated





# run model
#abun_blab_local_mod <- brm(bf(CoverChange ~ Moisture + GrazerIntensity + 
#                                DomGrazer + ElevationSubsite + permafrost + 
#                                SurveyedArea + StartBorAbun + (1|SiteSubsite)),
#                        data = abun.data.inc, family = 'gaussian',
#                        iter = 2000, chains = 4, warmup = 400, 
#                        file = "models/dec2024_abun_blab_local_mod")

#summary(abun_blab_local_mod) 
#conditional_effects(abun_blab_local_mod)

# this model has convergence issues: moisture_mixed and permafrost_discontinuous have a very small sample size
# they both belong to DISKO:DISTURBANCE (all plots) so let's remove this for this mod
abun.data.inc2 <- abun.data.inc %>% filter(SiteSubsite %notin% "DISKO:DISTURBANCE")

# run model
abun_blab_local_mod2 <- brm(bf(CoverChange ~ Moisture + GrazerIntensity + 
                                DomGrazer + ElevationSubsite + permafrost + 
                                SurveyedArea + StartBorAbun + (1|SiteSubsite)),
                           data = abun.data.inc2, family = 'gaussian',
                           iter = 2000, chains = 4, warmup = 400, 
                           file = "models/dec2024_abun_blab_local_mod2")

summary(abun_blab_local_mod2) # start bor negative
conditional_effects(abun_blab_local_mod2)
# where the starting B+LAB cover was greater, there were fewer increases in abundance



# start bor predictions
abun_blab_local_pred <- ggpredict(abun_blab_local_mod2, 
                                 terms = "StartBorAbun [0:100]")
colnames(abun_blab_local_pred) = c('StartBorAbun', 'fit', 'lwr', 'upr', 'dunno')




#### 03) Local model (all) ####
abun.data.all2 <- abun.data.clim2 %>% filter(SiteSubsite %notin% "DISKO:DISTURBANCE")

# run model
abun_blab_local_mod_all <- brm(bf(CoverChange ~ Moisture + GrazerIntensity + 
                                    DomGrazer + ElevationSubsite + permafrost + 
                                    SurveyedArea + StartBorAbun + (1|SiteSubsite)),
                               data = abun.data.all2, family = 'gaussian',
                               iter = 2000, chains = 4, warmup = 400, 
                               file = "models/dec2024_abun_blab_local_mod_all_complete")

print(summary(abun_blab_local_mod_all), digits = 4) # start bor negative, elevation positive
conditional_effects(abun_blab_local_mod_all)



# start bor predictions (full range model)
abun_blab_local_fullrange_pred <- ggpredict(abun_blab_local_mod_all, 
                                  terms = "StartBorAbun [0:100]")
colnames(abun_blab_local_fullrange_pred) = c('StartBorAbun', 'fit', 'lwr', 'upr', 'dunno')

# separate negative and zero values
abun.data.local.neg <- abun.data.all2 %>% filter(CoverChange <=0)


# FIGURE 2l 
(abun_blab_local_plot <- ggplot() +
    geom_point(data = abun.data.inc2, aes(x = StartBorAbun, y = CoverChange), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = abun.data.local.neg, aes(x = StartBorAbun, y = CoverChange), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = abun_blab_local_fullrange_pred, aes(x = StartBorAbun, y = fit), colour = "black", linewidth = 2) +
    geom_ribbon(data = abun_blab_local_fullrange_pred, aes(x = StartBorAbun, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    geom_line(data = abun_blab_local_pred, aes(x = StartBorAbun, y = fit), colour = "#2c6611", linewidth = 2) +
    geom_ribbon(data = abun_blab_local_pred, aes(x = StartBorAbun, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    xlab("Initial boreal status (%)") + ylab("Borealization (BAI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          #legend.position = c(0.87, 0.89), legend.key = element_blank(),
          legend.background = element_blank()))

ggsave(abun_blab_local_plot, filename = "figures/dec2024_fig2_bia_startborabun_new.png", 
       width = 20, height = 15, units = "cm")




# elevation predictions (inc only)
abun_blab_elev_pred <- ggpredict(abun_blab_local_mod2, terms = "ElevationSubsite [sample = 45]")
colnames(abun_blab_elev_pred) = c('ElevationSubsite', 'fit', 'lwr', 'upr', 'dunno')


# elevation predictions (full range)
abun_blab_elev_fullrange_pred <- ggpredict(abun_blab_local_mod_all, terms = "ElevationSubsite [sample = 47]")
colnames(abun_blab_elev_fullrange_pred) = c('ElevationSubsite', 'fit', 'lwr', 'upr', 'dunno')


# FIGURE 2k
(abun_blab_elev_plot <- ggplot() +
    geom_point(data = abun.data.inc2, aes(x = ElevationSubsite, y = CoverChange),
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") +
    geom_point(data = abun.data.local.neg, aes(x = ElevationSubsite, y = CoverChange),
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") +
    geom_line(data = abun_blab_elev_pred, aes(x = ElevationSubsite, y = fit), colour = "#2c6611", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = abun_blab_elev_pred, aes(x = ElevationSubsite, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = abun_blab_elev_fullrange_pred, aes(x = ElevationSubsite, y = fit), colour = "black", linetype = "solid", linewidth = 2) +
    geom_ribbon(data = abun_blab_elev_fullrange_pred, aes(x = ElevationSubsite, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Elevation (m a.s.l.)") + ylab("Borealization (BAI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.background = element_blank()))

ggsave(abun_blab_elev_plot, filename = "figures/dec2024_fig2_bia_elev_plot.png",
       width = 20, height = 15, units = "cm")





## COLONISERS MODS ----
master.class4.com # bic
master_perm_site # covariates
clima #climate
start.bor.spps # 

# join all datasets together
colo.data.full <- master.class4.com %>% 
  select(SiteSubsitePlot, BIC_B_LAB) %>%
  left_join(., master_perm_site, by = "SiteSubsitePlot") %>%
  left_join(., clima, by = "SiteSubsitePlot") %>%
  left_join(., start.bor.spps, by = "SiteSubsitePlot") %>%
  left_join(., start.bor.cov, by = "SiteSubsitePlot") %>%
  tidyr::replace_na(list(StartBorSpps = 0))

# check distribution
hist(colo.data.full$BIC_B_LAB) # zero-one inflated beta

# all correct now
str(colo.data.full)


## B+LAB / INCREASES

# keep increases in cover only
col.data.inc <- colo.data.full %>% filter(BIC_B_LAB > 0)
hist(col.data.inc$BIC_B_LAB)

mean(col.data.inc$BIC_B_LAB) # 0.77
min(col.data.inc$BIC_B_LAB) # 0.16
max(col.data.inc$BIC_B_LAB) # 1


#### 01) Biogeographic model (inc only) ####

# check correlations between variables
colo.data.biog <- col.data.inc %>% 
  select(LAT.x, Region, biome, DistanceTreelineCorrectedKmCentred, InterruptedCorrected) %>%
  mutate(Biome = case_when(biome == 6 ~ 0, biome == 11 ~ 1),
         Interrupted = case_when(InterruptedCorrected == "Uninterrupted" ~ 0,
                                 InterruptedCorrected == "Small water bodies" ~ 1,
                                 InterruptedCorrected == "Mountains" ~ 2,
                                 InterruptedCorrected == "Large water bodies" ~ 3),
         Region = case_when(Region == "Eurasia" ~ 1,
                            Region == "North America-West" ~ 2,
                            Region == "North America-East" ~ 3,
                            Region == "GreenIceLand" ~ 4)) %>%
  rename(Latitude = LAT.x, DistTreeline = DistanceTreelineCorrectedKmCentred) %>%
  select(-c(InterruptedCorrected, biome))



# use only pairwise-complete observations
pair.cor.col.mod1 = cor(colo.data.biog, method = "spearman")

# FIGURE S2a
# visualize correlogram
(pairwise.col.mod1.num <- corrplot(pair.cor.col.mod1, type = 'lower', method = "color", 
                                    diag = FALSE, tl.col = 'black', tl.srt = 45,
                                   addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.8,
                                   cl.cex = 1))

# we remove latitude (correlated with biome and interrupted) 
# we remove interrupted (correlated with latitude, dist treeline and biome)
# we leave region, distance to treeline and biome. 

# correct data so it fits a beta distribution
col.data.inc2 <- col.data.inc %>% 
  mutate(BIC_B_LAB_corrected = BIC_B_LAB - 0.0001) %>%
  mutate(biome = case_when(biome == 6 ~ "Alpine", biome == 11 ~ "Arctic"))


# model
col_blab_biogeo_mod <- brm(bf(BIC_B_LAB_corrected ~ DistanceTreelineCorrectedKmCentred +
                                biome + Region + (1|SiteSubsite)),
                            data = col.data.inc2, family = 'beta',
                            init = 0,
                            iter = 2000, chains = 4, warmup = 400, 
                            file = "models/dec2024_col_blab_biogeo_mod")

print(summary(col_blab_biogeo_mod), digits = 5) # all ns
conditional_effects(col_blab_biogeo_mod) # Eurasia > GreenIceland


# dist treeline predictions (col only)
col_blab_treeline_pred <- ggpredict(col_blab_biogeo_mod, 
                                    terms = "DistanceTreelineCorrectedKmCentred [sample = 50]")
colnames(col_blab_treeline_pred) = c('DistanceTreelineCorrectedKmCentred', 'fit', 'lwr', 'upr', 'dunno')




#### 01) Biogoegraphic model (full range) ####
colo.data.full 

# correct biome and offset from ZOIB for a zero-inflated distribution
col.data.full2 <- colo.data.full %>% 
  mutate(biome = case_when(biome == 6 ~ "Alpine", biome == 11 ~ "Arctic")) %>%
  mutate(BIC_B_LAB_corrected = case_when(BIC_B_LAB > 0 ~ BIC_B_LAB - 0.0001,
                                         BIC_B_LAB == 0 ~ 0))


## BINOMIAL COUNTS (for binomial model)

# calculate the Borealization Colonization Index
bic.bin <- master.class2.com %>% 
  filter(trend == "Coloniser") %>%
  group_by(SiteSubsitePlot, Class) %>% 
  summarise(n = n()) %>%
  #left_join(., duration.only, by = "SiteSubsitePlot") %>%
  ungroup() %>%
  group_by(SiteSubsitePlot) %>%
  mutate(TotalColonisers = sum(n)) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  left_join(., duration.only, by = "SiteSubsitePlot") %>%
  mutate(BorealCount = Boreal + Low_Arctic_Boreal)


# plots without colonisers (boreal or otherwise) are missing from this dataset - 
# need to fill them up with 0
# as they still count as zero borealization (even if they had no colonisations at all)
bic.bin2 <- plots.only %>% 
  left_join(bic.bin, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0) %>%
  #left_join(., duration.only, by = "SiteSubsitePlot") %>%
  left_join(., subsite.only, by = "SiteSubsitePlot") %>%
  left_join(., site.only, by = "SiteSubsitePlot")

bic.counts <- bic.bin2 %>% select(SiteSubsitePlot, TotalColonisers, BorealCount)

# fit model as binomial
col.data.full.counts <- col.data.full2 %>% left_join(., bic.counts, by = "SiteSubsitePlot")

col_blab_biogeo_mod_all_bino <- brm(BorealCount | trials(TotalColonisers) ~ 
                                      DistanceTreelineCorrectedKmCentred + biome + 
                                      Region + (1|SiteSubsite), 
                                    family = binomial("logit"),
                                    data = col.data.full.counts, 
                                    init = 0,
                                    iter = 2000, chains = 4, warmup = 400, 
                                    file = "models/dec2024_col_blab_biogeo_all_mod_bino")

print(summary(col_blab_biogeo_mod_all_bino), digits = 5) # dist treeline negative
conditional_effects(col_blab_biogeo_mod_all_bino) #ns 


## DISTANCE TO TREELINE 

# dist treeline predictions (full range binomial)
col_blab_treeline_fullrange_counts_pred <- ggpredict(col_blab_biogeo_mod_all_bino, 
                                                     terms = "DistanceTreelineCorrectedKmCentred [sample = 60]")
colnames(col_blab_treeline_fullrange_counts_pred) = c('DistanceTreelineCorrectedKmCentred', 'fit', 'lwr', 'upr', 'dunno')

# Mean value used for centering 
mean(treeline$DistanceTreelineCorrectedKm) #342.2577

# Mean value of the predicted values for Y (fit)
mean(col_blab_treeline_fullrange_counts_pred$fit) # 0.443882

# Back-centre dist treeline values
col_blab_treeline_pred_more <- col_blab_treeline_pred %>% 
  mutate(treeline_raw = DistanceTreelineCorrectedKmCentred + 342.2577)

# Back-centre (full range)
col_blab_treeline_pred_fullrange_counts_more <- col_blab_treeline_fullrange_counts_pred %>% 
  mutate(treeline_raw = DistanceTreelineCorrectedKmCentred + 342.2577)

colo.data.zero <- colo.data.full %>% filter(BIC_B_LAB == 0) 

# FIGURE 2b
(col_blab_treeline_plot_backt <- ggplot() +
    geom_point(data = col.data.inc2, aes(x = DistanceTreelineCorrectedKm, y = BIC_B_LAB), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = colo.data.zero, aes(x = DistanceTreelineCorrectedKm, y = BIC_B_LAB), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = col_blab_treeline_pred_more, aes(x = treeline_raw, y = fit), colour = "#2c6611", linetype = "dashed", linewidth = 2) +
    geom_ribbon(data = col_blab_treeline_pred_more, aes(x = treeline_raw, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_treeline_pred_fullrange_counts_more, aes(x = treeline_raw, y = fit), colour = "black", linetype = "solid", linewidth = 2) +
    geom_ribbon(data = col_blab_treeline_pred_fullrange_counts_more, aes(x = treeline_raw, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Distance to treeline (km)") + 
    ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(col_blab_treeline_plot_backt, filename = "figures/revised_fig2_bic_treeline.png", 
       width = 20, height = 15, units = "cm")


## REGION

# region predictions (inc only)
col_blab_region_pred <- ggpredict(col_blab_biogeo_mod, 
                                  terms = "Region") %>% select(-group) 
colnames(col_blab_region_pred) = c('Region', 'fit', 'lwr', 'upr')

col_blab_region_pred2 <- col_blab_region_pred %>% 
  mutate(Region = case_when(Region == "North America-West" ~ "WNA", 
                            Region == "North America-East" ~ "ENA", 
                            Region == "GreenIceLand" ~ "GI",
                            Region == "Eurasia" ~ "EA")) %>%
  mutate(Order = case_when(Region == "WNA" ~ 1, 
                           Region == "ENA" ~ 2, 
                           Region == "GI" ~ 3,
                           Region == "EA" ~ 4)) %>%
  mutate(ModelType = "Borealization only")



# region predictions (full range)
col_blab_region_fullrange_bino_pred <- ggpredict(col_blab_biogeo_mod_all_bino, 
                                                 terms = "Region") %>% select(-group) 

colnames(col_blab_region_fullrange_bino_pred) = c('Region', 'fit', 'lwr', 'upr')

col_blab_region_fullrange_bino_pred2 <- col_blab_region_fullrange_bino_pred %>% 
  mutate(Region = case_when(Region == "North America-West" ~ "WNA", 
                            Region == "North America-East" ~ "ENA", 
                            Region == "GreenIceLand" ~ "GI",
                            Region == "Eurasia" ~ "EA")) %>%
  mutate(Order = case_when(Region == "WNA" ~ 1, 
                           Region == "ENA" ~ 2, 
                           Region == "GI" ~ 3,
                           Region == "EA" ~ 4)) %>%
  mutate(ModelType = "Full range")


# merge the two datasets
col_blab_two_mods_bino <- bind_rows(col_blab_region_pred2, col_blab_region_fullrange_bino_pred2)



# FIGURE 2a
(col_blab_region_plot <- ggplot() +
    geom_errorbar(data = col_blab_two_mods_bino, aes(x = reorder(Region, Order), ymin = lwr, ymax = upr, fill = ModelType), 
                  size = 1, width = 0.3, position = position_dodge(0.5)) + 
    geom_point(data = col_blab_two_mods_bino, aes(x = reorder(Region, Order), y = fit, fill = ModelType), 
               size = 9, pch = 21, colour = "black", position = position_dodge(0.5)) + 
    scale_fill_manual(values = c("Borealization only" = "#2c6611", "Full range" = "#a19f9f")) +
    xlab("Biogeographic region") + 
    ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position = "none",
          legend.key = element_blank(),
          legend.title = element_blank(), 
          legend.text = element_blank(),
          legend.background = element_blank()))

ggsave(col_blab_region_plot, filename = "figures/revised_fig2_bic_region.png", 
       width = 20, height = 15, units = "cm")






#### 02) Climatic model (inc only) ####
col.data.clim <- col.data.inc2 %>% 
  select(WarmQSlope, PrecSlope, MinTempSlope, tmin_clim, warmq_clim, prec_clim)

# use only pairwise-complete observations
pair.cor.col.mod2 = cor(col.data.clim, method = "spearman", use = "pairwise.complete.obs")

# FIGURE S2c
# visualize correlogram (save to 1000 width)
(pairwise.col.mod2.num <- corrplot(pair.cor.col.mod2, type = 'lower', method = "color", 
                                    diag = FALSE, tl.col = 'black', tl.srt = 45,
                                   addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.8,
                                   cl.cex = 1))

# tmin_clim related to prec_clim, we remove tmin_clim as prec_clim is the only othe prec variable



# model
col_blab_clim_mod <- brm(bf(BIC_B_LAB_corrected ~ WarmQSlope + PrecSlope + MinTempSlope +
                               warmq_clim + prec_clim + (1|SiteSubsite)),
                          data = col.data.inc2, family = 'beta', init = 0,
                          file = "models/dec2024_col_blab_clim_mod")

print(summary(col_blab_clim_mod), digits = 4) # warmq slope negative, prec slope negative, temp positive
conditional_effects(col_blab_clim_mod)


# temp climatology predictions
col_blab_warmclim_pred <- ggpredict(col_blab_clim_mod, 
                                  terms = "warmq_clim [sample = 45]")
colnames(col_blab_warmclim_pred) = c('warmq_clim', 'fit', 'lwr', 'upr', 'dunno')




#### 02) Climatic model (full range) ####

# remove the INCLINE outlier as it's clearly erroneous
col.data.clim.fixed <- col.data.full2 %>% 
  mutate(PrecSlope = case_when(SiteSubsitePlot == "INCLINE_SKJ:SKJ:SKJ_7_1" ~ NA,
                               TRUE ~ PrecSlope))

# try removing the INCLINE plot instead of putting it as NA
#col.data.clim.fixed.xxx <- col.data.full2 %>% 
#  filter(SiteSubsitePlot != "INCLINE_SKJ:SKJ:SKJ_7_1")

col.data.clim.fixed.counts <- col.data.clim.fixed %>% left_join(., bic.counts, by = "SiteSubsitePlot")

col_blab_clim_all_mod_bino <- brm(BorealCount | trials(TotalColonisers) ~ 
                                    WarmQSlope + PrecSlope + MinTempSlope +
                                    warmq_clim + prec_clim + (1|SiteSubsite), 
                                  data = col.data.clim.fixed.counts, family = binomial("logit"),
                                  init = 0, iter = 2000, chains = 4, warmup = 400, 
                                  file = "models/dec2024_col_blab_clim_all_mod_bino")


print(summary(col_blab_clim_all_mod_bino), digits = 4) # warmq clim positive, prec clim positive, warm clim positive, prec clim positive
# prec slope ns (negative), temp slope ns (negative)
conditional_effects(col_blab_clim_all_mod_bino)
pp_check(col_blab_clim_all_mod_bino)



## PREC CLIMATOLOGY


# temp climatology predictions (increase only)
col_blab_precclim_pred <- ggpredict(col_blab_clim_mod, terms = "prec_clim [sample = 50]")
colnames(col_blab_precclim_pred) = c('prec_clim', 'fit', 'lwr', 'upr', 'dunno')

# temp climatology predictions (full range)
col_blab_precclim_fullrange_bino_pred <- ggpredict(col_blab_clim_all_mod_bino, 
                                                   terms = "prec_clim [sample = 70]")
colnames(col_blab_precclim_fullrange_bino_pred) = c('prec_clim', 'fit', 'lwr', 'upr', 'dunno')


# FIGURE 2d
(col_blab_precclim_plot <- ggplot() +
    geom_point(data = col.data.inc, aes(x = prec_clim, y = BIC_B_LAB), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = colo.data.zero, aes(x = prec_clim, y = BIC_B_LAB), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = col_blab_precclim_pred, aes(x = prec_clim, y = fit), colour = "#2c6611", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = col_blab_precclim_pred, aes(x = prec_clim, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_precclim_fullrange_bino_pred, aes(x = prec_clim, y = fit), colour = "black", linewidth = 2) +
    geom_ribbon(data = col_blab_precclim_fullrange_bino_pred, aes(x = prec_clim, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Precipitation (mm)") + ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(col_blab_precclim_plot, filename = "figures/revised_fig2_bic_precclim.png", 
       width = 20, height = 15, units = "cm")




## TEMP CLIMATOLOGY

# temp climatology predictions (full range)
col_blab_warmclim_fullrange_bino_pred <- ggpredict(col_blab_clim_all_mod_bino, 
                                                   terms = "warmq_clim [sample = 70]")
colnames(col_blab_warmclim_fullrange_bino_pred) = c('warmq_clim', 'fit', 'lwr', 'upr', 'dunno')


# FIGURE 2c
(col_blab_warmclim_plot <- ggplot() +
    geom_point(data = col.data.inc, aes(x = warmq_clim, y = BIC_B_LAB), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = colo.data.zero, aes(x = warmq_clim, y = BIC_B_LAB), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = col_blab_warmclim_pred, aes(x = warmq_clim, y = fit), colour = "#2c6611", linewidth = 2) +
    geom_ribbon(data = col_blab_warmclim_pred, aes(x = warmq_clim, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_warmclim_fullrange_bino_pred, aes(x = warmq_clim, y = fit), colour = "black", linewidth = 2) +
    geom_ribbon(data = col_blab_warmclim_fullrange_bino_pred, aes(x = warmq_clim, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Mean Annual Temperature (°C)") + ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(col_blab_warmclim_plot, filename = "figures/revised_fig2_bic_warmclim.png", 
       width = 20, height = 15, units = "cm")




## TEMP SLOPE

# temp slope predictions
col_blab_warmqs_pred <- ggpredict(col_blab_clim_mod, 
                                  terms = "WarmQSlope [sample = 45]")
colnames(col_blab_warmqs_pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'dunno')


# temp slope predictions (full range)
col_blab_warmqs_fullrange_bino_pred <- ggpredict(col_blab_clim_all_mod_bino, 
                                                 terms = "WarmQSlope [sample = 60]")
colnames(col_blab_warmqs_fullrange_bino_pred) = c('WarmQSlope', 'fit', 'lwr', 'upr', 'dunno')


# FIGURE 2e
(col_blab_warmslope_plot <- ggplot() +
    geom_point(data = col.data.inc, aes(x = WarmQSlope, y = BIC_B_LAB), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = colo.data.zero, aes(x = WarmQSlope, y = BIC_B_LAB), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = col_blab_warmqs_pred, aes(x = WarmQSlope, y = fit), colour = "#2c6611", linewidth = 2) +
    geom_ribbon(data = col_blab_warmqs_pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_warmqs_fullrange_bino_pred, aes(x = WarmQSlope, y = fit), colour = "black", linetype = "dashed", linewidth = 2) +
    geom_ribbon(data = col_blab_warmqs_fullrange_bino_pred, aes(x = WarmQSlope, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Temperature change (°C per year)") + ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(col_blab_warmslope_plot, filename = "figures/revised_fig2_bic_warmqslope.png", 
       width = 20, height = 15, units = "cm")




## PRECIPITATION SLOPE

# prec change predictions (full range)
col_bino_pred <- ggpredict(col_blab_clim_all_mod_bino, 
                           terms = "PrecSlope [sample = 50]")
colnames(col_bino_pred) = c('PrecSlope', 'fit', 'lwr', 'upr', 'dunno')


# prec slope predictions (increase only)
col_blab_precs_pred <- ggpredict(col_blab_clim_mod, 
                                 terms = "PrecSlope [sample = 55]")
colnames(col_blab_precs_pred) = c('PrecSlope', 'fit', 'lwr', 'upr', 'dunno')


# remove the outlier value
colo.data.inc.fixed <- col.data.inc %>% 
  mutate(PrecSlope = case_when(SiteSubsitePlot == "INCLINE_SKJ:SKJ:SKJ_7_1" ~ NA,
                               TRUE ~ PrecSlope))


# FIGURE 2f
(col_blab_precs_bino_plot <- ggplot() +
    geom_point(data = colo.data.inc.fixed, aes(x = PrecSlope, y = BIC_B_LAB), 
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") + 
    geom_point(data = colo.data.zero, aes(x = PrecSlope, y = BIC_B_LAB), 
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_line(data = col_blab_precs_pred, aes(x = PrecSlope, y = fit), colour = "#2c6611", linewidth = 2) +
    geom_ribbon(data = col_blab_precs_pred, aes(x = PrecSlope, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_bino_pred, aes(x = PrecSlope, y = fit), colour = "black", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = col_bino_pred, aes(x = PrecSlope, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Precipitation change (mm per year)") + ylab("Borealization (BCI)") +
    xlim(c(-1.3, 11.5)) +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))

ggsave(col_blab_precs_bino_plot, filename = "figures/revised_fig2_bic_precslope_bino.png", 
       width = 20, height = 15, units = "cm")


#### 03) Local model (inc only) ####
col.data.local <- col.data.inc2 %>% 
  select(Moisture, DomGrazer, GrazerIntensity, permafrost, 
         SurveyedArea, ElevationSubsite, StartBorSpps)

# convert categorical variables to numerical for correlogram
col.data.local.num <- col.data.local %>% mutate(Moisture = case_when(Moisture == "Dry" ~ 0,
                                                                       Moisture == "Moist" ~ 1,
                                                                       Moisture == "Mixed" ~ 2,
                                                                       Moisture == "Wet" ~ 3),
                                                  DomGrazer = case_when(DomGrazer == "None" ~ 0,
                                                                        DomGrazer == "Insects" ~ 1,
                                                                        DomGrazer == "Small" ~ 2,
                                                                        DomGrazer == "Birds" ~ 3,
                                                                        DomGrazer == "Mixed" ~ 4,
                                                                        DomGrazer == "Large" ~ 5),
                                                  GrazerIntensity = case_when(GrazerIntensity == "Low" ~ 1,
                                                                              GrazerIntensity == "Medium" ~ 2,
                                                                              GrazerIntensity == "High" ~ 3),
                                                  permafrost = case_when(permafrost == "none" ~ 0,
                                                                         permafrost == "sporadic" ~ 1,
                                                                         permafrost == "continuous" ~ 2))

# rename variables for better fit in figure
col.data.local.num2 <- col.data.local.num %>% rename(Elevation = ElevationSubsite,
                                                       GrazerInt = GrazerIntensity,
                                                       Permafrost = permafrost,
                                                       PlotSize = SurveyedArea,
                                                       InitialBor = StartBorSpps)

# use only pairwise-complete observations
pair.cor.col.mod3 = cor(col.data.local.num2, use="pairwise.complete.obs", method = "spearman")

# FIGURE S2e
# visualize correlogram (1000 save)
(pairwise.col.mod3.num <- corrplot(pair.cor.col.mod3, type = 'lower', method = "color", 
                                    diag = FALSE, tl.col = 'black', tl.srt = 45,
                                   addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.5,
                                   cl.cex = 1))
# nothing majorly correlated


# the model doesn't run with Disko in because it's the only site with one value of 
# discontinuous permafrost and mixed moisture - removing this for convergence
col.data.inc3 <- col.data.inc2 %>% filter(SiteSubsite %notin% "DISKO:DISTURBANCE")

col_blab_local_mod <- brm(bf(BIC_B_LAB_corrected ~ Moisture + GrazerIntensity + 
                                DomGrazer + ElevationSubsite + permafrost + 
                                SurveyedArea + StartBorSpps + (1|SiteSubsite)),
                           data = col.data.inc3, family = 'beta', init = 0,
                           iter = 2000, chains = 4, warmup = 400, 
                           file = "models/dec2024_col_blab_local_mod")

print(summary(col_blab_local_mod), digits = 5) # ns
conditional_effects(col_blab_local_mod)


# elevation predictions 
col_blab_elev_pred <- ggpredict(col_blab_local_mod,
                                terms = "ElevationSubsite [sample = 45]")
colnames(col_blab_elev_pred) = c('ElevationSubsite', 'fit', 'lwr', 'upr', 'dunno')


# plot size predictions 
col_blab_ps_pred <- ggpredict(col_blab_local_mod,
                                terms = "SurveyedArea [sample = 10]")
colnames(col_blab_ps_pred) = c('SurveyedArea', 'fit', 'lwr', 'upr', 'dunno')


#### 03) Local model (all) ####
col.data.full3 <- col.data.full2 %>% filter(SiteSubsite %notin% "DISKO:DISTURBANCE")

col.data.local.counts <- col.data.full3 %>% left_join(., bic.counts, by = "SiteSubsitePlot")


col_blab_local_all_mod_bino <- brm(BorealCount | trials(TotalColonisers) ~ 
                                   Moisture + GrazerIntensity + 
                                   DomGrazer + ElevationSubsite + permafrost + 
                                   SurveyedArea + StartBorSpps + (1|SiteSubsite), 
                              data = col.data.local.counts, family = binomial("logit"),
                              init = 0, iter = 2000, chains = 4, warmup = 400, 
                              file = "models/dec2024_col_blab_local_all_mod_bino")

print(summary(col_blab_local_all_mod_bino), digits = 4) # plot size positive, elevation positive
conditional_effects(col_blab_local_all_mod_bino) #ns


## ELEVATION

# elevation predictions (full range)
col_blab_elev_fullrange_bino_pred <- ggpredict(col_blab_local_all_mod_bino, terms = "ElevationSubsite [sample = 45]")
colnames(col_blab_elev_fullrange_bino_pred) = c('ElevationSubsite', 'fit', 'lwr', 'upr', 'dunno')

colo.data.zero.local <- col.data.full3 %>% filter(BIC_B_LAB == 0)

# FIGURE 2g
(col_blab_elev_plot <- ggplot() +
    geom_point(data = col.data.inc3, aes(x = ElevationSubsite, y = BIC_B_LAB),
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") +
    geom_point(data = colo.data.zero.local, aes(x = ElevationSubsite, y = BIC_B_LAB),
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") +
    geom_line(data = col_blab_elev_pred, aes(x = ElevationSubsite, y = fit), colour = "#2c6611", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = col_blab_elev_pred, aes(x = ElevationSubsite, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_elev_fullrange_bino_pred, aes(x = ElevationSubsite, y = fit), colour = "black", linetype = "solid", linewidth = 2) +
    geom_ribbon(data = col_blab_elev_fullrange_bino_pred, aes(x = ElevationSubsite, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab("Elevation (m a.s.l.)") + ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.background = element_blank()))

ggsave(col_blab_elev_plot, filename = "figures/revised_fig2_bic_elev_plot.png",
       width = 20, height = 15, units = "cm")



## PLOT SIZE

# plot size predictions (full range)
col_blab_ps_fullrange_bino_pred <- ggpredict(col_blab_local_all_mod_bino, terms = "SurveyedArea [sample = 10]")
colnames(col_blab_ps_fullrange_bino_pred) = c('SurveyedArea', 'fit', 'lwr', 'upr', 'dunno')

# plot - BIC ~ plot size
(col_blab_ps_plot <- ggplot() +
    geom_point(data = col.data.inc3, aes(x = SurveyedArea, y = BIC_B_LAB),
               alpha = 0.15, size = 10, pch = 20, colour = "#2c6611") +
    geom_point(data = colo.data.zero.local, aes(x = SurveyedArea, y = BIC_B_LAB),
               alpha = 0.2, size = 10, pch = 20, colour = "#a19f9f") +
    geom_line(data = col_blab_ps_pred, aes(x = SurveyedArea, y = fit), colour = "#2c6611", linewidth = 2, linetype = "dashed") +
    geom_ribbon(data = col_blab_ps_pred, aes(x = SurveyedArea, ymin = lwr, ymax = upr), fill = "#2c6611", alpha = 0.3) +
    geom_line(data = col_blab_ps_fullrange_bino_pred, aes(x = SurveyedArea, y = fit), colour = "black", linetype = "solid", linewidth = 2) +
    geom_ribbon(data = col_blab_ps_fullrange_bino_pred, aes(x = SurveyedArea, ymin = lwr, ymax = upr), fill = "black", alpha = 0.3) +
    xlab(expression(bold(paste("Plot size (m"^bold("2"), ")")))) + 
                           ylab("Borealization (BCI)") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.background = element_blank()))


ggsave(col_blab_ps_plot, filename = "figures/revised_fig2_bic_plotsize.png",
       width = 20, height = 15, units = "cm")







## FIGURE S4 ----


## ABUNDANCE WITH ALL (MAP)

bia.blab.site.all <- bia.labb2 %>% 
  left_join(., site.only, by = "SiteSubsitePlot") %>% 
  left_join(., site.long, by = "SITE") %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))

# map of average cover change for each site
bia.blab.site.mean.all <- bia.blab.site.all %>% group_by(SITE) %>% 
  summarise(MeanCoverChange = mean(CoverChange))

length(unique(bia.blab.site.mean.all$SITE)) #32

# plot only B_LAB together
(bic.abun.blab.plot.all <- ggplot() +
    geom_point(data = bia.blab.site.all, aes(x = reorder(SITE, LONG), y = CoverChange),
               alpha = 0.3, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_point(data = bia.blab.site.mean.all, aes(x = SITE, y = MeanCoverChange), 
               size = 5, pch = 21, stroke = 1.5, colour = "black") +
    xlab("Study area") + ylab("Borealization (BAI)\n") +
    theme(axis.text.x  = element_text(angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          axis.text.y  = element_text(size = 15, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 20),
          axis.title.y = element_text(face="bold", size = 20),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          legend.background = element_blank()))

ggsave(bic.abun.blab.plot.all, filename = "figures/dec2024_points_blab_bia_new_all.png", 
       width = 30, height = 20, units = "cm")




# map

# join back in
bia.blab.site.mean.coord.all <- bia.blab.site.mean.all %>% left_join(site.coords, by = "SITE")

# Convert coords to utm
site.coords.utm.abun.all <- transform_coord(x = bia.blab.site.mean.coord.all, lon = "LONG", lat = "LAT", 
                                   new.names = c("lon.utm", "lat.utm"), 
                                   verbose = FALSE, bind = TRUE)

# Plot map at site level - unique species per site (with labels)
(site.bia.blab.map.all <- basemap(limits = 50, land.col = "#c9c7c1", 
                                  grid.size = 1, grid.col = "#8f8f8f") + 
    geom_point(data = site.coords.utm.abun.all, 
               aes(x = lon.utm, y = lat.utm, colour = MeanCoverChange), 
               size = 15, alpha = 0.7) + 
    #scale_colour_viridis(option = "D") +
    scale_colour_gradient(low = "#f2e81d", high = "#2c6611") +
    labs(colour='Mean BAI') +
    geom_label_repel(data = site.coords.utm.abun.all, aes(lon.utm, lat.utm, label = SITE), 
                     color = "black", box.padding = 2, 
                     segment.color = "black", segment.size = 1.2, 
                     fill = "white", label.size = 0.7,  
                     size = 9, max.iter = 150000, max.overlaps = 50) +
    theme(legend.title = element_text(size = 20), 
          legend.text=element_text(size = 20), 
          legend.position="right", 
          legend.margin=margin(1,1,1,1),
          legend.box.margin=margin(10,10,10,10), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")))

ggsave(site.bia.blab.map.all, filename = "figures/dec2024_map_blab_bia_new_all.png", 
       width = 30, height = 30, units = "cm")



## COLONISATIONS WITH ALL (MAP) 

## figure per site
points.bic.all <- master.class4.com %>% 
  select(SiteSubsitePlot, BIC_B_LAB) 

points.bic.site.all <- left_join(points.bic.all, site.only, by = "SiteSubsitePlot") %>% 
  left_join(., site.long, by = "SITE") %>%
  mutate(SITE = case_when(SITE == "BARROW" ~ "UTQIAĠVIK", TRUE ~ SITE))

length(unique(points.bic.site.all$SITE)) #32

# map of average BIC for each site
bic.blab.site.mean.all <- points.bic.site.all %>% group_by(SITE) %>% 
  summarise(MeanBIC = mean(BIC_B_LAB)) %>%
  left_join(., site.coords, by = "SITE")


# plot per site and plots
(bic.points.all <- ggplot() +
    geom_point(data = points.bic.site.all, aes(x = reorder(SITE, LONG), y = BIC_B_LAB), 
               alpha = 0.3, size = 10, pch = 20, colour = "#a19f9f") + 
    geom_point(data = bic.blab.site.mean.all, aes(x = SITE, y = MeanBIC), 
               size = 5, pch = 21, stroke = 1.5, colour = "black") +
    #scale_fill_viridis(option = "D") + 
    xlab("Study area") + ylab("Borealization (BCI)\n") +
    theme(axis.text.x  = element_text(angle = 52, 
                                      vjust = 1, hjust = 1,
                                      size = 15, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 20),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.title.y = element_text(face="bold", size = 20),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          #legend.position = c(0.95, 0.85), legend.key = element_blank(),
          legend.background = element_blank()))

ggplot2::ggsave(bic.points.all, 
                filename = "figures/dec2024_points_blab_bic_new_all.png", 
                width = 30, height = 20, units = "cm")


# Convert coords to utm
site.coords.utm2.all <- transform_coord(x = bic.blab.site.mean.all, lon = "LONG", lat = "LAT", 
                                    new.names = c("lon.utm", "lat.utm"), 
                                    verbose = FALSE, bind = TRUE)

# Plot map at site level - unique species per site (with labels)
(site.bic.blab.map.all <- basemap(limits = 50, land.col = "#c9c7c1", 
                                  grid.size = 1, grid.col = "#8f8f8f") + 
    geom_point(data = site.coords.utm2.all, 
               aes(x = lon.utm, y = lat.utm, colour = MeanBIC), 
               size = 15, alpha = 0.7) + 
    scale_colour_gradient(low = "#f2e81d", high = "#2c6611") +
    labs(colour='Mean BCI') +
    geom_label_repel(data = site.coords.utm2.all, aes(lon.utm, lat.utm, label = SITE), 
                     color = "black", box.padding = 2, 
                     segment.color = "black", segment.size = 1.2, 
                     fill = "white", label.size = 0.7,  
                     size = 9, max.iter = 150000, max.overlaps = 50) +
    theme(legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20), 
          legend.position = "right", 
          legend.margin = margin(1,1,1,1),
          legend.box.margin = margin(10,10,10,10), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")))

ggsave(site.bic.blab.map.all, filename = "figures/dec2024_map_blab_bic_new_all.png", 
       width = 30, height = 30, units = "cm")

