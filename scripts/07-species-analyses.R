## Arctic plant borealization
## Mariana Garcia Criado
## Script 7. Species-level analyses
## May 2024

## LIBRARIES ----
library(tidyverse)
library(broom)
library(brms)
library(corrplot)
library(ggeffects)
library(ggpubr)


## DATA ----

# only species and AB class
species.class0 <- read.csv("data/species_ab_class_jan2025.csv") 

# ITEX data
load("data/boreal_master_perm.RData")

# join class to main database and removes morphospecies
master.class <- boreal_master_perm2 %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  tidyr::drop_na(Class)

unique(master.class$SPECIES_CLEAN) #287 species



## COLONISERS /  START-END ----

# categorise species trajectories
master.class2 <- master.class %>% group_by(SiteSubsitePlot) %>% 
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

# keep colonisers only
list.sp.times.col <- master.class2 %>% 
  filter(trend == "Coloniser") %>%
  group_by(SPECIES_CLEAN) %>% 
  summarise(TimesColonised = n()) %>% 
  ungroup()
  

# table of colonisers:
list.sp.times.col.table <- list.sp.times.col %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  select(SPECIES_CLEAN, TimesColonised, Class) %>%
  arrange(-TimesColonised) %>%
  rename(Species = SPECIES_CLEAN)
write.csv(list.sp.times.col.table, "data/dec2024_list_sp_time_col_table.csv")


## ABUNDANCES / START-END ----
list.sp.cover <- master.class %>% 
  select(SiteSubsitePlotYear, SiteSubsitePlot, SPECIES_CLEAN, YEAR, Class, RelCover, Duration) %>%
  group_by(SiteSubsitePlot) %>% 
  mutate(MinYear = min(YEAR)) %>% mutate(MaxYear = max(YEAR)) %>%
  filter(MinYear != MaxYear) %>%
  filter(YEAR == MinYear | YEAR == MaxYear) %>%
  ungroup() %>%
  mutate(Timepoint = case_when(YEAR == MinYear ~ "Start",
                               YEAR == MaxYear ~ "End")) %>%
  select(SiteSubsitePlot, SPECIES_CLEAN, RelCover, Timepoint, Duration) %>%
  mutate(row = row_number()) %>%
  group_by(SiteSubsitePlot, SPECIES_CLEAN) %>%
  pivot_wider(names_from = Timepoint, values_from = RelCover) %>%
  summarise_each(funs(first(.[!is.na(.)]))) %>%
  mutate_at(vars(Start:End), ~replace_na(., 0)) %>%
  ungroup() %>%
  select(-row) %>%
  mutate(CoverChangePerYear = (End - Start)/Duration)

# calculate n per species and confidence intervals
list.sp.cover.cin <- list.sp.cover %>% 
  group_by(SPECIES_CLEAN) %>%
  mutate(n = length(SPECIES_CLEAN)) %>%
  mutate(ConfidenceIntervals = (sd(CoverChangePerYear)/sqrt(n))*1.96) %>%
  ungroup()

# join in categorical and numerical traits
all.traits <- read.csv("data/species_traits_dataset_full.csv") 

list.sp.cover.all <- list.sp.cover.cin %>% 
  left_join(., all.traits, by = "SPECIES_CLEAN")

# calculate overall cover change  
list.sp.cover.mean <- list.sp.cover.all %>% 
  group_by(SPECIES_CLEAN) %>%
  mutate(MeanCoverChange = mean(CoverChangePerYear)) %>%
  ungroup() %>%
  distinct(SPECIES_CLEAN, .keep_all = T) %>%
  select(-SiteSubsitePlot) %>%
  left_join(., list.sp.times.col, by = "SPECIES_CLEAN") %>%
  mutate_at(vars(TimesColonised), ~replace_na(., 0)) %>%
  mutate(ClassShort = case_when(ClassNew == "Boreal specialist" ~"B",
                                ClassNew == "Boreal-tundra boundary" ~ "BTB",
                                ClassNew == "Arctic specialist" ~ "A",
                                ClassNew == "Ubiquitous" ~ "U"))

# NOTE: we have 272 species, and not 287 as in the original dataset.  
# Because there are species that only appear in between timepoints, so not present at start and/or end
# thus they are being filtered out when using start vs end

# save table with both cols and cover change
list.sp.cover.times.table <- list.sp.cover.mean %>% 
  select(SPECIES_CLEAN, TimesColonised, MeanCoverChange, ConfidenceIntervals, n, 
         ClassShort, FunctionalGroup, PlantHeight) %>%
  arrange(-TimesColonised)

# TABLE S4
write.csv(list.sp.cover.times.table, "data/dec2024_list_sp_time_col_cov_table.csv")


## Are the species that have colonised the most also the ones with greatest increases?
times.cover.mod <- brm(TimesColonised ~ MeanCoverChange, data = list.sp.cover.mean,
                   family = 'negbinomial',
                   prior = prior(gamma(0.01, 0.01), class = shape),
                   iter = 2000, chains = 4, warmup = 400, 
                   file = "models/dec2024_times_cov_species_mod")
summary(times.cover.mod) # positive (marginally significant)


# check distributions
hist(list.sp.cover.mean$TimesColonised, breaks = 30)
hist(list.sp.cover.mean$MeanCoverChange, breaks = 30)

mean(list.sp.cover.mean$TimesColonised) # 7.48
var(list.sp.cover.mean$TimesColonised) # 134.8

# fix a few things to enable convergence
list.sp.cover.mean2 <- list.sp.cover.mean %>% 
  mutate(Woodiness = case_when(Woodiness == 1 ~ "Woody", 
                               TRUE ~ "Not woody"),
         N_fixer = case_when(N_fixer == 1 ~ "Fixer",
                             TRUE ~ "Not fixer"))

write.csv(list.sp.cover.mean2, "data/dec2024_species_traits_db_for_mods.csv")



## TRAIT CORRELATIONS ----

# we should make two: one for times colonised and another for abun change, as they are 2 different data subsets
b.lab.inc <- list.sp.cover.mean2 %>% 
  filter(ClassNew %in% c("Boreal-tundra boundary", "Boreal specialist"), 
         MeanCoverChange > 0)

# convert variables with 2 categories for numerical to make a correlogram
cord.db <- b.lab.inc %>% 
  mutate(BerryCat = case_when(BerryCat == "Berry" ~ 1, TRUE ~ 0),
         Woodiness = case_when(Woodiness == "Not woody" ~ 0, TRUE ~ 1),
         N_fixer = case_when(N_fixer == "Not fixer" ~ 0, TRUE ~ 1)) %>%
  select(Woodiness, BerryCat, N_fixer, PlantHeight, SLA, SeedMass, LeafN, LeafCN)

str(cord.db)

cord.db.log <- cord.db %>% mutate(LogPlantHeight = log(PlantHeight),
                                  LogSLA = log(SLA),
                                  LogSeed = log(SeedMass),
                                  LogLeafN = log(LeafN),
                                  LogLeafCN = log(LeafCN)) %>%
  select(Woodiness, BerryCat, N_fixer, LogPlantHeight, LogSLA, LogSeed, LogLeafN, LogLeafCN)

# create matrix
pair.cor.traits = cor(cord.db.log, use = "pairwise.complete.obs", method = "spearman")

# FIGURE S3b
# visualize correlogram (1000 width)
(trait.plot.circle <- corrplot(pair.cor.traits, type = 'lower', method = "color",
                               diag = FALSE, tl.col = 'black', tl.srt = 45, na.label = "NA",
                               addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.5,
                               cl.cex = 1))

# leaf N correlated with leaf CN and seed mass, we remove leafN

## Times colonised pairwise correlations
onecol.blab <- list.sp.cover.mean2 %>%
  filter(ClassNew %in% c("Boreal-tundra boundary", "Boreal specialist"),
         TimesColonised > 0)

# convert variables with 2 categories for numerical to make a correlogram
cord.db.col <- onecol.blab %>% 
  mutate(BerryCat = case_when(BerryCat == "Berry" ~ 1, TRUE ~ 0),
         Woodiness = case_when(Woodiness == "Not woody" ~ 0, TRUE ~ 1),
         N_fixer = case_when(N_fixer == "Not fixer" ~ 0, TRUE ~ 1)) %>%
  select(Woodiness, BerryCat, N_fixer, PlantHeight, SLA, SeedMass, LeafN, LeafCN)

str(cord.db.col)

# Traits are included in the models in a log-scale, so let's check correlations with log-data
cord.db.col.log <- cord.db.col %>% mutate(LogPlantHeight = log(PlantHeight),
                                  LogSLA = log(SLA),
                                  LogSeed = log(SeedMass),
                                  LogLeafN = log(LeafN),
                                  LogLeafCN = log(LeafCN)) %>%
  select(Woodiness, BerryCat, N_fixer, LogPlantHeight, LogSLA, LogSeed, LogLeafN, LogLeafCN)

# create matrix
pair.cor.traits.col = cor(cord.db.col.log, use = "pairwise.complete.obs", method = "spearman")

# FIGURE S3a
(trait.plot.col.circle <- corrplot(pair.cor.traits.col, type = 'lower', method = "color",
                               diag = FALSE, tl.col = 'black', tl.srt = 45, na.label = "NA",
                               addCoef.col = 'black', number.cex = 1.8, tl.cex = 1.5,
                               cl.cex = 1))



# the only correlation is between Leaf N and Leaf CN which is to be expected
# we retain Leaf N which has a greater number of traits
# leaf CN also correlated with seed mass and SLA.
# removing woodiness because related to FG (shrub)
# family removed for now as some species only have one species, might be problematic
# removing berry because the majority are non-berry (graminoids, forbs)
# removing N_fixer because the group of N_fixer species is super small and it's brekaing the models
# all N_fixers are forbs (and fabaceae) so this will be implicit in the forb func group
# deciduousness not included as very small sample size (all other species are non dec)

family.sum <- list.sp.cover.mean2 %>% group_by(family) %>% summarise(n = n())




## ABUNDANCE MODELS ----

# 1) ABUNDANCE: ALL SPECIES, INCREASES ONLY (n = 129)
abun.inc <- list.sp.cover.mean2 %>% 
  filter(MeanCoverChange > 0)

# sample sizes per group
n.sizes.abun <- abun.inc %>% group_by(ClassNew) %>% summarise(n())


# Which group of species has increased the most in abundance over time?
abun.class <- brm(MeanCoverChange ~ ClassNew, 
                  data = abun.inc, family = 'gaussian',
                  iter = 2000, chains = 4, warmup = 400, 
                  file = "models/dec2024_species_abun_class_inc_mod")
summary(abun.class) # ns
conditional_effects(abun.class)

# Predictions
abun.class.pred <- as.data.frame(ggpredict(abun.class, terms = "ClassNew")) %>% select(-group)
colnames(abun.class.pred) = c('Class', 'fit', 'lwr', 'upr')

str(abun.class.pred)

# organise classes per distribution
abun.class.pred2 <- abun.class.pred %>% 
  mutate(Class = case_when(Class == "Boreal-tundra boundary" ~ "Boreal-Tundra", 
                           Class == "Boreal specialist" ~ "Boreal",
                           Class == "Arctic specialist" ~ "Arctic",
                           TRUE ~ Class)) %>%
  mutate(Order = case_when(Class == "Boreal" ~ 1,
                           Class == "Boreal-tundra" ~ 2,
                           Class == "Arctic" ~ 3,
                           Class == "Ubiquitous" ~ 4))


# FIGURE 3b
abun.inc2 <- abun.inc %>% 
  mutate(ClassShort = case_when(ClassNew == "Boreal-tundra boundary" ~ "Boreal-Tundra", 
                                ClassNew == "Boreal specialist" ~ "Boreal",
                                ClassNew == "Arctic specialist" ~ "Arctic",
                                TRUE ~ ClassNew)) %>%
  mutate(Order = case_when(ClassShort == "Boreal" ~ 1,
                           ClassShort == "Boreal-Tundra" ~ 2,
                           ClassShort == "Arctic" ~ 3,
                           ClassShort == "Ubiquitous" ~ 4))

# violin together with model predictions
(abun.class.violin <- ggplot() + 
    geom_violin(data = abun.inc2, aes(x = reorder(ClassShort, Order), y = MeanCoverChange), 
                size = 0, colour = "#f5eac9", fill = "#f5eac9") +
    geom_errorbar(data = abun.class.pred2, aes(x = reorder(Class, Order), ymin = lwr, ymax = upr), 
                  size = 1, width = 0.25) + 
    geom_point(data = abun.class.pred2, aes(x = reorder(Class, Order), y = fit), 
               size = 10, pch = 21, fill = "#a19f9f", colour = "black") + 
    xlab("\nClass") + ylab("Mean abundance increase (%)\n") + 
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.text.y  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"), 
          legend.position = "none", legend.key = element_blank(),
          legend.background = element_blank()))





# 2) MULTIVARIATE, ABUNDANCE: BOREAL + LAB SPECIES INCREASES ONLY
b.lab.inc <- list.sp.cover.mean2 %>% 
  filter(ClassNew %in% c("Boreal-tundra boundary", "Boreal specialist"), 
         MeanCoverChange > 0)
  
# Which traits are associated with B+LAB species increasing in abundance over time? (n = 24)
species.abun.full.blab.inc.mod <- brm(MeanCoverChange ~ FunctionalGroup + 
                                        log(PlantHeight) + log(SLA) + log(SeedMass) + log(LeafN), 
                                      data = b.lab.inc, family = 'gaussian',
                                      iter = 2000, chains = 4, warmup = 400, 
                                      file = "models/dec2024_species_abun_full_blab_inc_mod")

summary(species.abun.full.blab.inc.mod) # all ns
conditional_effects(species.abun.full.blab.inc.mod)



# 3) UNIVARIATE, ABUNDANCE: B+LAB INCREASE ONLY

# functional group (n = 82)
species.abun.blab.inc.fg.mod <- brm(MeanCoverChange ~ FunctionalGroup,
                                data = b.lab.inc, family = 'gaussian',
                                iter = 2000, chains = 4, warmup = 400,
                                file = "models/dec2024_species_abun_blab_inc_fg_mod")
summary(species.abun.blab.inc.fg.mod)
conditional_effects(species.abun.blab.inc.fg.mod) #ns


# height (n = 58)
species.abun.blab.inc.hei.mod <- brm(MeanCoverChange ~ log(PlantHeight),
                                   data = b.lab.inc, family = 'gaussian',
                                   iter = 2000, chains = 4, warmup = 400,
                                   file = "models/dec2024_species_abun_blab_inc_hei_mod")
summary(species.abun.blab.inc.hei.mod) #ns

# SLA (n = 51)
species.abun.blab.inc.sla.mod <- brm(MeanCoverChange ~ log(SLA),
                                 data = b.lab.inc, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400,
                                 file = "models/dec2024_species_abun_blab_inc_sla_mod")
summary(species.abun.blab.inc.sla.mod) #ns

# Seed (n = 27)
species.abun.inc.blab.seed.mod <- brm(MeanCoverChange ~ log(SeedMass),
                                 data = b.lab.inc, family = 'gaussian',
                                 iter = 2000, chains = 4, warmup = 400,
                                 file = "models/dec2024_species_abun_blab_inc_seed_mod")
summary(species.abun.inc.blab.seed.mod) #ns

# leaf N (n = 41)
species.abun.inc.blab.leafn.mod <- brm(MeanCoverChange ~ log(LeafN),
                                  data = b.lab.inc, family = 'gaussian',
                                  iter = 2000, chains = 4, warmup = 400,
                                  file = "models/dec2024_species_abun_blab_inc_leafn_mod")
summary(species.abun.inc.blab.leafn.mod) #ns

# All univariate models yield the same significance as the multivariate model - great!





## COLONIZATIONS MODELS ----


# 1) COLONISATIONS: ALL SPECIES TIMESCOLONISED > 0
onecol <- list.sp.cover.mean2 %>% filter(TimesColonised > 0)
hist(onecol$TimesColonised)

# Which group of species has colonised the most over time 
# (from those species that have colonised at least once? (n = 220)
col.class.mod2 <- brm(TimesColonised ~ ClassNew, data = onecol, 
                     family = 'negbinomial',
                     prior = prior(gamma(0.01, 0.01), class = shape),
                     iter = 2000, chains = 4, warmup = 400, 
                     file = "models/dec2024_species_col_class_inc_mod2")
summary(col.class.mod2) 
conditional_effects(col.class.mod2) # B have colonised fewer times than BT and Cosmopolitan, the others overlap



# sample sizes per group
n.sizes <- onecol %>% group_by(ClassNew) %>% summarise(n())



# Predictions
col.class.pred <- as.data.frame(ggpredict(col.class.mod2, terms = "ClassNew")) %>% select(-group)
colnames(col.class.pred) = c('Class', 'fit', 'lwr', 'upr')

str(col.class.pred)



# organise classes per distribution
col.class.pred2 <- col.class.pred %>% 
  mutate(Class = case_when(Class == "Boreal-tundra boundary" ~ "Boreal-Tundra", 
                           Class == "Boreal specialist" ~ "Boreal",
                           Class == "Arctic specialist" ~ "Arctic",
                           TRUE ~ Class)) %>%
  mutate(Order = case_when(Class == "Boreal" ~ 1,
                           Class == "Boreal-tundra" ~ 2,
                           Class == "Arctic" ~ 3,
                           Class == "Ubiquitous" ~ 4))


# FIGURE 3a
onecol2 <- onecol %>% 
  mutate(ClassShort = case_when(ClassNew == "Boreal-tundra boundary" ~ "Boreal-Tundra", 
                                ClassNew == "Boreal specialist" ~ "Boreal",
                                ClassNew == "Arctic specialist" ~ "Arctic",
                                TRUE ~ ClassNew)) %>%
  mutate(Order = case_when(ClassShort == "Boreal" ~ 1,
                           ClassShort == "Boreal-Tundra" ~ 2,
                           ClassShort == "Arctic" ~ 3,
                           ClassShort == "Ubiquitous" ~ 4))

# violin together with model predictions
(col.class.violin <- ggplot() + 
    geom_violin(data = onecol2, aes(x = reorder(ClassShort, Order), y = TimesColonised), 
                size = 0, colour = "#f5eac9", fill = "#f5eac9") +
    geom_errorbar(data = col.class.pred2, aes(x = reorder(Class, Order), ymin = lwr, ymax = upr), 
                  size = 1, width = 0.25) + 
    geom_point(data = col.class.pred2, aes(x = reorder(Class, Order), y = fit), 
               size = 10, pch = 21, fill = "#a19f9f", colour = "black") + 
    xlab("\nClass") + ylab("Plots colonised\n") + 
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.text.y  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"), 
          legend.position = "none", legend.key = element_blank(),
          legend.background = element_blank()))



# FIGURE 3
(class.violin.panel <- ggarrange(col.class.violin, abun.class.violin, 
                                 labels = c("a)", "b)"), align = "hv",
                                 nrow = 1, ncol = 2, font.label = list(size = 28, face = "plain")))

ggsave(class.violin.panel, filename = "figures/june2025_class_panel.png", 
       width = 55, height = 22, units = "cm")



# 2) COLONISATIONS, MULTIVARIATE: B+LAB SPECIES (TIMES COLONISED > 0)
onecol.blab <- onecol %>% filter(ClassNew %in% c("Boreal specialist", "Boreal-tundra boundary"))

onecol.blab.log <- onecol.blab %>% mutate(LogPlantHeight = log(PlantHeight),
                                          LogSLA = log(SLA),
                                          LogSeed = log(SeedMass),
                                          LogLeafN = log(LeafN))


# Same model with log() as function
 species.col.full.blab.mod2 <- brm(TimesColonised ~ FunctionalGroup + 
                                     log(PlantHeight) + log(SLA) + log(SeedMass) + log(LeafN), 
                                   data = onecol.blab, family = 'negbinomial',
                                   prior = prior(gamma(0.01, 0.01), class = shape), 
                                   iter = 2000, chains = 4, warmup = 400, 
                                   file = "models/dec2024_species_col_full_blab_mod2")
  
  summary(species.col.full.blab.mod2) # height negative
  conditional_effects(species.col.full.blab.mod2) # forb have colonised fewer times than shrubs, small overlap with grams

  # sample size per category of func group
  n.fg <- onecol.blab %>% group_by(FunctionalGroup) %>% summarise(n = n())


# func group predictions
sp.fg.pred <- as.data.frame(ggpredict(species.col.full.blab.mod2, terms = "FunctionalGroup")) %>% select(-group)
colnames(sp.fg.pred) = c('FunctionalGroup', 'fit', 'lwr', 'upr')


# FIGURE 4b
(sp.fg.violin <- ggplot() + 
    geom_violin(data = onecol.blab, aes(x = FunctionalGroup, y = TimesColonised, 
                                        fill = FunctionalGroup), colour = NA, size = 0) +
    geom_errorbar(data = sp.fg.pred, aes(x = FunctionalGroup, ymin = lwr, ymax = upr, colour = FunctionalGroup), 
                  size = 1, width = 0.25) + 
    geom_point(data = sp.fg.pred, aes(x = FunctionalGroup, y = fit, colour = FunctionalGroup), 
               size = 8) + 
    scale_colour_manual(values = c("#95558b", "#dc9537", "#5fbf33"), 
                        breaks = c("Forb", "Graminoid", "Shrub"), 
                        name = "Functional Group") + 
    scale_fill_manual(values = c("#95558b20", "#dc953720", "#5fbf3320"), 
                      breaks = c("Forb", "Graminoid", "Shrub"), 
                      name = "Functional Group") + 
    xlab("\nFunctional group") + ylab("Plots colonised\n") + 
    theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.text.y  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(size = 24, face = "bold"), 
          axis.title.y = element_text(size = 24, face = "bold"), 
          legend.position = "none", legend.key = element_blank(),
          legend.background = element_blank()))

ggsave(sp.fg.violin, filename = "figures/june2025_sp_fg_plot.png", 
       height = 15, width = 20, units = "cm")


# 3) COLONISATIONS, UNIVARIATE MODELS: B+LAB ONLY

# Univariate, FG
species.col.blab.fg.mod2 <- brm(TimesColonised ~ FunctionalGroup, 
                               data = onecol.blab, family = 'negbinomial',
                               prior = prior(gamma(0.01, 0.01), class = shape), 
                               iter = 2000, chains = 4, warmup = 400, 
                               file = "models/dec2024_species_col_blab_fg_mod2")
summary(species.col.blab.fg.mod2)
conditional_effects(species.col.blab.fg.mod2) # shrub colonised more than forbs


# Univariate, height (without 0s)
species.col.blab.hei.mod5 <- brm(TimesColonised ~ log(PlantHeight), 
                                 data = onecol.blab.log, family = 'negbinomial',
                                 prior = prior(gamma(0.01, 0.01), class = shape), 
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/dec2024_species_col_blab_height_mod5")

summary(species.col.blab.hei.mod5) # height negative

# predictions
sp_hei_uni_pred5 <- ggpredict(species.col.blab.hei.mod5, 
                             terms = "PlantHeight [all]", back.transform = TRUE)
colnames(sp_hei_uni_pred5) = c('PlantHeight', 'fit', 'lwr', 'upr', 'dunno')

# FIGURE 4a
(sp_hei_plot_bt <- ggplot() +
    geom_point(data = onecol.blab.log, aes(x = PlantHeight, y = TimesColonised, fill = FunctionalGroup), 
               alpha = 0.5, size = 6, pch = 21, colour = "black", stroke = 1) + 
    geom_line(data = sp_hei_uni_pred5, aes(x = PlantHeight, y = fit), colour = "black", linewidth = 1.2) +
    geom_ribbon(data = sp_hei_uni_pred5, aes(x = PlantHeight, ymin = lwr, ymax = upr), fill = "black", alpha = 0.2) +
    scale_fill_manual(values = c("#95558b", "#dc9537", "#5fbf33"), name = "Functional Group") + 
    xlab("\nPlant height (m)") + ylab("Plots colonised\n") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 25),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 25),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position = c(0.85, 0.85),
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20),
          legend.key=element_blank(),
          legend.background = element_blank()))

ggsave(sp_hei_plot_bt, filename = "figures/dec2024_sp_hei_uni_plot_backt.png", 
       width = 20, height = 15, units = "cm")



# Univariate, SLA (without 0s)
species.col.blab.sla.mod2 <- brm(TimesColonised ~ log(SLA), 
                                data = onecol.blab, family = 'negbinomial',
                                prior = prior(gamma(0.01, 0.01), class = shape), 
                                iter = 2000, chains = 4, warmup = 400, 
                                file = "models/dec2024_species_col_blab_sla_mod2")
summary(species.col.blab.sla.mod2) # SLA negative



# Univariate, seed (without 0s)
species.col.blab.seed.mod2 <- brm(TimesColonised ~ log(SeedMass), 
                                 data = onecol.blab, family = 'negbinomial',
                                 prior = prior(gamma(0.01, 0.01), class = shape), 
                                 iter = 2000, chains = 4, warmup = 400, 
                                 file = "models/dec2024_species_col_blab_seed_mod2")
summary(species.col.blab.seed.mod2) #ns



# Univariate, leaf N (without 0s)
species.col.blab.leafn.mod2 <- brm(TimesColonised ~ log(LeafN), 
                                  data = onecol.blab, family = 'negbinomial',
                                  prior = prior(gamma(0.01, 0.01), class = shape), 
                                  iter = 2000, chains = 4, warmup = 400, 
                                  file = "models/dec2024_species_col_blab_leafn_mod2")
summary(species.col.blab.leafn.mod2) # ns


## UNIVARIATE MODS WITHOUT 0: 
## GRAMS+SHRUBS COLONISED MORE, MULTIVARITE MODS: SHRUBS COLONISED MORE
## SLA UNIV NEGATIVE, SLA MULTIV NS


## UNIVARIATE MODS: no difference between with and without 0s

## MULTIVARIATE WITH 0S: SLA POSITIVE, MULTIVARIATE WITHOUT 0S NS, UNIVARIATE WITH AND WITHOUT 0S SLA NEGATIVE



## CHI-SQUARE TEST ----

# Chi-square test of independence
# To test if the proportions per FG and trend reflect that of the main database

## CHI: COLONISING SPECIES ##
col.chi <- onecol.blab %>% 
  select(SPECIES_CLEAN, FunctionalGroup) %>% 
  group_by(FunctionalGroup) %>% summarise(Col_Blab_Count = n())

# overall database (all categories, all trends)
all.sp.chi <- read.csv("data/dec2024_list_sp_time_col_cov_table.csv") %>%
  select(SPECIES_CLEAN, FunctionalGroup) %>%
  group_by(FunctionalGroup) %>% summarise(All_Count = n())

# Merge into one dataframe
compare.chi <- left_join(all.sp.chi, col.chi, by = "FunctionalGroup") %>% 
  column_to_rownames(var = "FunctionalGroup")

# Convert to contingency table
trends.matrix <- as.table(as.matrix(compare.chi))

chisq <- chisq.test(trends.matrix, correct=FALSE)
chisq # p = 0.26, the rows and column are not statistically associated

prop.test(trends.matrix)
# p > 0.05 so the proportions of the two groups are similar between the global database and the groups

chisq$observed
chisq$expected
round(chisq$residuals, 3)


## CHI: ABUNDANCE CHANGE ##
abun.chi <- b.lab.inc %>% 
  select(SPECIES_CLEAN, FunctionalGroup) %>% 
  group_by(FunctionalGroup) %>% summarise(Abun_Blab_Count = n())

# Merge into one dataframe
compare.abun.chi <- left_join(all.sp.chi, abun.chi, by = "FunctionalGroup") %>% 
  column_to_rownames(var = "FunctionalGroup")

# Convert to contingency table
trends.abun.matrix <- as.table(as.matrix(compare.abun.chi))

chisq.abun <- chisq.test(trends.abun.matrix, correct=FALSE)
chisq.abun # p = 0.148, the rows and column are not statistically associated

prop.test(trends.abun.matrix)
# p > 0.05 so the proportions of the two groups are similar between the global database and the groups

chisq$observed
chisq$expected
round(chisq$residuals, 3)





## COLONISERS PRESENT IN SUBSITES AT START OF MONITORING? ----
col.sp <- master.class2 %>% filter(trend == "Coloniser") %>% 
  mutate(Species_Plot = paste0(SPECIES_CLEAN, "_", SiteSubsite)) %>%
  select(SPECIES_CLEAN, SiteSubsitePlot, Species_Plot)
  
# list of species that were colonisers within their subsites
col.sp.vector <- col.sp$Species_Plot  
  
# were the plot colonisers present initially in the subsite?  
col.subs.db <- master.class2 %>% 
  mutate(Species_Plot = paste0(SPECIES_CLEAN, "_", SiteSubsite)) %>%
  filter(Species_Plot %in% col.sp.vector) %>%
  mutate(PresentAtSubsite = case_when(YEAR == MinYear ~ "Present at start",
                                      YEAR == MaxYear ~ "Colonisation point")) %>%
  mutate(SubsitePresenceNumber = case_when(PresentAtSubsite == "Present at start" ~ 1, 
                                           TRUE ~ 0)) %>%
  group_by(Species_Plot) %>%
  mutate(SubsitePresence = sum(SubsitePresenceNumber)) %>%
  ungroup() %>%
  distinct(Species_Plot, .keep_all = TRUE)
  
# total colonisers at all subsites (same species will be repeated) = 839
# number of colonisers that were present in the subsite previously = 540
col.present <- col.subs.db %>% filter(SubsitePresence > 0)
# 540/839 = 64.4%


