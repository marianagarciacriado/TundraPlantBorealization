## Arctic plant borealization
## Mariana Garcia Criado 
## Script 9. Net change vs BCI/Times colonised
## Jan 2025

# This script compares BCI/times colonised (as used in the manuscript) with net species change 
# (another way of calculating change, refer to paper Methods for rationale behind this)


## PACKAGES ----
library(tidyverse)
library(brms)
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

# save smaller reference datasets
duration.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, Duration)

site.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SITE)

subsite.only <- master.class.com %>% distinct(SiteSubsitePlot, .keep_all = T) %>% 
  select(SiteSubsitePlot, SiteSubsite)

plots.only <- master.class.com %>% distinct(SiteSubsitePlot)

site.long <- master.class.com %>% distinct(SITE, .keep_all = T) %>%
  select(SITE, LONG)




## BCI VS NET ----
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

# calculate the Borealization Colonization Index (as originally done)
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

# plots without colonisers (boreal or otherwise) are missing from this dataset - 
# need to fill them up with 0
# as they still count as zero borealization (even if they had no colonisations at all)
master.class4.com <- plots.only %>% 
  left_join(master.class3.com, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0) %>%
  left_join(., subsite.only, by = "SiteSubsitePlot")

## Calculate as net change

# extinct first
extinct <- master.class2.com %>% 
  filter(trend == "Extinct") %>%
  group_by(SiteSubsitePlot, Class) %>% 
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(SiteSubsitePlot) %>%
  mutate(TotalExtinct = sum(n)) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  rename(Boreal_Extinct = Boreal,
         Cosmopolitan_Extinct = Cosmopolitan,
         LAB_Extinct = Low_Arctic_Boreal,
         HLA_Extinct = High_Low_Arctic,
         Low_Arctic_Extinct = Low_Arctic)

# net change calculation
net.change <- master.class4.com %>%
  left_join(extinct, by = "SiteSubsitePlot") %>%
  replace(is.na(.), 0) %>%
  mutate(BorealNetChange = (Low_Arctic_Boreal + Boreal) - (LAB_Extinct + Boreal_Extinct)) %>%
  mutate(BorealNetChangeRelative = BorealNetChange / (TotalColonisers - TotalExtinct))


# Are shorter/longer timeseries picking up greater rates of borealization?
(xyz <- ggplot(net.change, aes(x = BIC_B_LAB, y = BorealNetChange)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))

# model
xyz.mod <- brm(BorealNetChange ~ BIC_B_LAB,
                      data = net.change, 
                      iter = 2000, chains = 4, warmup = 400, 
                      file = "models/net_vs_bic")
summary(xyz.mod) # positive

# model with subsite as random effect
xyz.mod.re <- brm(BorealNetChange ~ BIC_B_LAB + (1|SiteSubsite),
               data = net.change, 
               iter = 2000, chains = 4, warmup = 400, 
               file = "models/net_vs_bic_re")
summary(xyz.mod.re) # positive


# predictions
com.net.mod <- ggpredict(xyz.mod.re, terms = "BIC_B_LAB")
colnames(com.net.mod) = c('BIC_B_LAB', 'fit', 'lwr', 'upr', 'dunno')


# FIGURE S1a
(com_net_plot <- ggplot() +
    geom_point(data = net.change, aes(x = BIC_B_LAB, y = BorealNetChange), 
               alpha = 0.15, size = 10, pch = 20, colour = "#525453") + 
    geom_line(data = com.net.mod, aes(x = BIC_B_LAB, y = fit), colour = "#525453", linewidth = 2, linetype = "solid") +
    geom_ribbon(data = com.net.mod, aes(x = BIC_B_LAB, ymin = lwr, ymax = upr), fill = "#525453", alpha = 0.3) +
    xlab("\nBCI") + ylab("Boreal net change\n") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))



## TIMES COLONISED VS NET ----

# only species and AB class
species.class0 <- read.csv("data/species_ab_class_jan2025.csv") 

# ITEX data
load("data/boreal_master_perm.RData")

# join class to main database and removes morphospecies
master.class <- boreal_master_perm2 %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  tidyr::drop_na(Class)

unique(master.class$SPECIES_CLEAN) #287 species


# categorise species trajectories
sp.class2 <- master.class %>% group_by(SiteSubsitePlot) %>% 
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

length(unique(sp.class2$SPECIES_CLEAN)) #272

# keep colonisers only
list.sp.times.col <- sp.class2 %>% 
  filter(trend == "Coloniser") %>%
  group_by(SPECIES_CLEAN) %>% 
  summarise(TimesColonised = n()) %>% 
  ungroup()

# keep extincts only
list.sp.times.ext <- sp.class2 %>% 
  filter(trend == "Extinct") %>%
  group_by(SPECIES_CLEAN) %>% 
  summarise(TimesExtinct = n()) %>% 
  ungroup()

# calculate net change
net.change.spps <- list.sp.times.col %>%
  full_join(list.sp.times.ext, by = "SPECIES_CLEAN") %>%
  replace(is.na(.), 0) %>%
  mutate(SpeciesNetChange = TimesColonised - TimesExtinct)


# model
sp.mod <- brm(SpeciesNetChange ~ TimesColonised,
               data = net.change.spps, 
               iter = 2000, chains = 4, warmup = 400, 
               file = "models/sp_net_vs_times")
summary(sp.mod) # positive

# plot
(xyz.sp <- ggplot(net.change.spps, aes(x = TimesColonised, y = SpeciesNetChange)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))



# boreal species only
net.change.sp.class <- net.change.spps %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  filter(ClassNew %in% c("Boreal specialist", "Boreal-tundra boundary"))

# model
sp.mod.bor <- brm(SpeciesNetChange ~ TimesColonised,
              data = net.change.sp.class, 
              #family = 'zero_one_inflated_beta',
              #prior = prior(gamma(0.01, 0.01), class = phi), inits = 0,
              iter = 2000, chains = 4, warmup = 400, 
              file = "models/sp_net_vs_times_bor")
summary(sp.mod.bor) # positive


# plot
(xyz.sp.bor <- ggplot(net.change.sp.class, aes(x = TimesColonised, y = SpeciesNetChange)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))



## including zeroes for all species ##
sp.only <- sp.class2 %>% distinct(SPECIES_CLEAN)

sp.all <- sp.only %>% 
  left_join(., list.sp.times.col, by = "SPECIES_CLEAN") %>%
  left_join(., list.sp.times.ext, by = "SPECIES_CLEAN") %>%
  replace(is.na(.), 0) %>%
  mutate(SpeciesNetChange = TimesColonised - TimesExtinct)


# extinct vs coloniser species 
ext.col.final <- brm(TimesExtinct ~ TimesColonised,
                  data = sp.all, 
                  iter = 2000, chains = 4, warmup = 400, 
                  file = "models/ext_vs_col_all")
summary(ext.col.final) # positive
conditional_effects(ext.col.final)

# net change vs colonisations
sp.mod.all <- brm(SpeciesNetChange ~ TimesColonised,
              data = sp.all, 
              iter = 2000, chains = 4, warmup = 400, 
              file = "models/sp_net_vs_times_all")
summary(sp.mod.all) # positive

# plot
(xyz.sp.all <- ggplot(sp.all, aes(x = TimesColonised, y = SpeciesNetChange)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))



# boreal species only
net.change.sp.class.all <- sp.all %>% 
  left_join(., species.class0, by = "SPECIES_CLEAN") %>%
  filter(ClassNew %in% c("Boreal specialist", "Boreal-tundra boundary"))

# model
sp.mod.bor.all <- brm(SpeciesNetChange ~ TimesColonised,
                  data = net.change.sp.class.all, 
                  iter = 2000, chains = 4, warmup = 400, 
                  file = "models/sp_net_vs_times_bor_all")
summary(sp.mod.bor.all) #positive

# predictions
sp.net.mod <- ggpredict(sp.mod.bor.all, terms = "TimesColonised")
colnames(sp.net.mod) = c('TimesColonised', 'fit', 'lwr', 'upr', 'dunno')


# plot
(xyz.sp.bor.all <- ggplot(net.change.sp.class.all, aes(x = TimesColonised, y = SpeciesNetChange)) +
    geom_jitter(alpha = 0.1) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth"))


# FIGURE S1b
(sp_net_plot <- ggplot() +
    geom_point(data = net.change.sp.class.all, aes(x = TimesColonised, y = SpeciesNetChange), 
               alpha = 0.15, size = 10, pch = 20, colour = "#525453") + 
    geom_line(data = sp.net.mod, aes(x = TimesColonised, y = fit), colour = "#525453", linewidth = 2, linetype = "solid") +
    geom_ribbon(data = sp.net.mod, aes(x = TimesColonised, ymin = lwr, ymax = upr), fill = "#525453", alpha = 0.3) +
    xlab("\nTimes colonised") + ylab("Species net change\n") +
    theme(axis.text.x  = element_text(size = 20, colour = "black"), 
          axis.title.x = element_text(face="bold", size = 24),
          axis.text.y = element_text(size = 20, colour = "black"),
          axis.title.y = element_text(face="bold", size = 24),
          panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
          panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.background = element_blank()))


# FIGURE S1 panel
(net.panel <- ggarrange(com_net_plot, sp_net_plot, 
                          labels = c("a)", "b)"), align = "hv",
                          nrow = 1, ncol = 2, font.label = list(size = 28, face = "plain")))

ggsave(net.panel, filename = "figures/net_change_corr_panel.png", 
       width = 55, height = 22, units = "cm")




## SP LOST WITH TRAITS ----
all.traits <- read.csv("data/species_traits_dataset_full.csv") 

# species that were lost at least once
lost.traits <- list.sp.times.ext %>% left_join(., all.traits) %>% 
  filter(ClassNew %in% c("Boreal specialist", "Boreal-tundra boundary"))


# model with traits
lost.traits.mod <- brm(TimesExtinct ~ FunctionalGroup + 
                                    log(PlantHeight) + log(SLA) + log(SeedMass) + log(LeafN), 
                                  data = lost.traits, family = 'negbinomial',
                                  prior = prior(gamma(0.01, 0.01), class = shape), 
                                  iter = 2000, chains = 4, warmup = 400, 
                                  file = "models/lost_traits_mod")
summary(lost.traits.mod) # shorter species are also being lost more often
conditional_effects(lost.traits.mod) # shrubs and graminoids are also lost more often than forbs



# models (all)
sp.extmod <- brm(SpeciesNetChange ~ TimesExtinct,
              data = net.change.spps, 
              iter = 2000, chains = 4, warmup = 400, 
              file = "models/sp_net_vs_times_ext")
summary(sp.extmod) # negative significant 
conditional_effects(sp.extmod)

# model (boreal only)
sp.extmod.bor <- brm(SpeciesNetChange ~ TimesExtinct,
                 data = lost.traits.netch, 
                 iter = 2000, chains = 4, warmup = 400, 
                 file = "models/sp_net_vs_times_ext_bor")
summary(sp.extmod.bor) # negative significant 
conditional_effects(sp.extmod.bor)

# times col vs lost
sp.extcol.mod <- brm(TimesColonised ~ TimesExtinct,
                 data = net.change.spps, 
                 iter = 2000, chains = 4, warmup = 400, 
                 file = "models/sp_col_vs_ext")
summary(sp.extcol.mod) # positive significant 
conditional_effects(sp.extcol.mod)




# species that were lost at least once
lost.traits.netch <- net.change.spps %>% left_join(., all.traits, by = "SPECIES_CLEAN") %>% 
  filter(ClassNew %in% c("Boreal specialist", "Boreal-tundra boundary"))

# model with traits
netchange.traits.mod <- brm(SpeciesNetChange ~ FunctionalGroup + 
                         log(PlantHeight) + log(SLA) + log(SeedMass) + log(LeafN), 
                       data = lost.traits.netch, 
                       iter = 2000, chains = 4, warmup = 400, 
                       file = "models/netchange_traits_mod")
summary(netchange.traits.mod) #ns
conditional_effects(netchange.traits.mod) # ns
