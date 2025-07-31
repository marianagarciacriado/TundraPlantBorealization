## Arctic plant borealization
## Mariana Garcia Criado & Sarah Elmendorf
## Script 8. Table S1 classification
## May 2024


## LIBRARIES ----
library(tidyverse)
library(grid)
library(gridExtra)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)

df <- read.csv("data/species_ab_class_jan2025.csv") %>% distinct(Acronym, .keep_all = T)

# sets how many rows per plot
species_per_group <- 25


df <- read.csv("data/species_ab_class_jan2025.csv") %>% 
  distinct(Acronym, .keep_all = T) %>%
  # pivot longer
  dplyr::select(-c(X, SPECIES_CLEAN)) %>%
  pivot_longer(cols = c(Zone_A:BorealSpecies), names_to = "Zone", values_to = "value") %>%
  # relevel the factors of zone
  mutate(Zone = fct_relevel(
    .$Zone, "BorealSpecies", "Zone_E", "Zone_D", "Zone_C", "Zone_B", "Zone_A")) %>%
  mutate(Zone = fct_relevel(
      .$Zone, "BorealSpecies", "Zone_E","Zone_D", "Zone_C", "Zone_B", "Zone_A" 
    ),
    #SPECIES_CLEAN = fct_rev(SPECIES_CLEAN),
    value = fct_relevel(.$value, "No", "X", "?", "**", "*", "R", "b", "S", "F")
  ) %>%
  filter(Acronym %notin% "NA_NA_NA_NA_NA_NA")

# clean up and order
df2 <- df %>% mutate(ZoneShort = case_when(Zone == "BorealSpecies" ~ "Bor",
                                            Zone == "Zone_E" ~ "E",
                                            Zone == "Zone_D" ~ "D",
                                            Zone == "Zone_C" ~ "C",
                                            Zone == "Zone_B" ~ "B",
                                            Zone == "Zone_A" ~ "A")) %>%
  mutate(ZoneShort = fct_relevel(.$ZoneShort, "Bor", "E", "D", "C", "B", "A")) 


# use a named vector for colors
colour_code_vector <- c("F" = "#fde725",
                 "S" = "#b8de29",
                 "b" = "#56c667",
                 "R" = "#3cbc75",
                 "*" = "#1f968b",
                 "**" = "#1f968b",
                 "?" = "#277d8e",
                 "X" = "#2d708e",
                 "No" = "#6b5fb0")



## BOREAL ----
df.bor <- df2 %>% filter(ClassNew == "Boreal specialist")

(df.bor.plot <- ggplot(df.bor, aes(x = ZoneShort, y = Acronym, fill = value)) +
  geom_tile(color = "white", width = 1, height = 1) +
  geom_text(aes(label = value), color = "black", size = 8) + # Add cell values as text
  scale_fill_manual(values = colour_code_vector) +
  theme_minimal() +
  coord_fixed() + 
  scale_x_discrete(position = "top") +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 20, colour = "black")
  ))

ggsave(df.bor.plot, filename = "figures/table_s3_b.png", 
       width = 10, height = 30, units = "cm")


## BOREAL-TUNDRA ----

df.btb <- df2 %>% filter(ClassNew == "Boreal-tundra boundary")

(df.btb.plot <- ggplot(df.btb, aes(x = ZoneShort, y = Acronym, fill = value)) +
    geom_tile(color = "white", width = 1, height = 1) +
    geom_text(aes(label = value), color = "black", size = 8) + # Add cell values as text
    scale_fill_manual(values = colour_code_vector) +
    theme_minimal() +
    coord_fixed() + 
    scale_x_discrete(position = "top") +
    theme(
      legend.position = "none", axis.title.x = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
      axis.text.x = element_text(size = 20, colour = "black")
    ))

ggsave(df.btb.plot, filename = "figures/table_s3_btb.png", 
       width = 10, height = 40, units = "cm")

## ARCTIC ----
df.a <- df2 %>% filter(ClassNew == "Arctic specialist")

(df.a.plot <- ggplot(df.a, aes(x = ZoneShort, y = Acronym, fill = value)) +
    geom_tile(color = "white", width = 1, height = 1) +
    geom_text(aes(label = value), color = "black", size = 8) + # Add cell values as text
    scale_fill_manual(values = colour_code_vector) +
    theme_minimal() +
    coord_fixed() + 
    scale_x_discrete(position = "top") +
    theme(legend.position = "none", axis.title.x = element_blank(),
  axis.title.y = element_blank(), axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
  axis.text.x = element_text(size = 20, colour = "black")))

ggsave(df.a.plot, filename = "figures/table_s3_arct.png", 
       width = 10, height = 40, units = "cm")


## UBIQUITOUS ----
df.ub <- df2 %>% filter(ClassNew == "Ubiquitous")

(df.ub.plot <- ggplot(df.ub, aes(x = ZoneShort, y = Acronym, fill = value)) +
    geom_tile(color = "white", width = 1, height = 1) +
    geom_text(aes(label = value), color = "black", size = 8) + # Add cell values as text
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = colour_code_vector) +
    theme_minimal() +
    coord_fixed() + 
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 20, colour = "black")))

ggsave(df.ub.plot, filename = "figures/table_s3_ubiq.png", 
       width = 10, height = 40, units = "cm")


