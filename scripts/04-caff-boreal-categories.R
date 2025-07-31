## Arctic plant borealization
## Mariana Garcia Criado
## Script 4. Species classification
## May 2024


## LIBRARIES ----
library(tidyverse)



## SPECIES NAMES ----
caff0 <- read.csv("data/ABA_2013_Appendix_9_1.csv", fileEncoding = "latin1")

# Imports the "-" symbol as a error - fix here and remove digits from species names
caff <- caff0 %>% 
  mutate(across(c(2:7), ~ ifelse(.=="\x96", "No", .))) %>%
  rename(Species_name = X) %>%
  mutate(Species_name = str_remove_all(Species_name, "[:digit:]"))

write.csv(caff, "data/ABA_2013_Appendix_9_1_nonumbers.csv")

# Up next is including the full name of the species instead of the abbreviated genus 
# Table format is tricky to work in R, so I've manually fixed the genus for all
# Now need to fix the right name for each species
# file name: ABA_2013_Appendix_9_1_nonumbers_genus
caff.genus <- read.csv("data/ABA_2013_Appendix_9_1_nonumbers_genus.csv")

# remove everything before first period
caff.genus$Species_name_new <- gsub("^.*?\\.", "", caff.genus$Species_name)

# combine genus and species into one column
caff.genus2 <- caff.genus %>% 
  mutate(GENUS_SPECIES = paste0(Genus, Species_name_new)) %>%
  mutate(SPECIES_CLEAN = str_trim(GENUS_SPECIES))

# standardise to species (remove ssp, var)
caff.genus2$SPECIES_CLEAN2 <- gsub("ssp.*", "", caff.genus2$SPECIES_CLEAN)
caff.genus2$SPECIES_CLEAN3 <- gsub("var.*", "", caff.genus2$SPECIES_CLEAN2)
caff.genus2$SPECIES_CLEAN4 <- gsub("\\(.*", "", caff.genus2$SPECIES_CLEAN3)
caff.genus2$SPECIES_CLEAN5 <- gsub("agg.*", "", caff.genus2$SPECIES_CLEAN4)

sort(unique(caff.genus2$SPECIES_CLEAN5))

# clean manually 
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Equisetum "] <- "Equisetum variegatum"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Aruncus dioicus s. lat."] <- "Aruncus dioicus"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Boechera di"] <- "Boechera divaricarpa"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Boschniakia"] <- "Boschniakia rossica"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Carex halophila auct."] <- "Carex halophila"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Oxytropis"] <- "Oxytropis varians"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Papaver"] <- "Papaver variegatum"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Ranunculus auricomus s. lat."] <- "Ranunculus auricomus"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Ranunculus monophyllos s. lat."] <- "Ranunculus monophyllos"
caff.genus2$SPECIES_CLEAN5[caff.genus2$SPECIES_CLEAN5 == "Stellaria longipes s. lat."] <- "Stellaria longipes"

# all ok now
sort(unique(caff.genus2$SPECIES_CLEAN5))

# clean introduced white spaces
caff.genus3 <- caff.genus2 %>% mutate(SPECIES_CLEAN6 = str_trim(SPECIES_CLEAN5)) 

# sort hyphens (they appear as ? but inside a symbol so can't use str_replace)
options(max.print=999999)
sort(unique(caff.genus3$SPECIES_CLEAN6))

caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Alisma plantago�aquatica"] <- "Alisma plantago-aquatica"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Arctostaphylos uva�ursi"] <- "Arctostaphylos uva-ursi"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Athyrium filix�femina"] <- "Athyrium filix-femina"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Braya thorild�wulffii"] <-  "Braya thorild-wulffii"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Capsella bursa�pastoris"] <-  "Capsella bursa-pastoris"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Draba \"pseudo�oxycarpa\""] <-  "Draba pseudo-oxycarpa"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Dryopteris filix�mas"] <-  "Dryopteris filix-mas"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Impatiens noli�tangere"] <-  "Impatiens noli-tangere"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Papaver \"lenaense\""] <-  "Papaver lenaense"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Papaver \"murrayi\""] <-  "Papaver murrayi"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Pedicularis novaiae�zemliae"] <-  "Pedicularis novaiae-zemliae"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Pedicularis sceptrum�carolinum"] <-  "Pedicularis sceptrum-carolinum"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Vaccinium vitis�idaea"] <-  "Vaccinium vitis-idaea"
caff.genus3$SPECIES_CLEAN6[caff.genus3$SPECIES_CLEAN6 == "Salix uva�ursi"] <-  "Salix uva-ursi"

# confirm that all good
sort(unique(caff.genus3$SPECIES_CLEAN6))


# keep relevant columns only, add classification
caff.genus4 <- caff.genus3 %>% 
  dplyr::select(SPECIES_CLEAN6,Arctic.herb.subzone, northern.Arctic.dwarf.shrub.subzone, 
                middle.Arctic.dwarf.shrub.subzone, southern.Arctic.dwarf.shrub.subzone,
                Arctic.shrub.subzone, Adjacent.non.Arctic..Boreal.or.Boreal.alpine..area,
                Species.occurring.in.the.Arctic, Arctic.endemic.species..A.E.,
                Borderline.species..bEN.) %>%
  rename(SPECIES_CLEAN = SPECIES_CLEAN6,
         BorealSpecies = Adjacent.non.Arctic..Boreal.or.Boreal.alpine..area,
         ArcticSpecies = Species.occurring.in.the.Arctic,
         ArcticEndemic = Arctic.endemic.species..A.E.,
         BorderlineSpecies = Borderline.species..bEN.,
         Zone_A = Arctic.herb.subzone,
         Zone_B = northern.Arctic.dwarf.shrub.subzone,
         Zone_C = middle.Arctic.dwarf.shrub.subzone,
         Zone_D = southern.Arctic.dwarf.shrub.subzone,
         Zone_E = Arctic.shrub.subzone) %>%
  #distinct(SPECIES_CLEAN, .keep_all = TRUE) %>%
  mutate(Synonym = SPECIES_CLEAN)

# merge with our own species list
load("data/boreal_master_perm.RData")

sp.list <- boreal_master_perm2 %>% distinct(SPECIES_CLEAN) #351 unique species

# join databse with species occurrences
sp.list.ab <- sp.list %>% left_join(., caff.genus4, by = "SPECIES_CLEAN") # 455 species (incl. dups)

# remove the rows that have exactly the same values across zones
# only 5 records removed - 450 left
sp.list.ab.dis <- sp.list.ab %>% 
  distinct(SPECIES_CLEAN, Zone_A, Zone_B, Zone_C, Zone_D, Zone_E, BorealSpecies, .keep_all = TRUE)



## DUPLICATE SPECIES CATEGORIES ----

# Sort out duplicates
# i.e., those taxa that have the same name at species levels because they had different var/subsp/etc
# we keep the broader distribution for each at the species level

# 67 species to check
sp.dup.n <- sp.list.ab.dis %>% group_by(SPECIES_CLEAN) %>% summarise(n())

# individually choosing the most abundant category per zone
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Achillea millefolium"] <-  "**"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Achillea millefolium"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Achillea millefolium"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Androsace chamaejasme"] <-  "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Androsace chamaejasme"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Androsace chamaejasme"] <-  "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Androsace chamaejasme"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria alpina"] <-  "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria alpina"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria alpina"] <-  "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria alpina"] <-  "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "?"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Antennaria friesiana"] <-  "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Arctagrostis latifolia"] <-  "R"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Arctagrostis latifolia"] <-  "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Arctagrostis latifolia"] <-  "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Arctagrostis latifolia"] <-  "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Arctagrostis latifolia"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Astragalus alpinus"] <-  "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Astragalus alpinus"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Astragalus alpinus"] <-  "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Betula nana"] <-  "S"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Betula pubescens"] <-  "?"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Betula pubescens"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Calamagrostis canadensis"] <-  "?"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Calamagrostis canadensis"] <-  "S"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Calamagrostis canadensis"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Calamagrostis lapponica"] <-  "?"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Calamagrostis lapponica"] <-  "S"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Campanula rotundifolia"] <-  "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Campanula rotundifolia"] <-  "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Campanula rotundifolia"] <-  "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex bigelowii"] <-  "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex bigelowii"] <-  "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex brunnescens"] <- "S"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Carex brunnescens"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Carex brunnescens"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex capillaris"] <- "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex capillaris"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Carex capillaris"] <- "F"

sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Carex microchaeta"] <- "S"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Carex nardina"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex nardina"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex nardina"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Carex nardina"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex parallela"] <- "S"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex parallela"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex saxatilis"] <- "S"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex saxatilis"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Carex scirpoidea"] <- "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Carex scirpoidea"] <- "S"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium alpinum"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium alpinum"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium alpinum"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium alpinum"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium alpinum"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium fontanum"] <- "**"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium fontanum"] <- "*"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Cerastium fontanum"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Diapensia lapponica"] <- "R"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Dryas integrifolia"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Dryas integrifolia"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Dryas integrifolia"] <- "F"

sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Dryas integrifolia"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum angustifolium"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum angustifolium"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum angustifolium"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum angustifolium"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum angustifolium"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum scheuchzeri"] <- "R"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum scheuchzeri"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum scheuchzeri"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum scheuchzeri"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum scheuchzeri"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Eriophorum vaginatum"] <- "S"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Festuca rubra"] <- "S"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Festuca rubra"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Festuca rubra"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Festuca rubra"] <- "F"

sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Gentianella propinqua"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Gentianella propinqua"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Hedysarum boreale"] <- "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Hedysarum boreale"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Hedysarum hedysaroides"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Hedysarum hedysaroides"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Juniperus communis"] <- "R"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Juniperus communis"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Lagotis glauca"] <- "?"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Lagotis glauca"] <- "S"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Lagotis glauca"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Lagotis glauca"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Lagotis glauca"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Linnaea borealis"] <- "?"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Luzula arcuata"] <- "R"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Luzula arcuata"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Luzula multiflora"] <- "S"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Luzula multiflora"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Luzula multiflora"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Luzula parviflora"] <- "S"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Luzula parviflora"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Luzula parviflora"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes hieraciifolia"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes hieraciifolia"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes hieraciifolia"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes nelsoniana"] <- "R"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes nelsoniana"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes nelsoniana"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Micranthes nelsoniana"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Orthilia secunda"] <- "R"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Orthilia secunda"] <- "S"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Orthilia secunda"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Oxytropis deflexa"] <- "?"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Oxytropis deflexa"] <- "R"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Oxytropis deflexa"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Oxytropis deflexa"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Oxytropis maydelliana"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Papaver lapponicum"] <- "S"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Papaver lapponicum"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Papaver lapponicum"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Papaver lapponicum"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Parrya nudicaulis"] <- "R"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Parrya nudicaulis"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Parrya nudicaulis"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Parrya nudicaulis"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Parrya nudicaulis"] <- "F"

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Petasites frigidus"] <- "?"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Petasites frigidus"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Petasites frigidus"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Petasites frigidus"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Petasites frigidus"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Poa alpina"] <- "R"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Poa alpina"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Poa alpina"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Poa alpina"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "S"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Poa arctica"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "R"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Poa pratensis"] <- "F"

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Potentilla hyparctica"] <- "S"
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Potentilla hyparctica"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Potentilla hyparctica"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Potentilla hyparctica"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Potentilla hyparctica"] <- "F"

sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Prunella vulgaris"] <- "b" #b over R

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Pyrola grandiflora"] <- "S" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Pyrola grandiflora"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Pyrola grandiflora"] <- "F"

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Rhododendron lapponicum"] <- "S" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Rhododendron lapponicum"] <- "F" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Rhododendron tomentosum"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Rhododendron tomentosum"] <- "F" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Rumex acetosa"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Rumex acetosa"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Rumex acetosa"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Rumex acetosa"] <- "F" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Rumex aquaticus"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Rumex aquaticus"] <- "S" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Rumex aquaticus"] <- "S" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Salix glauca"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Salix glauca"] <- "S" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Salix glauca"] <- "F" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Salix rotundifolia"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Salix rotundifolia"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Salix rotundifolia"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Salix rotundifolia"] <- "S" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Saussurea angustifolia"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Saussurea angustifolia"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Saussurea angustifolia"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Saussurea angustifolia"] <- "F" 

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F" 
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F"
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F"
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F"
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Saxifraga cespitosa"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Scorzoneroides autumnalis"] <- "?" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Scorzoneroides autumnalis"] <- "F" 

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "R" 
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "F" 
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "F" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Silene involucrata"] <- "F" 

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "S" 
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "F" 
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "F" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "S" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "F"
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Silene uralensis"] <- "F"

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Solidago virgaurea"] <- "R" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Solidago virgaurea"] <- "F" 

sp.list.ab.dis$Zone_A[sp.list.ab.dis$SPECIES_CLEAN == "Trisetum spicatum"] <- "R" 
sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Trisetum spicatum"] <- "F" 
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Trisetum spicatum"] <- "F" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Trisetum spicatum"] <- "F" 

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium uliginosum"] <- "?" 
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium uliginosum"] <- "S" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium uliginosum"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium uliginosum"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium uliginosum"] <- "F" 

sp.list.ab.dis$Zone_B[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium vitis-idaea"] <- "?" 
sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium vitis-idaea"] <- "R" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium vitis-idaea"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Vaccinium vitis-idaea"] <- "F" 

sp.list.ab.dis$Zone_C[sp.list.ab.dis$SPECIES_CLEAN == "Veronica alpina"] <- "S" 
sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Veronica alpina"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Veronica alpina"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Veronica alpina"] <- "F" 

sp.list.ab.dis$Zone_D[sp.list.ab.dis$SPECIES_CLEAN == "Viola biflora"] <- "F" 
sp.list.ab.dis$Zone_E[sp.list.ab.dis$SPECIES_CLEAN == "Viola biflora"] <- "F" 
sp.list.ab.dis$BorealSpecies[sp.list.ab.dis$SPECIES_CLEAN == "Viola biflora"] <- "F" 

# leave one row per species (based on the most abundant subsp/var of that species, as per above)
sp.list.ab.dis.one <- sp.list.ab.dis %>% 
  distinct(SPECIES_CLEAN, Zone_A, Zone_B, Zone_C, Zone_D, Zone_E, BorealSpecies, .keep_all = T)

#back to 351 - perfect!



## SPECIES WITH NO CATEGORIES ----

# first look if it's in the CAFF list under a different name (synonyms)
# if not look in other sources
na.ab.sp <- sp.list.ab.dis.one %>% filter(is.na(Zone_A)) %>% 
  filter(str_detect(SPECIES_CLEAN, 'XXX|xxx', negate = T)) 

# Synonyms checked manually against the Kew plant database 
# first synonym in the CAFF list is taken, also check against Fjallflora
na.ab.sp.fill <- na.ab.sp %>% 
  mutate(Synonym = case_when(SPECIES_CLEAN == "Aconitum septentrionale" ~ "Aconitum lycoctonum", 
                             SPECIES_CLEAN == "Alopecurus magellanicus" ~ "Alopecurus borealis",
                             SPECIES_CLEAN == "Alopecurus alpinus" ~ "Alopecurus borealis",
                             SPECIES_CLEAN == "Androsace ochotensis" ~ "Douglasia ochotensis",
                             SPECIES_CLEAN == "Antennaria pulchella" ~ "Antennaria media",
                             SPECIES_CLEAN == "Arabidopsis lyrata" ~ "Arabidopsis petraea", # can also be ssp Arabidopsis kamtschatica but they are the same category
                             SPECIES_CLEAN == "Braya purpurascens" ~ "Braya glabella",
                             SPECIES_CLEAN == "Calamagrostis stricta" ~ "Calamagrostis neglecta",
                             SPECIES_CLEAN == "Cardamine digitalis" ~ "Cardamine digitata",
                             SPECIES_CLEAN == "Carex fimbriata" ~ "Carex pilulifera",
                             SPECIES_CLEAN == "Carex microcarpa" ~ "Carex aquatilis",
                             SPECIES_CLEAN == "Carex × turfosa" ~ "Carex nigra",
                             SPECIES_CLEAN == "Cerastium glabratum" ~ "Cerastium alpinum",
                             SPECIES_CLEAN == "Dactylorhiza viridis" ~ 	"Coeloglossum viride",
                             SPECIES_CLEAN == "Deschampsia flexuosa" ~ "Avenella flexuosa",
                             SPECIES_CLEAN == "Doronicum grandiflorum" ~ "Not available",
                             SPECIES_CLEAN == "Epilobium latifolium" ~ "Chamerion latifolium",
                             SPECIES_CLEAN == "Eriophorum chamissonis" ~ "Eriophorum russeolum",
                             SPECIES_CLEAN == "Festuca richardsonii" ~"Festuca rubra",
                             SPECIES_CLEAN == "Gnaphalium supinum" ~ "Omalotheca supina",
                             SPECIES_CLEAN == "Hedysarum alpinum" ~ "Not available",
                             SPECIES_CLEAN == "Lagotis minor" ~ "Lagotis glauca",
                             SPECIES_CLEAN == "Lycopodium complanatum" ~ "Diphasiastrum complanatum",
                             SPECIES_CLEAN == "Lysimachia europaea" ~ "Trientalis europaea",
                             SPECIES_CLEAN == "Oreojuncus trifidus" ~ "Juncus trifidus",
                             SPECIES_CLEAN == "Oreomecon radicata" ~ 	"Papaver radicatum",
                             SPECIES_CLEAN == "Oxytropis campestris" ~ "Not available",
                             SPECIES_CLEAN == "Pedicularis sudetica" ~ "Not available",
                             SPECIES_CLEAN == "Persicaria bistorta" ~ "Not available",
                             SPECIES_CLEAN == "Persicaria vivipara" ~ "Bistorta vivipara",
                             SPECIES_CLEAN == "Pilosella officinarum" ~ "Not available",
                             SPECIES_CLEAN == "Potentilla fruticosa" ~ "Dasiphora fruticosa",
                             SPECIES_CLEAN == "Pyrola rotundifolia" ~ "Not available",
                             SPECIES_CLEAN == "Ranunculus lapponicus" ~ "Coptidium lapponicum",
                             SPECIES_CLEAN == "Ranunculus pallasii" ~ "Coptidium pallasii",
                             SPECIES_CLEAN == "Rumex alpestris" ~ "Rumex acetosa",
                             SPECIES_CLEAN == "Salix daphnoides" ~ "Salix pulchra",
                             SPECIES_CLEAN == "Salix lanata x reticulata" ~ "Not available",
                             SPECIES_CLEAN == "Salix uva-ursi" ~ "Salix arbuscula", # could also be myrsinites but same values
                             SPECIES_CLEAN == "Saxifraga flagellaris" ~ "Saxifraga platysepala",
                             SPECIES_CLEAN == "Silene apetala" ~ "Not available",
                             SPECIES_CLEAN == "Stellaria crassipes" ~ "Stellaria longipes",
                             SPECIES_CLEAN == "Vaccinium microcarpum" ~ "Oxycoccus microcarpus"))

# join again with CAFF to extract info for the synonym
na.fill.synonym <- na.ab.sp.fill %>% 
  left_join(caff.genus4, by = "Synonym", relationship = "many-to-many") 

# multiple matches available because of ssp/var, need to manually clean just like above

# do a distinct to remove the rows that have exactly the same values across zones
# 2 rows removed
na.fill.synonym2 <- na.fill.synonym %>% 
  distinct(SPECIES_CLEAN.x, Zone_A.y, Zone_B.y, Zone_C.y, Zone_D.y, Zone_E.y, BorealSpecies.y, .keep_all = TRUE)

# clean manually remaining ones
na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Antennaria pulchella"] <- "R" 
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Antennaria pulchella"] <- "S" 
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Antennaria pulchella"] <- "S"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Antennaria pulchella"] <- "S" 

na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Arabidopsis lyrata"] <- "R" 
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Arabidopsis lyrata"] <- "F" 
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Arabidopsis lyrata"] <- "F" 

na.fill.synonym2$Zone_A.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "R" 
na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "F"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "F"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "F"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Braya purpurascens"] <- "S"

na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Calamagrostis stricta"] <- "?"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Calamagrostis stricta"] <- "S"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Calamagrostis stricta"] <- "F"

na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex × turfosa"] <- "F"

na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex microcarpa"] <- "R"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex microcarpa"] <- "F"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex microcarpa"] <- "F"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex microcarpa"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Carex microcarpa"] <- "F"

na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Cerastium glabratum"] <- "R"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Cerastium glabratum"] <- "F"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Cerastium glabratum"] <- "F"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Cerastium glabratum"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Cerastium glabratum"] <- "F"

na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Festuca richardsonii"] <- "S"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Festuca richardsonii"] <- "F"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Festuca richardsonii"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Festuca richardsonii"] <- "F"

na.fill.synonym2$Zone_A.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lagotis minor"] <- "?"
na.fill.synonym2$Zone_B.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lagotis minor"] <- "S"
na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lagotis minor"] <- "F"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lagotis minor"] <- "F"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lagotis minor"] <- "F"

na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lycopodium complanatum"] <- "?"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Lycopodium complanatum"] <- "R"

na.fill.synonym2$Zone_C.y[na.fill.synonym2$SPECIES_CLEAN.x == "Rumex alpestris"] <- "R"
na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Rumex alpestris"] <- "F"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Rumex alpestris"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Rumex alpestris"] <- "F"

na.fill.synonym2$Zone_D.y[na.fill.synonym2$SPECIES_CLEAN.x == "Salix daphnoides"] <- "R"
na.fill.synonym2$Zone_E.y[na.fill.synonym2$SPECIES_CLEAN.x == "Salix daphnoides"] <- "F"
na.fill.synonym2$BorealSpecies.y[na.fill.synonym2$SPECIES_CLEAN.x == "Salix daphnoides"] <- "F"

# remove the rows that have exactly the same values across zones
# back to 42 species - perfect
na.fill.synonym3 <- na.fill.synonym2 %>% 
  distinct(SPECIES_CLEAN.x, Zone_A.y, Zone_B.y, Zone_C.y, Zone_D.y, Zone_E.y, BorealSpecies.y, .keep_all = TRUE)

# save dataset with synonyms for script 05-trait-extraction
species.list.with.synos <- na.fill.synonym %>% 
  dplyr::select(SPECIES_CLEAN.x, Synonym) %>% 
  rename(SPECIES_CLEAN = SPECIES_CLEAN.x)

write.csv(species.list.with.synos, "data/species_list_with_synos.csv")

# dataset with names and categories
dataset1 <- sp.list.ab.dis.one %>% filter(!is.na(Zone_A)) %>% dplyr::select(-Synonym) # all filled in without NAs

# dataset with most categories (except for 9 species, which will be filled below)
dataset2 <- na.fill.synonym3 %>% dplyr::select(-c("Zone_A.x", "Zone_B.x", "Zone_C.x", "Zone_D.x", "Zone_E.x", 
                                          "BorealSpecies.x", "ArcticSpecies.x", "ArcticEndemic.x",
                                          "BorderlineSpecies.x", "SPECIES_CLEAN.y", "Synonym")) %>% 
  rename(Zone_A = Zone_A.y, Zone_B = Zone_B.y, Zone_C = Zone_C.y, Zone_D= Zone_D.y, Zone_E = Zone_E.y,
         BorealSpecies = BorealSpecies.y, ArcticSpecies = ArcticSpecies.y, ArcticEndemic = ArcticEndemic.y,
         BorderlineSpecies = BorderlineSpecies.y, SPECIES_CLEAN = SPECIES_CLEAN.x)

# check that they have the same variables
names(dataset1)
names(dataset2)

# bind together
bind.all <- rbind(dataset1, dataset2)



## BOREAL/ARCTIC CLASSIFICATION ----

# create acronym based on species occurrence in each zone
bind.acron <- bind.all %>% 
  mutate(Acronym = paste(Zone_A, Zone_B, Zone_C, Zone_D, Zone_E, BorealSpecies, sep = "_"))

sort(unique(bind.acron$Acronym))

# group acronyms into vectors
low.arctic <- c("No_No_No_X_X_No", "No_No_R_F_F_R", "No_No_R_S_F_R")

high.low.arctic <- c("F_F_F_F_R_No", "F_F_F_F_S_R", "No_R_S_S_S_R", "R_F_F_F_S_R", 
                     "S_F_F_F_S_No", "S_F_F_S_R_R", "F_F_F_F_S_No",
                     "No_No_S_F_F_R", "No_R_F_F_F_R")
    
boreal <- c("No_No_No_No_*_F", "No_No_No_No_*_S", "No_No_No_No_b_S", "No_No_No_No_R_F", 
            "No_No_No_No_R_S", "No_No_No_R_R_F", "No_No_?_R_R_F", "No_No_No_No_b_F")
  
high.low.arctic.boreal <- c("?_R_F_F_F_S", "?_S_F_F_F_S", "?_R_F_F_F_F", "?_S_F_F_F_F",
                            "F_F_F_F_F_F",  "F_F_F_F_F_S",
                            "No_No_F_F_F_F", "No_?_F_F_F_F", 
                            "No_No_S_F_F_F", "No_R_S_F_F_F", "No_R_S_F_F_S", "No_R_S_S_S_S",
                            "No_S_F_F_F_F", "No_S_F_F_F_S", "No_No_S_F_F_S", "No_R_F_F_F_F", "No_R_F_F_F_S",  
                            "No_?_S_F_F_F", 
                            "R_F_F_F_F_S", "R_S_F_F_F_S", "R_F_F_F_F_F", 
                            "S_F_F_F_F_F", "S_F_F_F_F_S", "S_F_F_S_F_F")

low.arctic.boreal <- c("No_No_**_**_F_F", "No_No_**_*_F_F", "No_No_**_R_F_F", "No_No_?_?_F_F", "No_No_?_F_F_F",
                       "No_No_?_S_F_F", "No_No_?_S_S_S", "No_No_No_?_F_F", "No_No_No_F_F_F", "No_No_No_No_F_F",
                       "No_No_No_R_F_F", "No_No_No_S_F_F", "No_No_R_F_F_F", "No_No_R_R_F_F", "No_No_R_S_F_F",
                       "No_No_R_S_S_F", "No_No_R_S_S_S", "No_No_No_R_F_S", "No_No_*_*_S_F", "No_No_No_*_S_F", 
                       "No_No_No_?_S_F", "No_No_No_No_S_F", "No_No_No_R_S_F", "No_No_No_R_S_S", "No_No_R_F_F_S", 
                       "No_No_R_R_F_S", "No_No_No_No_R_R", "No_R_R_F_F_S", "No_?_R_F_F_F", "No_No_**_F_F_F", 
                       "No_No_?_R_F_F")
  
#  "NA_NA_NA_NA_NA_NA" species are fixed manually below



# classify species into categories based on the acronyms
bind.acron2 <- bind.acron %>% mutate(Class = case_when(Acronym %in% low.arctic ~ "Low_Arctic",
                                                       Acronym %in% high.low.arctic ~ "High_Low_Arctic",
                                                       Acronym %in% boreal ~ "Boreal",
                                                       Acronym %in% high.low.arctic.boreal ~ "Cosmopolitan",
                                                       Acronym %in% low.arctic.boreal ~ "Low_Arctic_Boreal",
                                                       # fill in NAs from checking distribution in GBIF
                                                       SPECIES_CLEAN == "Persicaria bistorta" ~ "Low_Arctic_Boreal",
                                                       SPECIES_CLEAN == "Doronicum grandiflorum" ~ "Boreal",
                                                       SPECIES_CLEAN == "Hedysarum alpinum" ~ "Low_Arctic_Boreal",
                                                       SPECIES_CLEAN == "Oxytropis campestris" ~ "Low_Arctic_Boreal",
                                                       SPECIES_CLEAN == "Pedicularis sudetica" ~ "Cosmopolitan",
                                                       SPECIES_CLEAN == "Pilosella officinarum" ~ "Boreal",
                                                       SPECIES_CLEAN == "Pyrola rotundifolia" ~ "Low_Arctic_Boreal",
                                                       SPECIES_CLEAN == "Salix lanata x reticulata" ~ "Low_Arctic_Boreal",
                                                       SPECIES_CLEAN == "Silene apetala" ~ "Cosmopolitan"))

# summary per group                                                       
acron.summary <- bind.acron2 %>% group_by(Class) %>% summarise(n = n())

# updated classes
bind.acron3 <- bind.acron2 %>% 
  mutate(ClassNew = case_when(Class %in% c("Low_Arctic", "High_Low_Arctic") ~ "Arctic specialist",
                              Class == "Boreal" ~ "Boreal specialist",
                              Class == "Cosmopolitan" ~ "Ubiquitous",
                              Class == "Low_Arctic_Boreal" ~ "Boreal-tundra boundary"))

# save
write.csv(bind.acron3, "data/species_ab_class_jan2025.csv")

# summary per group                                                       
acron.summ.upd <- bind.acron3 %>% group_by(ClassNew) %>% summarise(n = n())



