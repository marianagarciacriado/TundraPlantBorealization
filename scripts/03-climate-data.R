## Arctic plant borealization
## Mariana Garcia Criado 
## Script 3. Climate data
## May 2024

# Raw climate files are too large to upload into Zenodo
# Climatic data can be downloaded from https://chelsa-climate.org/ or using the R packages as below
# The cleaned up version of the climate data ready for analysis is produced at the end of this script 
# ('climate_data.csv.csv')


## PACKAGES ----
#devtools::install_github("HelgeJentsch/ClimDatDownloadR")
library(ClimDatDownloadR)
library(raster)
library(broom)


## FUNCTIONS ----
`%notin%` <- Negate(`%in%`)


## DOWNLOAD ----

# this needs a very strong internet connection to download
# better use eduroam
options(timeout = max(4000, getOption("timeout")))

# mean temp (tas)
# errors because temp variables are missing data for 01/1979,
# so skipping to start in 1980 directly
chelsadown <- Chelsa.timeseries.download(
  parameter = "temp",
  start.year.var = 1980,
  start.month.var = 1,
  end.year.var = 2019,
  end.month.var = 12,
  include.month.var = c(1:12),
  version.var = c("2.1"),
  clipping = FALSE,
  clip.shapefile = NULL,
  buffer = 0,
  clip.extent = c(-180, 180, -90, 90))


# tasmin (tmin in package)
chelsadown.tmin <- Chelsa.timeseries.download(
  parameter = "tmin",
  start.year.var = 1980,
  start.month.var = 1,
  end.year.var = 2019,
  end.month.var = 12,
  include.month.var = c(1:12),
  version.var = c("2.1"),
  clipping = FALSE,
  clip.shapefile = NULL,
  buffer = 0,
  clip.extent = c(-180, 180, -90, 90))


# precipitation
# 1979 is available but for consistency with temp let's start in 1980 too.
chelsadown.prec <- Chelsa.timeseries.download(
  save.location = "D:/CHELSA-V2-1-APRIL2024",
  parameter = "prec",
  start.year.var = 1980,
  start.month.var = 1,
  end.year.var = 2018,
  end.month.var = 12,
  version.var = c("2.1"),
  clipping = FALSE,
  clip.shapefile = NULL,
  buffer = 0, 
  clip.extent = c(-180, 180, -90, 90))




## WARMEST QUARTER TEMP ----

## Extracting Warm Quarter data per plot

# defining the filepath
folderpath_mat <- "E:/CHELSA-V2-1-APRIL2024/temp/ChelsaV2.1Timeseries"
filenames_mat <- list.files(folderpath_mat, pattern = "*.tif")
filepath_mat = paste0(folderpath_mat, "/", filenames_mat)

# create raster stack
mat_stack <- stack(filepath_mat)

# add CRS
crs(mat_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract the coordinates as a spatial points object 
load("data/boreal_master_perm.RData")

latlong <- boreal_master_perm2 %>% dplyr::select(LONG, LAT) %>% 
  drop_na() %>% distinct() %>% SpatialPoints()
crs(latlong) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract MAT values for each pair of coordinates
master <- mat_stack %>% raster::extract(., latlong, df = TRUE)

# convert the SpatialPoints object into a dataframe 
latlong2 <- as.data.frame(latlong)

# reassign the ID to the latlong and assign it to the mastersheet too
latlong2$ID <- row.names(latlong2)

# merge the two dataframes
master_mat <- merge(master, latlong2, by = c("ID"))

# reshape from wide to long format and convert to celsius
master_mat2 <- master_mat %>% 
  tidyr::pivot_longer(names_to = "warmq_year", values_to = "value", cols = 2:481) %>%
  separate(warmq_year, c("CHELSA", "variable", "month", "year")) %>% 
  dplyr::select(., -CHELSA) %>%
  mutate(celsius = (value/10) - 273.15) 

check <- master_mat2 %>% filter(celsius < -200) %>% distinct(month, year)

hist(master_mat2$celsius)

# 01-12 1980-1982, and 01-04/1983 have problematic values.
# let's remove the years with problematic values 1980-1983
# and keep warmest quarter only (Jun-Jul-Aug)
master_warmq_year_ts <- master_mat2 %>% 
  filter(year %notin% c(1980:1983)) %>%
  filter(month %in% c("06", "07", "08")) %>%
    dplyr::select(ID, LAT, LONG, variable, year, month, celsius)

write.csv(master_warmq_year_ts, "data/chelsa_warmq_timeseries.csv")

# save as climatology (mean of all years/months per site)
master_warmq_climatology <- master_warmq_year_ts %>%
  group_by(ID) %>%
  mutate(warmq_clim = mean(celsius)) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T) %>%
  dplyr::select(ID, LAT, LONG, warmq_clim)

write.csv(master_warmq_climatology, "data/chelsa_warmq_climatology.csv")



## PRECIPITATION ----

# defining the filepath
folderpath_map <- "E:/CHELSA-V2-1-APRIL2024/prec/ChelsaV2.1Timeseries"
filenames_map <- list.files(folderpath_map, pattern = "*.tif")
filepath_map = paste0(folderpath_map, "/", filenames_map)

# create raster stack
map_stack <- stack(filepath_map)

# add CRS
crs(map_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract MAP values for each pair of coordinates
master.map <- map_stack %>% raster::extract(., latlong, df = TRUE)

# merge the two dataframes
master_map <- merge(master.map, latlong2, by = c("ID"))

# reshape from wide to long format 
master_map2 <- master_map %>% 
  tidyr::pivot_longer(names_to = "prec_year", values_to = "value", cols = 2:481) %>%
  separate(prec_year, c("CHELSA", "variable", "month", "year")) %>% 
  dplyr::select(., -CHELSA) %>%
  mutate(prec = value/100)

hist(master_map2$prec)

# save as yearly timeseries
master_map_year <- master_map2 %>%
  group_by(ID, year) %>%
  mutate(annual_prec = sum(prec)) %>%
  #mutate(mean_annual_prec = mean(prec)) %>%
  ungroup() %>%
  distinct(ID, year, .keep_all = T) %>%
  dplyr::select(ID, LAT, LONG, year, annual_prec)

hist(master_map_year$annual_prec)

# save dataframe as a file
write.csv(master_map_year, "data/chelsa_annual_prec_timeseries.csv")

#master_map_year <- read.csv("data/chelsa_annual_prec_timeseries.csv")

# save as climatology (mean of all years per site)
master_map_climatology <- master_map_year %>%
  group_by(ID) %>%
  mutate(prec_clim = mean(annual_prec)) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T) %>%
  dplyr::select(ID, LAT, LONG, prec_clim)

# save 
write.csv(master_map_climatology, "data/chelsa_prec_climatology.csv")




## MINIMUM TEMP ----

## Extracting min temp data per plot

# defining the filepath
folderpath_min <- "E:/CHELSA-V2-1-APRIL2024/tmin/ChelsaV2.1Timeseries"
filenames_min <- list.files(folderpath_min, pattern = "*.tif")
filepath_min = paste0(folderpath_min, "/", filenames_min)

# create raster stack
min_stack <- stack(filepath_min)

# add CRS
crs(min_stack) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# extract min values for each pair of coordinates
master_min <- min_stack %>% raster::extract(., latlong, df = TRUE)

# merge the two dataframes
master_min_full <- merge(master_min, latlong2, by = c("ID"))

# reshape from wide to long format and convert to celsius
master_min_full_long <- master_min_full %>% 
  tidyr::pivot_longer(names_to = "tmin_year", values_to = "value", cols = 2:481) %>%
  separate(tmin_year, c("CHELSA", "variable", "month", "year")) %>% 
  dplyr::select(., -CHELSA) %>%
  mutate(celsius = (value/10) - 273.15) 

check2 <- master_min_full_long %>% filter(celsius < -200) %>% distinct(month, year)

hist(master_min_full_long$celsius)

# same problem with 4-12/2005, let's chuck this year
# removing the whole year since taking a few months would result in biased minimum temps
master_tmin_year_ts <- master_min_full_long %>% 
  filter(year != 2005) %>%
  dplyr::select(ID, LAT, LONG, variable, month, year, celsius)

write.csv(master_tmin_year_ts, "data/chelsa_tmin_timeseries.csv")


# save as climatology (mean of all years/months per site)
master_tmin_climatology <- master_tmin_year_ts %>%
  group_by(ID) %>%
  mutate(tmin_clim = mean(celsius)) %>%
  ungroup() %>%
  distinct(ID, .keep_all = T) %>%
  dplyr::select(ID, LAT, LONG, tmin_clim)

write.csv(master_tmin_climatology, "data/chelsa_tmin_climatology.csv")




## CHANGE OVER TIME ----
master_warmq_year_ts$year <- as.numeric(master_warmq_year_ts$year)
master_tmin_year_ts$year <- as.numeric(master_tmin_year_ts$year)
master_map_year$year <- as.numeric(master_map_year$year)

# Calculate slopes (annual change in warmest quarter) for each plot 
warmq.slopes <- master_warmq_year_ts %>%
  nest_by(ID) %>%  # 
  mutate(mod = list(lm(celsius ~ year, data = data))) %>% #
  summarise(tidy(mod)) %>%
  filter(term == "year") %>%
  dplyr::select(ID, estimate) %>% 
  rename(WarmQSlope = estimate) %>% unnest()

# Calculate slopes (annual change in min temp) for each plot 
mintemp.slopes <- master_tmin_year_ts %>%
  nest_by(ID) %>%  
  mutate(mod = list(lm(celsius ~ year, data = data))) %>% 
  summarise(tidy(mod)) %>%
  filter(term == "year") %>%
  dplyr::select(ID, estimate) %>% 
  rename(MinTempSlope = estimate) %>% unnest()

# Calculate slopes (annual change in precipitation) for each plot 
prec.slopes <- master_map_year %>%
  nest_by(ID) %>%  
  mutate(mod = list(lm(annual_prec ~ year, data = data))) %>% 
  summarise(tidy(mod)) %>%
  filter(term == "year") %>%
  dplyr::select(ID, estimate) %>% 
  rename(PrecSlope = estimate) %>% unnest()


## CLIMATE FILE ----

# Join with plot data
latlong2$ID <- as.numeric(latlong2$ID)

# First with ID and lat/long
id.slopes <- left_join(latlong2, warmq.slopes, by = "ID")
id.slopes2 <- left_join(id.slopes, prec.slopes, by = "ID")
id.slopes3 <- left_join(id.slopes2, mintemp.slopes, by = "ID")

# Now with plot name data
load("data/boreal_master_perm.RData")

plot <- boreal_master_perm2 %>% dplyr::select(SiteSubsitePlot, LAT, LONG) %>% 
  distinct(SiteSubsitePlot, .keep_all = T)

clim.all <- left_join(plot, id.slopes3, by = c("LAT", "LONG"))

# include climatologies
clim.all2 <- left_join(clim.all, master_tmin_climatology, 
                       by = c("ID", "LAT", "LONG"))

clim.all3 <- left_join(clim.all2, master_warmq_climatology, 
                       by = c("ID", "LAT", "LONG"))

clim.all4 <- left_join(clim.all3, master_map_climatology, 
                       by = c("ID", "LAT", "LONG"))


# save dataframe as a file
write.csv(clim.all4, "data/climate_data.csv")




## Temp change vs latitude ----

# Have locations in higher latitudes in our dataset experienced greater rates of warming? 
# These models can be run by loading the files below

# Prep data for analysis
load("data/boreal_master_perm.RData")

subsite.only <- boreal_master_perm2 %>% 
  distinct(SiteSubsitePlot, .keep_all = T) %>%
  select(SiteSubsitePlot, SiteSubsite)

# all plots within a subsite will have the same values of change
# since (we only have coordinates at the subsite level
clim.all4 <- read.csv("data/climate_data.csv") %>% 
  left_join(., subsite.only)

# fit model at the plot level with subsite random effect
tempch.lat <- brm(bf(WarmQSlope ~ LAT + (1|SiteSubsite)),
                  data = clim.all4, 
                  family = 'beta', 
                  init = 0, cores = 4, 
                  iter = 5000, chains = 4, warmup = 400,
                  control = list(max_treedepth = 12),
                  file = "models/tempch_lat_mod")

# this converged when upping iterations and max_treedepth, though Rhat 1.01
summary(tempch.lat) # positive significant
plot(tempch.lat)
conditional_effects(tempch.lat)

# compare with unique values at subsite level 
clim.unique <- clim.all4 %>% distinct(SiteSubsite, .keep_all = TRUE)

# subsite-level model (smaller sample size, but all plots within a subsite have the same climatological values)
tempch.lat.sb <- brm(bf(WarmQSlope ~ LAT),
                     data = clim.unique, 
                     family = 'beta', 
                     init = 0, 
                     iter = 2000, chains = 4, warmup = 400,
                     file = "models/tempch_lat_sb_mod")

summary(tempch.lat.sb) #positive significant 
