# TundraPlantBorealization
Data and code for the article: García Criado et al. (2025). Borealization of plant communities in the Arctic is driven by boreal-tundra species. Ecology Letters.

# Authors
Mariana García Criado, Isabel C. Barrio, James D. M. Speed, Anne D. Bjorkman, Sarah C. Elmendorf, Isla H. Myers-Smith, Rien Aerts, Juha M. Alatalo, Katlyn R. Betway-May, Robert G. Björk, Mats P. Björkman, Daan Blok, Elisabeth J. Cooper, J. Hans C. Cornelissen, William A. Gould, Ragnhild Gya, Greg H. R. Henry, Luise Hermanutz, Robert D. Hollister, Annika K. Jägerbrand, Ingibjörg S. Jónsdóttir, Elina Kaarlejärvi, Olga Khitun, Simone I. Lang, Petr Macek, Jeremy L. May, Anders Michelsen, Signe Normand, Siri L. Olsen, Eric Post, Riikka Rinnan, Niels Martin Schmidt, Sofie Sjogersten, Anne Tolvanen, Joachim P. Töpper, Andrew Trant, Vigdis Vandvik and Tage Vowles

Contact: Mariana García Criado, mariana.garcia.criado@gmail.com

# Data use guidelines
Data output from this manuscript are publicly available using a Creative Commons Attribution 4.0 International copyright (see license.txt). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis. 

# Data
Climate data from CHELSA can be accessed at https://chelsa-climate.org/ 
Permafrost data from Obu et al. (2019) can be accessed at https://www.sciencedirect.com/science/article/pii/S0012825218305907
Trait data from TRY can be accessed at https://www.try-db.org/TryWeb/Home.php

Plant composition data is included here in the 'data' folder. Script 01-borealiz-data-prep shows the process of generating the main input file ('boreal_master_dec2024.RData', which loads in R as 'itex.full16'), but the original raw files including the whole of the ITEX+ database are not provided in this repository. The full ITEX+ database will be available as part of an upcoming data paper, which will be available here: https://github.com/annebj/ITEX30_VegComp

When a raw data file is not included in this repository, this is specified in the script as "this file is not available in the repo as it contains raw data". There are two reasons for this: either the data have not been made publicly available yet (see above for the ITEX+ database) or the files are too large to upload and can be downloaded directly via their own website (see above for CHELSA, TRY and permafrost data). In all cases, the summarised data files are provided in this repository as input files and enable the reproducibility of all figures and analyses of this manuscript.

# Scripts
All the analyses undertaken for this manuscript are split between multiple R scripts within the scripts folder. They can be followed in a sequential order (i.e., 1 to 11), with the '01-borealiz-data-prep' script showing the process of generating the summarised input files.

# Software requirements
R version 4.2.0. or greater.

R packages: arrow, brms, broom, ClimDatDownloadR, corrplot, ggOceanMaps, ggeffects, ggnewscale, ggpubr, ggrepel, grid, gridExtra, raster, tidyverse, vegan, viridis, WorldFlora

# Data description
An overview of the columns of the main data file of plant composition and associated categories ('boreal_master_dec2024.RData') is included below. Each row is a plant record.
- SiteSubsitePlotYear: the name of the specific plot, within its subsite and site, in that particular year.
- SiteSubsitePlot: the name of the specific plot, within its subsite and site.
- SiteSubsite: the name of the specific subsite within its site.
- SITE: the name of the specific site. 
- SUBSITE: the name of the specific subsite. 
- PLOT: the name of the specific plot.
- YEAR: the year when the plot was surveyed.
- STATUS: indicates whether the identified plant species was alive or dead at the time of surveyed.
- TISSUE: if available, it indicates the part of the plant that was hit during a point-framing survey.
- ValueType: type of surveying, including percent cover and point-framing with and without XY coordinates.
- ORIGINAL_NAME: name with which the species was originally identified.
- GENUS: genus of the species with which it was originally identified.
- GFNARROWwalker: not standardised functional group, including the values SEVER, RUSH, FORB, SDECI, GRASS, SHRUBU, SHRUB, GRAMINOIDU, SEDGE, GRAMU, WOODYU.
- SPECIES_NAME: name of the species after preliminary taxonomic standardisation.
- RelCover: relative % cover of the species abundance within that particular plot and year.
- Moisture: categorical variable indicating the moisture class at the subsite level.
- LAT: latitudinal coordinates at the subsite level.
- LONG: longitudinal coordinates at the subsite level.
- ELEV: elevation above sea level, in meters.
- AZONE: Arctic zone.
- SurveyedArea: plot size, in m2.
- DataContributor: PI of the ITEX subsite.
- DomGrazer: type of dominant grazer.
- GrazerIntensity: intensity of grazing, according to PI expertise.
- family: taxonomic family of the species.
- SPECIES_CLEAN: final species name after full taxonomic check.
- FunctionalGroup: standardised functional group to higher level (shrub, graminoid, forb).
- Woodiness: whether a species is woody or not (1 = woody, 0 = not woody).
- Deciduousness: whether a species is deciduous or not (Evergreen, Deciduous, Not Applicable).
- BerryCat: whether a species produces berries or not (Not berry, Berry).
- N_fixer: whether a species is a Nitrogen fixer or not (0 = not a fixer, 1 = fixer).
- lat_grid: rounded latitudinal value to the nearest decimal.
- lon_grid: rounded longitudinal value to the nearest decimal.
- gridcell: grid cell defined by coordinates.
- Region: biogeographical region identified based on glaciation history.
- Duration: duration of the monitoring per plot, in years (from first to last time point).


