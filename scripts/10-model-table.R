## Arctic plant borealization
## Mariana Garcia Criado 
## Model output table
## June 2025

# This script creates Table S3 in the manuscript. All model objects need to be loaded into the session. 


## LIBRARIES ----
library(dplyr)
library(broom)
library(stargazer)
library(stringr)


## FUNCTION ----

# Jonathan Chang's function
p_summarize <- function(model) {
  brms::posterior_summary(model) %>% 
    as_tibble(rownames = "parameter")
}



## SCRIPT 7 MODEL TABLE ----

# Add model objects to list
models.list7 <- list(abun_blab_biogeo_mod_dec2024, 
                     abun_blab_biogeo_mod_all_dec2024,
                     abun_blab_clim_mod_dec2024, 
                     abun_blab_clim_mod_all_dec2024,
                     abun_blab_local_mod2, 
                     abun_blab_local_mod_all,
                     
                     col_blab_biogeo_mod, 
                     col_blab_biogeo_mod_all_bino,
                     col_blab_clim_mod, 
                     col_blab_clim_all_mod_bino,
                     col_blab_local_mod, 
                     col_blab_local_all_mod_bino, 
                     
                     abun.class,
                     species.abun.full.blab.inc.mod,
                     species.abun.blab.inc.fg.mod,
                     species.abun.blab.inc.hei.mod,
                     species.abun.blab.inc.sla.mod,
                     species.abun.inc.blab.seed.mod,
                     species.abun.inc.blab.leafn.mod,
                     
                     col.class.mod2,
                     species.col.full.blab.mod2,
                     species.col.blab.fg.mod2,
                     species.col.blab.hei.mod5,
                     species.col.blab.sla.mod2,
                     species.col.blab.seed.mod2,
                     species.col.blab.leafn.mod2)

# model names
models.name7 <- c("BAI Biogeographical (increases)", 
                     "BAI Biogeographical (full range)",
                     "BAI Climate (increases)", 
                     "BAI Climate (full range)",
                     "BAI Local (increases)", 
                     "BAI Local (full range)",
                     
                     "BCI Biogeographical (increases)", 
                     "BCI Biogeographical (full range)",
                     "BCI Climate (increases)", 
                     "BCI Climate (full range)",
                     "BCI Local (increases)", 
                     "BCI Local (full range)", 
                     
                     "Species abundance change per class",
                     "Species abundance change (multivariate)",
                     "Species abundance change vs functional group (univariate)",
                     "Species abundance change vs height (univariate)",
                     "Species abundance change vs SLA (univariate)",
                     "Species abundance change vs seed mass (univariate)",
                     "Species abundance change vs leaf N (univariate)",
                     
                     "Times colonised per class",
                     "Times clonised (multivariate)",
                     "Times clonised vs functional group (univariate)",
                     "Times clonised vs height (univariate)",
                     "Times clonised vs SLA (univariate)",
                     "Times clonised vs seed mass (univariate)",
                     "Times clonised vs leaf N (univariate)")


# Bind them together
model_number <- 1:26
mod.df <- data.frame(model_number, models.name7)

# Extract parameters
mod.table <- lapply(models.list7, p_summarize) %>% 
  bind_rows(.id = "model_number") 

# Add model name to table
mod.table$model_number <- as.integer(mod.table$model_number)
mod.table.final <- left_join(mod.table, mod.df, by = "model_number")

# Clean model parameters
mod.table.final2 <- mod.table.final %>% 
  dplyr::filter(parameter != "lp__") %>% 
  dplyr::filter(!str_detect(parameter, "^r_"))

mod.table.final2$models.name7[duplicated(mod.table.final2$models.name7)] <- "  "
mod.table.final2$models.name7 <- as.character(mod.table.final2$models.name7)
mod.table.final2$model_number[duplicated(mod.table.final2$model_number)] <- "  "

colnames(mod.table.final2) <- c("Model number", "Term", "Estimate", "Std. error", "Lower 95% CI", "Upper 95% CI", "Model name")
mod.table.final3 <- mod.table.final2[, c(1, 7, 2, 3, 4, 5, 6)]

#  Round to 3 decimals only because not working on stargazer function
mod.table.final4 <- mod.table.final3 %>% mutate_if(is.numeric, round, digits = 3)

# Save in csv
write.csv(mod.table.final4, "models/table_mod_outputs.csv")

