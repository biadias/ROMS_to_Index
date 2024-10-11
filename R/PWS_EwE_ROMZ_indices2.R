#-----------------------------------------------------------
# Authors: Original A.Rovellini and G. Adams, code modify by B. Ferris and B. Dias
# The code was originally developed for derive indices on average surface and bottom temperature for NMFS mgmt areas 610 to 630 between February and April, since my NMFS areas of interest are 630 and 640.
# Using ROMS NEP subset to the 300 m isocline
##refer to ROMS index generation_BF.qmd for more detailed explanation of code
##this copy generates indices for PWS and GWA areas (region 630,640)
## see variables names in  "NEP variables.csv" 


####Variables in ROMS NEP

library(kableExtra)
nep_vars <- read.csv("Data/NEP_variable_names.csv")
kable(nep_vars)

####Developing an index in R

####Load packages and functions

# install.packages("pacman")
pacman::p_load(mgcv, dplyr, lubridate, ggplot2, tidyr)
source("R/Delta_correction.R")

####Load data
#234769 km2 Bridget Area
#218404 km2 ROMS 

#the discrepancy between the areas is because WGOA area is from 0 to 1000m depth
# Check the PhL PhS Bottom (include or not include)

# areas m^2 of 0-300 m shelf by NMFS area
# 610 57225003746
# 620 62597059226
# 630 98582220025
# 640 32560976631
# 650 36726651409

# Load data averaged across depth and strata
nep_hind <- read.csv("Data/NEP_10k_revised_indices/nep_avg_hind_300.csv")
nep_hind$simulation = "hindcast"

nep_hist <- read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_hist_300.csv")
nep_hist$simulation = "historical"

nep_ssp126 <- read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp126_300.csv")
nep_ssp126$simulation = "ssp126"

nep_585 <- read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp585_300.csv")
nep_585$simulation = "ssp585"

# Combine in list
roms_avg_data <- do.call(rbind, list(nep_hind, nep_hist, nep_ssp126, nep_585))

# Combine in list
roms_sum_data <- do.call(rbind, list(nep_hind, nep_hist, nep_ssp126, nep_585))


# Add time and date information
roms_avg_data <- roms_avg_data %>%
  mutate(
    date = lubridate::as_date(date),
    month = lubridate::month(date),
    year = lubridate::year(date))

# Add time and date information
roms_sum_data <- roms_sum_data %>%
  mutate(
    date = lubridate::as_date(date),
    month = lubridate::month(date),
    year = lubridate::year(date))

#### Run bias correction for all variables 
# - rbinds bias-corrected projection to historical run
# - SSP126
ssp126_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp126"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)

# - SSP585 AVERAGE
ssp585_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)

# - SSP585 SUM
ssp585_sum_biascorrected <- delta_correction(
  hindcast = roms_sum_data %>% filter(simulation == "hindcast"),
  historical = roms_sum_data %>% filter(simulation == "historical"),
  projection = roms_sum_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)

# - SSP126 SUM
ssp126_sum_biascorrected <- delta_correction(
  hindcast = roms_sum_data %>% filter(simulation == "hindcast"),
  historical = roms_sum_data %>% filter(simulation == "historical"),
  projection = roms_sum_data %>% filter(simulation == "ssp126"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)
#######################################################################################
#### SST indices : Monthly sea surface temperature across 610, 620, 630 for SSP126 and SSP585

# Extract variables we want
goa_SST_ssp126 <- ssp126_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("610", "620", "630")) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620_630 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226 + value_dc_630 * 98582220025) / (57225003746 + 62597059226 + 98582220025)) # Take area weighted mean

goa_SST_ssp585 <- ssp585_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("610", "620", "630")) %>% 
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620_630 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226 + value_dc_630 * 98582220025) / (57225003746 + 62597059226 + 98582220025)) # Take area weighted mean


# Plot it
goa_SST_610_620_630 <- rbind(goa_SST_ssp126, goa_SST_ssp585)
ggplot(goa_SST_610_620_630, aes(year, value_dc_610_620_630, colour = simulation)) + geom_line() +
  ylab("SST (Celcius)") + xlab("Year") + facet_wrap(~simulation)

#######################################################################################
#### Mean monthly seasonal small copepod biomass spatial strata 610 and 620 for SSP585 across the water column

# Extract small copepod for desired area
#chose sum (total biomass) vs average
#chose "All" depths instead of surface or other depths
#chose WGOA 3 management areas
goa_cop_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "Cop" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_cop_ssp585 <- goa_cop_ssp585 %>%
  mutate(Cop_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_cop_ssp585 <- goa_cop_ssp585 %>%
  group_by(year, month) %>%
  summarise(Cop_biom = sum(Cop_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_cop_ssp585, file="GOA_NEP_ROMZ_CopBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_cop_ssp585, aes(year, Cop_biom, colour = month)) + geom_line() +
  ylab("Cop biomass (tonnes)") + xlab("Year")




# Extract small copepod for desired area
goa_cop_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "Cop" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_cop_ssp126 <- goa_cop_ssp126 %>%
  mutate(Cop_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_cop_ssp126 <- goa_cop_ssp126 %>%
  group_by(year, month) %>%
  summarise(Cop_biom = sum(Cop_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_cop_ssp126, file="GOA_NEP_ROMZ_CopBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_cop_ssp126, aes(year, Cop_biom, colour = month)) + geom_line() +
  ylab("Cop biomass (tonnes)") + xlab("Year")


#######################################################################################
#### Mean monthly seasonal large copepod biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extractlarge copepod for desired area
goa_NCa_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "NCa" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_NCa_ssp585 <- goa_NCa_ssp585 %>%
  mutate(NCa_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_NCa_ssp585 <- goa_NCa_ssp585 %>%
  group_by(year, month) %>%
  summarise(NCa_biom = sum(NCa_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_NCa_ssp585, file="GOA_NEP_ROMZ_NCaBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_NCa_ssp585, aes(year, NCa_biom, colour = month)) + geom_line() +
  ylab("NCa biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract large copepod for desired area
goa_NCa_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "NCa" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of NCa are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_NCa_ssp126 <- goa_NCa_ssp126 %>%
  mutate(NCa_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_NCa_ssp126 <- goa_NCa_ssp126 %>%
  group_by(year, month) %>%
  summarise(NCa_biom = sum(NCa_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_NCa_ssp126, file="GOA_NEP_ROMZ_NCaBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_NCa_ssp126, aes(year, NCa_biom, colour = month)) + geom_line() +
  ylab("NCa biomass (tonnes)") + xlab("Year")

#######################################################################################
#### Mean monthly seasonal euphausiid biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extract euphausiid for desired area
goa_Eup_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "Eup" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_Eup_ssp585 <- goa_Eup_ssp585 %>%
  mutate(Eup_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_Eup_ssp585 <- goa_Eup_ssp585 %>%
  group_by(year, month) %>%
  summarise(Eup_biom = sum(Eup_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_Eup_ssp585, file="GOA_NEP_ROMZ_EupBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_Eup_ssp585, aes(year, Eup_biom, colour = month)) + geom_line() +
  ylab("Eup biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract  euphausiid for desired area
goa_Eup_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "Eup" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of Eup are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_Eup_ssp126 <- goa_Eup_ssp126 %>%
  mutate(Eup_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_Eup_ssp126 <- goa_Eup_ssp126 %>%
  group_by(year, month) %>%
  summarise(Eup_biom = sum(Eup_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_Eup_ssp126, file="GOA_NEP_ROMZ_EupBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_Eup_ssp126, aes(year, Eup_biom, colour = month)) + geom_line() +
  ylab("Eup biomass (tonnes)") + xlab("Year")

#######################################################################################
#### Mean monthly seasonal large microzooplankton biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extract large microzooplankton for desired area
goa_MZL_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "MZL" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_MZL_ssp585 <- goa_MZL_ssp585 %>%
  mutate(MZL_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_MZL_ssp585 <- goa_MZL_ssp585 %>%
  group_by(year, month) %>%
  summarise(MZL_biom = sum(MZL_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_MZL_ssp585, file="GOA_NEP_ROMZ_MZLBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_MZL_ssp585, aes(year, MZL_biom, colour = month)) + geom_line() +
  ylab("MZL biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract large microzooplankton for desired area
goa_MZL_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "MZL" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of MZL are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_MZL_ssp126 <- goa_MZL_ssp126 %>%
  mutate(MZL_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_MZL_ssp126 <- goa_MZL_ssp126 %>%
  group_by(year, month) %>%
  summarise(MZL_biom = sum(MZL_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_MZL_ssp126, file="GOA_NEP_ROMZ_MZLBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_MZL_ssp126, aes(year, MZL_biom, colour = month)) + geom_line() +
  ylab("MZL biomass (tonnes)") + xlab("Year")

#######################################################################################
#### Mean monthly seasonal small microzooplankton biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extract small microzooplankton for desired area
goa_MZS_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "MZS" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_MZS_ssp585 <- goa_MZS_ssp585 %>%
  mutate(MZS_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_MZS_ssp585 <- goa_MZS_ssp585 %>%
  group_by(year, month) %>%
  summarise(MZS_biom = sum(MZS_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_MZS_ssp585, file="GOA_NEP_ROMZ_MZSBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_MZS_ssp585, aes(year, MZS_biom, colour = month)) + geom_line() +
  ylab("MZS biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract small microzooplankton for desired area
goa_MZS_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "MZS" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of MZS are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_MZS_ssp126 <- goa_MZS_ssp126 %>%
  mutate(MZS_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_MZS_ssp126 <- goa_MZS_ssp126 %>%
  group_by(year, month) %>%
  summarise(MZS_biom = sum(MZS_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_MZS_ssp126, file="GOA_NEP_ROMZ_MZSBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_MZS_ssp126, aes(year, MZS_biom, colour = month)) + geom_line() +
  ylab("MZS biomass (tonnes)") + xlab("Year")

#######################################################################################
#### Mean monthly seasonal small phytoplankton biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extract small phytoplankton for desired area
goa_PhS_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "PhS" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_PhS_ssp585 <- goa_PhS_ssp585 %>%
  mutate(PhS_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_PhS_ssp585 <- goa_PhS_ssp585 %>%
  group_by(year, month) %>%
  summarise(PhS_biom = sum(PhS_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_PhS_ssp585, file="GOA_NEP_ROMZ_PhSBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_PhS_ssp585, aes(year, PhS_biom, colour = month)) + geom_line() +
  ylab("PhS biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract small phytoplankton  for desired area
goa_PhS_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "PhS" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of PhS are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_PhS_ssp126 <- goa_PhS_ssp126 %>%
  mutate(PhS_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_PhS_ssp126 <- goa_PhS_ssp126 %>%
  group_by(year, month) %>%
  summarise(PhS_biom = sum(PhS_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_PhS_ssp126, file="GOA_NEP_ROMZ_PhSBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_PhS_ssp126, aes(year, PhS_biom, colour = month)) + geom_line() +
  ylab("PhS biomass (tonnes)") + xlab("Year")

#######################################################################################
#### Mean monthly seasonal large phytoplankton biomass spatial strata 610, 620, 630 for SSP585 and ssp126 across the water column ta 300m

### ssp585
# Extract large phytoplankton for desired area
goa_PhL_ssp585 <- ssp585_sum_biascorrected %>%
  filter(varname == "PhL" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp585")

# Convert to biomass
# - Units of Cop are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_PhL_ssp585 <- goa_PhL_ssp585 %>%
  mutate(PhL_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_PhL_ssp585 <- goa_PhL_ssp585 %>%
  group_by(year, month) %>%
  summarise(PhL_biom = sum(PhL_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_PhL_ssp585, file="GOA_NEP_ROMZ_PhLBiom_monthly_610t0630_ssp585_300.csv")

# Plot it
ggplot(goa_PhL_ssp585, aes(year, PhL_biom, colour = month)) + geom_line() +
  ylab("PhL biomass (tonnes)") + xlab("Year")

### - SSP126

# Extract large phytoplankton  for desired area
goa_PhL_ssp126 <- ssp126_sum_biascorrected %>%
  filter(varname == "PhL" & depthclass == "All" & NMFS_AREA %in% c(610, 620, 630)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of PhL are mg C m^-2 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram
goa_PhL_ssp126 <- goa_PhL_ssp126 %>%
  mutate(PhL_biom = case_when(
    NMFS_AREA == "610" ~ 57225003746 * value * 1e-9,
    # 610 area = 63986698621 m^2 for 1000m
    NMFS_AREA == "620" ~ 62597059226 * value * 1e-9, # 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "630" ~ 98582220025 * value * 1e-9 # 
  ))

# Monthly average 
goa_PhL_ssp126 <- goa_PhL_ssp126 %>%
  group_by(year, month) %>%
  summarise(PhL_biom = sum(PhL_biom))  # Sum biomass across 610 and 620 and 630

# Save results
write.csv(goa_PhL_ssp126, file="GOA_NEP_ROMZ_PhLBiom_monthly_610t0630_ssp126_300.csv")

# Plot it
ggplot(goa_PhL_ssp126, aes(year, PhL_biom, colour = month)) + geom_line() +
  ylab("PhL biomass (tonnes)") + xlab("Year")

