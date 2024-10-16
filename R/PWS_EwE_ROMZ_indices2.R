# --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --- --- --- --- --- --- ---
# Authors: Original A.Rovellini and G. Adams, code modify by B. Ferris and B. Dias
# The code was originally developed for derive indices on average surface and bottom temperature for NMFS mgmt areas 610 to 630 between February and April, since my NMFS areas of interest are 630 and 640.
# Using ROMS NEP subset to the 300 m isocline
##refer to ROMS index generation_BF.qmd for more detailed explanation of code
##this copy generates indices for PWS and GWA areas (region 630,640)
## see variables names in  "NEP variables.csv" 


#Variables in ROMS NEP

library(kableExtra)
nep_vars <- read.csv("Data/NEP_variable_names.csv")
kable(nep_vars)

#Developing an index in R

#Load packages and functions

# install.packages("pacman")
pacman::p_load(mgcv, dplyr, lubridate, ggplot2, tidyr)
source("R/Delta_correction.R")

#Load data
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

nep_ssp245 <- read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp245_300.csv")
nep_ssp245$simulation = "ssp245"

nep_585 <- read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp585_300.csv")
nep_585$simulation = "ssp585"

# Combine in list
roms_avg_data <- do.call(rbind, list(nep_hind, nep_hist, nep_ssp126,nep_ssp245, nep_585))

# Combine in list
roms_sum_data <- do.call(rbind, list(nep_hind, nep_hist, nep_ssp126,nep_ssp245, nep_585))


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

# Run bias correction for all variables 
# - rbinds bias-corrected projection to historical run
# - SSP126
ssp126_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp126"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)

# - SSP245 AVERAGE
ssp245_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp245"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)

# - SSP585 AVERAGE
ssp585_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE)


# --- --- --- --- --- --- ---  --- --- --- --- --- --- --- --- --- --- --- --- 
#### SST indices ####
#Monthly sea surface temperature across 630 and 640 for SSP126 SSP245 and SSP585 

# Extract variables we want
CGOA_SST_ssp126 <- ssp126_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("630","640")) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_630_640 = (value_dc_630 * 98582220025+ value_dc_640*32560976631) / ( 98582220025+32560976631)) # Take area weighted mean

CGOA_SST_ssp245 <- ssp245_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("630","640")) %>%
  mutate(simulation = "ssp245") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_630_640 = (value_dc_630 * 98582220025+ value_dc_640*32560976631) / ( 98582220025+32560976631)) # Take area weighted mean

CGOA_SST_ssp585 <- ssp585_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("630","640")) %>%
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_630_640 = (value_dc_630 * 98582220025+ value_dc_640*32560976631) / ( 98582220025+32560976631)) # Take area weighted mean


# Plot it
CGOA_SST_630_640 <- rbind(CGOA_SST_ssp126, CGOA_SST_ssp245,CGOA_SST_ssp585)
ggplot(goa_SST_630_640, aes(year, value_dc_630_640, colour = simulation)) + geom_line() +
  ylab("SST (Celcius)") + xlab("Year") + facet_wrap(~simulation)

# --- --- --- --- --- --- ---  --- --- --- --- --- --- ---  --- --- --- --- --- --- ---
# Copepod ####
#Mean monthly seasonal small copepod biomass spatial strata 630 and 640 for SSP585 across the water column

# Extract small copepod for desired area
#chose sum (total biomass) vs average
#chose "All" depths instead of surface or other depths
#chose CGOA management areas

#Surface (from 0 to 10 m depth)
#Bottom (deepest ROMS point assumed to be representative of bottom conditions)
#Midwater (from 10 m depth to just above the bottom layer - deepest ROMS point)

#formula to convert mgC/km3 to mt/km2
#metric tonnes/km2=(mgC*m^−3)* h/1000



###### Small Cop 126 ####

CGOA_cop_ssp126_df<- ssp126_biascorrected %>%
  filter(varname == "Cop" & depthclass == "All" & NMFS_AREA %in% c(630, 640)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass
# - Units of Cop are mg C m^-3 of water column
# - Area is in m^2
# - Final units will be tonne 1e-9 tonnes in a gram

#Bia Calculations
#formula to convert mgC/km3 to mt/km2
#metric tonnes/km2=(mgC*m^−3)* h/1000 h = depth (300)

CGOA_cop_ssp126 <- CGOA_cop_ssp126_df %>%
  mutate(Cop_biom = case_when(
    NMFS_AREA == "630" ~  value * (300/1000),# 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "640" ~ value * (300/1000) 
  ))


#Bridget calculations
#CGOA_cop_ssp126_2 <- CGOA_cop_ssp126_df %>%
#  mutate(Cop_biom = case_when(
#    NMFS_AREA == "630" ~ 98582220025 *value * 1e-9,# 620 area = 69583703140 m^2 for 1000m
#    NMFS_AREA == "640" ~ 32560976631 *value * 1e-9 
#  ))


# Monthly average 
CGOA_cop_ssp126 <- CGOA_cop_ssp126 %>%
  group_by(year, month) %>%
  summarise(Cop_biom = sum(Cop_biom))# Sum biomass across 630_640 using Bia calculations

CGOA_cop_ssp126_year <- CGOA_cop_ssp126 %>%
  group_by(year) %>%
  summarise(Cop_biom_mean = mean(Cop_biom))
  
# Save results
write.csv(CGOA_cop_ssp126, file="Output/CGOA_NEP_ROMZ_CopBiom_monthly_630_640_ssp126_300_mt_km2.csv")
write.csv(CGOA_cop_ssp126_year, file="Output/CGOA_NEP_ROMZ_CopBiom_year_630_640_ssp126_300_mt_km2.csv")


# Plot it
ggplot(CGOA_cop_ssp126, aes(year, Cop_biom, colour = month)) + geom_line() +
  ylab("Cop biomass (tonnes*km^-2)") + xlab("Year")+
  geom_smooth(method='lm')

ggplot(CGOA_cop_ssp126_year, aes(year, Cop_biom_mean)) + geom_line() +
  ylab("Cop biomass (tonnes*km^-2)") + xlab("Year")+
  geom_smooth(method='lm')



###### Large Cop 126 ####
CGOA_Lcop_ssp126_df <- ssp126_biascorrected %>%
  filter(varname == "NCa" & depthclass == "All" & NMFS_AREA %in% c(630, 640)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass # See formula above

CGOA_Lcop_ssp126 <- CGOA_Lcop_ssp126_df %>%
  mutate(Cop_biom = case_when(
    NMFS_AREA == "630" ~  value * (300/1000),# 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "640" ~ value * (300/1000) 
  ))

# Monthly  
CGOA_Lcop_ssp126 <- CGOA_Lcop_ssp126 %>%
  group_by(year, month) %>%
  summarise(Cop_biom = sum(Cop_biom))# Sum biomass across 630_640 using Bia calculations

CGOA_Lcop_ssp126_year <- CGOA_Lcop_ssp126 %>%
  group_by(year) %>%
  summarise(Cop_biom_mean = mean(Cop_biom))

# Save results
write.csv(CGOA_Lcop_ssp126, file="CGOA_NEP_ROMZ_LCopBiom_monthly_630_640_ssp126_300_mt_km2.csv")
write.csv(CGOA_Lcop_ssp126_year, file="CGOA_NEP_ROMZ_LCopBiom_yearly_630_640_ssp126_300_mt_km2.csv")

# Plot it
ggplot(CGOA_Lcop_ssp126_year, aes(year, Cop_biom_mean)) + geom_line() +
  ylab("Cop biomass (tonnes)") + xlab("Year")+
  geom_smooth(method="lm")



#--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
### Phytoplankton ####
#Mean monthly phytoplankton biomass spatial strata 630 640 for ssp126 across the water column ta 300m

##### Small Phytoplankton ####
# - SSP126

# Extract small phytoplankton  for desired area
CGOA_PhS_ssp126_df <- ssp126_biascorrected %>%
  filter(varname == "PhS" & depthclass == "All" & NMFS_AREA %in% c(630, 640)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass mt/km2
CGOA_PhS_ssp126 <- CGOA_PhS_ssp126_df %>%
  mutate(PhS_biom = case_when(
    NMFS_AREA == "630" ~  value * (300/1000),# 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "640" ~ value * (300/1000) 
  ))

# Monthly average 
CGOA_PhS_ssp126 <- CGOA_PhS_ssp126 %>%
  group_by(year, month) %>%
  summarise(PhS_biom = sum(PhS_biom))  # Sum biomass across 630 and 640

# Yearly average 
CGOA_PhS_ssp126_year <- CGOA_PhS_ssp126 %>%
  group_by(year) %>%
  summarise(PhS_biom_mean = mean(PhS_biom))

# Save results
write.csv(CGOA_PhS_ssp126, file="CGOA_NEP_ROMZ_PhSBiom_monthly_630_640_ssp126_300_mt_km2.csv")
write.csv(CGOA_PhS_ssp126_year, file="CGOA_NEP_ROMZ_PhSBiom_yearly_630_640_ssp126_300_mt_km2.csv")

# Plot it
ggplot(CGOA_PhS_ssp126_year, aes(year, PhS_biom_mean)) + geom_line() +
  ylab("PhS biomass (tonnes)") + xlab("Year")

##### Large Phytoplankton ####
# - SSP126

# Extract large phytoplankton  for desired area
CGOA_PhL_ssp126_df <- ssp126_biascorrected %>%
  filter(varname == "PhL" & depthclass == "All" & NMFS_AREA %in% c(630,640)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass

CGOA_PhL_ssp126 <- CGOA_PhL_ssp126_df%>%
  mutate(PhL_biom = case_when(
    NMFS_AREA == "630" ~  value * (300/1000),# 620 area = 69583703140 m^2 for 1000m
    NMFS_AREA == "640" ~ value * (300/1000) 
  ))

# Monthly average 
CGOA_PhL_ssp126 <- CGOA_PhL_ssp126 %>%
  group_by(year, month) %>%
  summarise(PhL_biom = sum(PhL_biom))  

#yearly 
CGOA_PhL_ssp126_year <- CGOA_PhL_ssp126 %>%
  group_by(year) %>%
  summarise(PhL_biom_year = mean(PhL_biom))


# Save results
write.csv(CGOA_PhL_ssp126, file="CGOA_NEP_ROMZ_PhLBiom_monthly_630_640_ssp126_300_mt_km2.csv")

# Plot it
ggplot(CGOA_PhL_ssp126, aes(year, PhL_biom, colour = month)) + geom_line() +
  ylab("PhL biomass (tonnes)") + xlab("Year")

### Productivity ####
##### Small Phytoplankton ####

# Extract small phytoplankton  for desired area
CGOA_prod_PhS_ssp126_df <- ssp126_biascorrected %>%
  filter(varname == "prod_PhS" & depthclass == "Surface" & NMFS_AREA %in% c(630, 640)) %>%
  mutate(simulation = "ssp126")

# Convert to biomass mt/km2
## Calculate tonnes km^-2 month^-1
#tonnes_km2_month <- (mg_C_m3_d * z * 1e-9) * 30 * 1e6
# 1) mg_C_m3_d: This is the carbon concentration in mg C m−3 d−1.
# 2) z: The depth in meters ( 10 meters, since I am only using the surface layer 0-10m depth).
# 3) 1e-9: Converts mg to metric tonnes.
# 4) 365: Represents the number of days in a year.
# 5) 1e6: Converts m² to km².


CGOA_prod_PhS_ssp126 <- CGOA_prod_PhS_ssp126_df %>%
  mutate(prod_PhS_biom = case_when(
    NMFS_AREA == "630" ~  (value * 10*1e-9)*365*1e6,
    NMFS_AREA == "640" ~ (value * 10*1e-9)*365*1e6
  ))

# Monthly average 
CGOA_prod_PhS_ssp126 <- CGOA_prod_PhS_ssp126 %>%
  group_by(year, month) %>%
  summarise(prod_PhS_biom_month = sum(prod_PhS_biom))


##### Large Phytoplankton ####

# Extract large phytoplankton production for desired area
CGOA_prod_PhL_ssp126_df <- ssp126_biascorrected %>%
  filter(varname == "prod_PhL" & depthclass == "Surface" & NMFS_AREA %in% c(630, 640)) %>%
  mutate(simulation = "ssp126")
# Convert to biomass mt/km2
CGOA_prod_PhL_ssp126 <- CGOA_prod_PhL_ssp126_df %>%
  mutate(prod_PhL_biom = case_when(
    NMFS_AREA == "630" ~  (value * 10*1e-9)*365*1e6,
    NMFS_AREA == "640" ~ (value * 10*1e-9)*365*1e6
  ))
CGOA_prod_PhL_ssp126 <- CGOA_prod_PhL_ssp126 %>%
  group_by(year, month) %>%
  summarise(prod_PhL_biom_month = sum(prod_PhL_biom))

CGOA_prod_Ph_ssp126 <- left_join(CGOA_prod_PhL_ssp126,CGOA_prod_PhS_ssp126, by= c("year", "month")) %>% 
  rowwise() %>% 
  mutate(CGOA_prod_Ph= mean(c(prod_PhL_biom_month, prod_PhS_biom_month)))


write.csv(CGOA_prod_Ph_ssp126, file="Output/CGOA_NEP_ROMZ_PP_monthly_630_640_ssp126_300_surface_mt_km2.csv")


