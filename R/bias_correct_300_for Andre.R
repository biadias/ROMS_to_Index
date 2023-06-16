# prepare 0-300 m
library(tidyverse)

pacman::p_load(mgcv, dplyr, lubridate, ggplot2, tidyr, sf)
source("R/Delta_correction.R")

# calc areas
nmfs <- st_read('Data/Depth trimmed NMFS shapefiles 300/NMFS610-650.shp')
areas <- nmfs %>% rowwise() %>% mutate(areas = st_area(geometry)) %>% ungroup()

# area 650 for the 300 m mask is split into two, due to a deep channel. Combine it.
areas <- areas %>% group_by(NMFS_AREA) %>% summarise(geometry = st_union(geometry), areas = sum(areas))

# areas m^2 of 0-300 m shelf by NMFS area
# 610 57225003746
# 620 62597059226
# 630 98582220025
# 640 32560976631
# 650 36726651409

# areas m^2 of 0-1000 m shelf by NMFS area
# 610 63986698621 [m^2]
# 620 69583703140 [m^2]
# 630 105918077937 [m^2]
# 640 37270389681 [m^2]
# 650 43952466109 [m^2]

# Load data averaged across depth and strata
nep_hind <- read.csv("Data/output_from_loon/nep_avg_hind_300.csv")
nep_hind$simulation = "hindcast"

nep_hist <- read.csv("Data/output_from_loon/nep_avg_wb_hist_300.csv")
nep_hist$simulation = "historical"

nep_ssp126 <- read.csv("Data/output_from_loon/nep_avg_wb_ssp126_300.csv")
nep_ssp126$simulation = "ssp126"

nep_ssp585 <- read.csv("Data/output_from_loon/nep_avg_wb_ssp585_300.csv")
nep_ssp585$simulation = "ssp585"

# Combine in list
roms_avg_data <- do.call(rbind, list(nep_hind, nep_hist, nep_ssp126, nep_ssp585))

# Add time and date information
roms_avg_data <- roms_avg_data %>%
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
  lognormal = FALSE,
  include_hindcast = FALSE)

# - SSP585
ssp585_biascorrected <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014, # Overlap years for historical and hindcast ROMS
  lognormal = FALSE,
  include_hindcast = FALSE)

# # plot differences between the bias-corrected and raw runs
# ssp585_biascorrected %>%
#   select(-unit, -value) %>%
#   left_join(nep_ssp585 %>% 
#               mutate(year = year(date), month = month(date)) %>%
#               select(NMFS_AREA, depthclass, value, varname, year, month),
#             by = c('NMFS_AREA', 'depthclass', 'varname', 'year', 'month')) %>%
#   filter(between(year, 2035, 2036)) %>%
#   select(-year, -month) %>%
#   pivot_longer(c(value, value_dc)) %>%
#   filter(varname == 'temp' & depthclass %in% c('Surface', 'Bottom') & NMFS_AREA == 610) %>%
#   ggplot(aes(x = date, y = value, color = name))+
#   geom_line()+
#   geom_point()+
#   theme_bw()+
#   facet_wrap(depthclass~NMFS_AREA)

# Extract variables we want
# surface
goa_SST_ssp126 <- ssp126_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("610", "620")) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226) / (57225003746 + 62597059226)) # Take area weighted mean

goa_SST_ssp585 <- ssp585_biascorrected %>%
  filter(varname == "temp" & depthclass == "Surface" & NMFS_AREA %in% c("610", "620")) %>% 
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226) / (57225003746 + 62597059226)) # Take area weighted mean


# Plot it
goa_SST_610_620_surface <- rbind(goa_SST_ssp126, goa_SST_ssp585)
ggplot(goa_SST_610_620_surface, aes(date, value_dc_610_620, colour = simulation)) + geom_line() +
  ylab("SST (Celsius)") + xlab("Year") + facet_wrap(~simulation)

write.csv(goa_SST_610_620_surface, 'goa_temperature_610_620_surface_300m_NoHindcast.csv', row.names = F)

# bottom
goa_SST_ssp126 <- ssp126_biascorrected %>%
  filter(varname == "temp" & depthclass == "Bottom" & NMFS_AREA %in% c("610", "620")) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226) / (57225003746 + 62597059226)) # Take area weighted mean

goa_SST_ssp585 <- ssp585_biascorrected %>%
  filter(varname == "temp" & depthclass == "Bottom" & NMFS_AREA %in% c("610", "620")) %>% 
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value), names_from = NMFS_AREA) %>%
  mutate(value_dc_610_620 = (value_dc_610 * 57225003746 + value_dc_620 * 62597059226) / (57225003746 + 62597059226)) # Take area weighted mean


# Plot it
goa_SST_610_620_bottom <- rbind(goa_SST_ssp126, goa_SST_ssp585)
ggplot(goa_SST_610_620_bottom, aes(date, value_dc_610_620, colour = simulation)) + geom_line() +
  ylab("SST (Celcius)") + xlab("Year") + facet_wrap(~simulation)

write.csv(goa_SST_610_620_bottom, 'goa_temperaure_610_620_bottom_300m_NoHindcast.csv', row.names = F)


# Check against Grant's SST file ------------------------------------------

# dat_grant <- read.csv('610_620_SST.csv')
# dat_albi <- read.csv('goa_SST_610_620_surface_1000.csv') # try full data first
# 
# test <- data.frame('g' = dat_grant$value_dc_610_620, 'a' = dat_albi$value_dc_610_620)
# test %>% mutate(diff = g - a) %>% pull(diff) %>% abs() %>% max() # 4.9919e-09
# 
# ggplot(dat_albi %>% filter(year %in% c(2014:2016)), aes(x = as.Date(date), y = value_dc_610_620))+
#   geom_line()+
#   facet_wrap(~simulation)


# Check with and without hindcast baked in --------------------------------

dat_hind <- read.csv('goa_temperaure_610_620_bottom_300m_Hindcast.csv')
dat_no_hind <- read.csv('goa_temperaure_610_620_bottom_300m_NoHindcast.csv')

dat_hind <- dat_hind %>% mutate(hind = 'yes')
dat_no_hind <- dat_no_hind %>% mutate(hind = 'no')

datall <- rbind(dat_hind, dat_no_hind)

datall %>%
  filter(between(year, 2010, 2020) & simulation == 'ssp585') %>%
  ggplot()+
  geom_line(aes(x = as.Date(date), y = value_dc_610_620, color = hind))+
  theme_bw()
  
