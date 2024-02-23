# Alberto Rovellini and Grant Adams
# Code modify by B. Dias
# The code was originally developed for derive indices on average surface and bottom temperature for NMFS mgmt areas 610 to 630 between February and April, since my NMFS areas of interest are 630 and 640.
# Using ROMS NEP subset to the 300 m isocline

pacman::p_load(mgcv, dplyr, lubridate, ggplot2, tidyr, sf, tidyverse)
source("R/Delta_correction.R")

# ------------------------------------------
# GENERATE SHAPEFILE FOR AREAS 640 AND 630
# ------------------------------------------

#library(rgdal)
#library(raster)
#
##read in first polygon shapefile
#NMFS630 <- readOGR('Data/Depth trimmed NMFS shapefiles 300/NMFS630.shp')
##read in second polygon shapefile
#NMFS640 <- readOGR('Data/Depth trimmed NMFS shapefiles 300/NMFS640.shp')
#NMFS630_640 <- union(NMFS630, NMFS640)
#plot(NMFS630_640, col='brown')
#writeOGR(NMFS630_640, dsn = 'Data/Depth trimmed NMFS shapefiles 300/', layer = 'NMFS630_640',
#         driver = "ESRI Shapefile")
#

# ------------------------------------------
# CALCULATE AREA PER NMFS MGMT AREA
# ------------------------------------------

#make sure raster pkg is not loaded
nmfs <-
  st_read('Data/Depth trimmed NMFS shapefiles 300/NMFS610-650.shp')
plot(nmfs, col = 'brown')

areas <-
  nmfs %>% rowwise() %>% mutate(areas = st_area(geometry)) %>% ungroup()

# area 650 for the 300 m mask is split into two, due to a deep channel. Combine it.
areas <-
  areas %>% group_by(NMFS_AREA) %>% summarise(geometry = st_union(geometry), areas = sum(areas))

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


# ------------------------------------------
# Load data averaged across depth and strata
# ------------------------------------------
nep_hind <-
  read.csv("Data/NEP_10k_revised_indices/nep_avg_hind_300.csv")
nep_hind$simulation = "hindcast"

nep_hist <-
  read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_hist_300.csv")
nep_hist$simulation = "historical"

nep_ssp126 <-
  read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp126_300.csv")
nep_ssp126$simulation = "ssp126"

nep_ssp585 <-
  read.csv("Data/NEP_10k_revised_indices/nep_avg_wb_ssp585_300.csv")
nep_ssp585$simulation = "ssp585"

# Combine in list
roms_avg_data <-
  do.call(rbind, list(nep_hind, nep_hist, nep_ssp126, nep_ssp585))

# Add time and date information
roms_avg_data <- roms_avg_data %>%
  mutate(
    date = lubridate::as_date(date),
    month = lubridate::month(date),
    year = lubridate::year(date)
  )


# ---------------------------------------------------------
# Run bias correction for all variables  (hind and no hind)
# ---------------------------------------------------------
# - - 1) NO HINDCAST SPLICED IN
# - SSP126
ssp126_biascorrected_nohind <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp126"),
  ref_yrs = 2000:2014,
  # Overlap years for historical and hindcast ROMS
  lognormal = FALSE,
  include_hindcast = FALSE
)

# - SSP585
ssp585_biascorrected_nohind <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014,
  # Overlap years for historical and hindcast ROMS
  lognormal = FALSE,
  include_hindcast = FALSE
)

# - - 2) HINDCAST SPLICED IN
# - SSP126
ssp126_biascorrected_hind <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp126"),
  ref_yrs = 2000:2014,
  # Overlap years for historical and hindcast ROMS
  lognormal = FALSE,
  include_hindcast = TRUE
) # Splice hindcast in

# - SSP585
ssp585_biascorrected_hind <- delta_correction(
  hindcast = roms_avg_data %>% filter(simulation == "hindcast"),
  historical = roms_avg_data %>% filter(simulation == "historical"),
  projection = roms_avg_data %>% filter(simulation == "ssp585"),
  ref_yrs = 2000:2014,
  # Overlap years for historical and hindcast ROMS
  lognormal = FALSE,
  include_hindcast = TRUE
) # Splice hindcast in


# - Plot differences between the bias-corrected and raw runs
ssp585_biascorrected_nohind %>%
  select(-unit,-value) %>%
  left_join(
    nep_ssp585 %>%
      mutate(year = year(date), month = month(date)) %>%
      select(NMFS_AREA, depthclass, value, varname, year, month),
    by = c('NMFS_AREA', 'depthclass', 'varname', 'year', 'month')
  ) %>%
  filter(between(year, 2020, 2036)) %>%
  select(-year,-month) %>%
  pivot_longer(c(value, value_dc)) %>%
  filter(varname == 'temp' &
           depthclass %in% c('Surface', 'Bottom') & NMFS_AREA == 630) %>%
  ggplot(aes(x = date, y = value, color = name)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_wrap(depthclass ~ NMFS_AREA)

ssp585_biascorrected_nohind %>%
  select(-unit,-value) %>%
  left_join(
    nep_ssp585 %>%
      mutate(year = year(date), month = month(date)) %>%
      select(NMFS_AREA, depthclass, value, varname, year, month),
    by = c('NMFS_AREA', 'depthclass', 'varname', 'year', 'month')
  ) %>%
  filter(between(year, 2020, 2036)) %>%
  select(-year,-month) %>%
  pivot_longer(c(value, value_dc)) %>%
  filter(varname == 'temp' &
           depthclass %in% c('Surface', 'Bottom') & NMFS_AREA == 640) %>%
  ggplot(aes(x = date, y = value, color = name)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  facet_wrap(depthclass ~ NMFS_AREA)


# Extract variables we want
# ---------------------------------------
# Surface and bottom temperature
# ---------------------------------------
# -- 1) No hindcast
# areas m^2 of 0-300 m shelf by NMFS area
# 610 57225003746
# 620 62597059226
# 630 98582220025
# 640 32560976631
# 650 36726651409

goa_temp_ssp126_nohind <- ssp126_biascorrected_nohind %>%
  filter(
    varname == "temp" &
      depthclass %in% c("Surface", "Bottom") &
      NMFS_AREA %in% c("630", "640")
  ) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value),
              names_from = NMFS_AREA) %>%
  mutate(
    value_dc_630_to_640 = (
      value_dc_630 * 98582220025 + value_dc_640 * 32560976631
    ) / (98582220025 + 32560976631)
  ) # Take area weighted mean

goa_temp_ssp585_nohind <- ssp585_biascorrected_nohind %>%
  filter(
    varname == "temp" &
      depthclass %in% c("Surface", "Bottom")  &
      NMFS_AREA %in% c("630", "640")
  ) %>%
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value),
              names_from = NMFS_AREA) %>%
  mutate(
    value_dc_630_to_640 = (
      value_dc_630 * 98582220025 + value_dc_640 * 32560976631
    ) / (98582220025 + 32560976631)
  ) # Take area weighted mean

goa_temp_630_to_640_nohind <-
  rbind(goa_temp_ssp126_nohind, goa_temp_ssp585_nohind)
goa_temp_630_to_640_nohind <-
  goa_temp_630_to_640_nohind %>% mutate(hind = 'no')


# -- 2) Including hindcast
goa_temp_ssp126_hind <- ssp126_biascorrected_hind %>%
  filter(
    varname == "temp" &
      depthclass %in% c("Surface", "Bottom") &
      NMFS_AREA %in% c("630", "640")
  ) %>%
  mutate(simulation = "ssp126") %>%
  pivot_wider(values_from = c(value_dc, value),
              names_from = NMFS_AREA) %>%
  mutate(
    value_dc_630_to_640 = (
      value_dc_630 * 98582220025 + value_dc_640 * 32560976631
    ) / (98582220025 + 32560976631)
  ) # Take area weighted mean

goa_temp_ssp585_hind <- ssp585_biascorrected_hind %>%
  filter(
    varname == "temp" &
      depthclass %in% c("Surface", "Bottom") &
      NMFS_AREA %in% c("630", "640")
  ) %>%
  mutate(simulation = "ssp585") %>%
  pivot_wider(values_from = c(value_dc, value),
              names_from = NMFS_AREA) %>%
  mutate(
    value_dc_630_to_640 = (
      value_dc_630 * 98582220025 + value_dc_640 * 32560976631
    ) / (98582220025 + 32560976631)
  ) # Take area weighted mean

goa_temp_630_to_640_hind <-
  rbind(goa_temp_ssp126_hind, goa_temp_ssp585_hind)
goa_temp_630_to_640_hind <-
  goa_temp_630_to_640_hind %>% mutate(hind = 'yes')


# -- Combine and plot
goa_temp_630_to_640 <-
  rbind(goa_temp_630_to_640_hind, goa_temp_630_to_640_nohind)

ggplot(goa_temp_630_to_640,
       aes(date, value_dc_630_to_640, colour = hind, alpha= hind)) + geom_line(linewidth=0.7) +
  scale_color_manual(values=c("#56B4E9", "#D55E00"))+
  scale_alpha_manual(values=c(0.8, 0.5))+
  ylab("SST (Celsius)") + xlab("Year") + facet_wrap( ~depthclass+simulation)

goa_temp_630_to_640  <- goa_temp_630_to_640  %>%
  arrange(depthclass, simulation, hind, year, month) %>%
  relocate(hind, .before = varname) %>%
  relocate(simulation, .before = varname)

write.csv(goa_temp_630_to_640 ,
          'Output/goa_temperature_630_to_640_300m.csv',
          row.names = F)


# ---------------------------------------
# Calc mean feb to april temperature
# ---------------------------------------
goa_temp_630_to_640_FebApril <- goa_temp_630_to_640 %>%
  filter(month %in% c(2:4)) %>%
  group_by(year, depthclass, varname, simulation, hind) %>%
  summarise(mean_value_dc_630_to_640 = mean(value_dc_630_to_640)) %>%
  mutate(varname = "Mean temp Feb to April") %>%
  arrange(depthclass, simulation, hind, year) %>%
  select(depthclass,
         hind,
         simulation,
         varname,
         year,
         mean_value_dc_630_to_640)

# -- Plot

ggplot(
  goa_temp_630_to_640_FebApril %>% filter(depthclass == "Surface"),
  aes(year, mean_value_dc_630_to_640, colour = hind, alpha=hind)
) +
  geom_line(size = 2) +
  scale_color_manual(values=c("#56B4E9", "#D55E00"))+
  scale_alpha_manual(values=c(0.8, 0.5))+
  ylab("SST (Celsius)") +
  xlab("Year Surface") +
  facet_wrap( ~simulation)

ggplot(
  goa_temp_630_to_640_FebApril %>% filter(depthclass == "Bottom"),
  aes(year, mean_value_dc_630_to_640, colour = hind, alpha=hind)
) +
  geom_line(size = 2) +
  scale_color_manual(values=c("#56B4E9", "#D55E00"))+
  scale_alpha_manual(values=c(0.8, 0.5))+
  ylab("Bottom temp (Celsius)") +
  xlab("Year Bottom") +
  facet_wrap( ~ simulation)

# -- Save
write.csv(
  goa_temp_630_to_640_FebApril,
  'Output/goa_temp_630_to_640_spring_300M.csv',
  row.names = F
)
