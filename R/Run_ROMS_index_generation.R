# Purpose
# This code pulls static variables from ROMS NetCDF output and map it to a desired model area. 

## Load files and dependencies
# install.packages("pacman")
pacman::p_load(tidyverse, tidync, sf, rnaturalearth, raster, data.table, maps, mapdata, angstroms, viridis, tabularaster, dtplyr)
select <- dplyr::select

## Read in shape files of desired area.
mask <- st_read("Data/Depth trimmed NMFS shapefiles/NMFS610-650.shp")
st_area(mask)/1000^2 # NMFS management area km^2

## Import ROMS data
# For GOA, we have grid information stored in a grid file, and the variables stored in the netCDF files.
romsfile_vars <- "data/ROMS/monthly_averages/nep_hind_moave_2007_01.nc" # read one ROMS file for depth information - depth matching will be done once and not at every time step
romsfile_grid <- "data/ROMS/NEP_grid_5a.nc"

source("R/ROMS_coordinate_mapping.R") # Maps ROMS coordiantes to NMFS management areas
source("R/ROMS_to_index_functions.R") # Functions to derive indices

netcdf_files <- list.files('data/ROMS/monthly_averages/', full.names = TRUE)


## Define variables to pull from ROMS
average_vars <- c('temp', "NH4", "Det", "Eup", "prod_Eup", "Iron", "NCa", "prod_NCa", "MZL", "prod_MZL", "PhL", "prod_PhL", "frat_PhL", "NO3", "salt", "Cop", "prod_Cop", "MZS", "prod_MZS", "PhS", "prod_PhS", "frat_PhS" ) # The variables are averages across depth and space

summed_vars <- average_vars[which(!average_vars %in% c("temp", "salt", "frat_PhS", "frat_PhL"))] # There variables are summed across depths within a depth layer, then averaged across space

# Double check variables are in the ROMS
sum(!average_vars %in% roms_variables$name) # Good if 0

## Summarize 
goa_averaged_vals <- roms_to_goa(netcdf_files = all_files, 
                        variables = average_vars,
                        min_depth = 0, max_depth = -1000,
                        average = TRUE)
# NOTE: about 10% of the ROMS cells are deeper than 1000 m, after masking?

goa_summed_vals <- roms_to_goa(netcdf_files = all_files, 
                                 variables = summed_vars,
                                 min_depth = 0, max_depth = -1000,
                                 average = TRUE)

# plot for testing
goa_averaged_vals %>%
  filter(varname == "temp") %>%
  ggplot(aes(x=date,y=value, group = depthclass, color = depthclass))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  facet_wrap(~NMFS_AREA, scales='free')
