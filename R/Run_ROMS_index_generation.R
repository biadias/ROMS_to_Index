# Purpose
# This code pulls static variables from ROMS NetCDF output and map it to a desired model area. 

## Load files and dependencies
# install.packages("pacman")
pacman::p_load(tidyverse, tidync, sf, rnaturalearth, raster, data.table, maps, mapdata, angstroms, viridis, tabularaster, dtplyr)
select <- dplyr::select

## Read in shape files of desired area.
mask <- st_read("Data/NMFS management area shapefiles/gf95_nmfs.shp")

# subset to NMFS areas of interest
mask <- mask %>% filter(NMFS_AREA %in% c(610,620,630,640,650))%>% # subset to 610-650
  filter(GF95_NMFS1 %in% c(186,194,259,585,870)) %>% # Removes inter-coastal in SE AK
  select(NMFS_AREA, AREA) # area here seems to be in m2

## Import ROMS data
# For GOA, we have grid information stored in a grid file, and the variables stored in the netCDF files.
romsfile_vars <- "Data/ROMS/monthly_averages/nep_hind_moave_2017_01.nc" # read one ROMS file for depth information - depth matching will be done once and not at every time step
romsfile_grid <- "Data/ROMS/NEP_grid_5a.nc"

source("R/ROMS_coordinate_mapping.R") # Maps ROMS coordiantes to NMFS management areas
source("R/ROMS_to_index_functions.R") # Functions to derive indices

netcdf_files <- list.files('Data/ROMS/monthly_averages/', full.names = TRUE)
netcdf_files <- netcdf_files[!grepl('.tmp', netcdf_files)] # gets rid of weird files on loon


## Define variables to pull from ROMS
average_vars <- c('temp', "NH4", "Det", "Eup", "prod_Eup", "Iron", "NCa", "prod_NCa", "MZL", "prod_MZL", "PhL", "prod_PhL", "frat_PhL", "NO3", "salt", "Cop", "prod_Cop", "MZS", "prod_MZS", "PhS", "prod_PhS", "frat_PhS" ) # The variables are averages across depth and space

summed_vars <- average_vars[which(!average_vars %in% c("temp", "salt", "frat_PhS", "frat_PhL"))] # There variables are summed across depths within a depth layer, then averaged across space

# Double check variables are in the ROMS
sum(!average_vars %in% roms_variables$name) # Good if 0

## Summarize 
goa_averaged_vals <- roms_to_goa(netcdf_files = netcdf_files, 
                        variables = average_vars,
                        min_depth = 0, max_depth = -1000,
                        average = TRUE)

goa_summed_vals <- roms_to_goa(netcdf_files = netcdf_files, 
                                 variables = summed_vars,
                                 min_depth = 0, max_depth = -1000,
                                 average = FALSE)

write.csv(goa_averaged_vals, 'nep_avg.csv', row.names = F)
write.csv(goa_summed_vals, 'nep_sum.csv', row.names = F)


# # plot for testing
# goa_averaged_vals %>%
#   filter(varname == "temp") %>%
#   ggplot(aes(x=date,y=value, group = depthclass, color = depthclass))+
#   geom_point()+
#   geom_line()+
#   theme_minimal()+
#   facet_wrap(~NMFS_AREA, scales='free')
