# Plot spatial objects for GOA REOMS indices

pacman::p_load(dplyr, sf, raster, devtools, ggplot2, maps, mapdata)
select <- dplyr::select

# NMFS areas original
nmfs <- st_read("Data/NMFS management area shapefiles/gf95_nmfs.shp")
nmfs <- nmfs %>% filter(NMFS_AREA %in% c(610,620,630,640,650)) %>% # subset to 610-650
  filter(GF95_NMFS1 %in% c(186,194,259,585,870)) # Removes inter-coastal in SE AK

# NMFS areas clipped to 1000 m isobath
nmfs_clipped_1000 <- st_read("Data/Depth trimmed NMFS shapefiles/NMFS610-650.shp")

# NMFS areas clipped to 300 m isobath
# This is what we currently use for the indices
nmfs_clipped_300 <- st_read("Data/Depth trimmed NMFS shapefiles 300/NMFS610-650.shp")

# coast
coast <- maps::map(database = 'worldHires', regions = c('USA','Canada'), plot = FALSE, fill=TRUE)

coast_sf <- coast %>% 
  st_as_sf(crs = 4326) %>% 
  st_transform(crs = st_crs(nmfs_clipped_1000)) %>% 
  st_combine() %>%
  st_crop(nmfs_clipped_1000 %>% st_bbox())

# view
ggplot()+
  geom_sf(data = nmfs, color = 'red', fill = NA)+
  geom_sf(data = nmfs_clipped_1000, color = 'blue', fill = NA)+
  geom_sf(data = nmfs_clipped_300, color = 'orange', fill = NA)+
  geom_sf(data = coast_sf, fill = 'grey')+
  theme_bw()


# Filter for area 650
nmfs_clipped_300_650 <- nmfs_clipped_300 %>%
  filter(NMFS_AREA == 650)

nmfs_clipped_1000_650 <- nmfs_clipped_1000 %>%
  filter(NMFS_AREA == 650)

nmfs_650 <- nmfs %>%
  filter(NMFS_AREA == 650)

# view
ggplot()+
  # geom_sf(data = nmfs_650, color = 'red', fill = NA)+
  #geom_sf(data = nmfs_clipped_1000_650, color = 'blue', fill = NA)+
  geom_sf(data = nmfs_clipped_300_650, color = 'orange', fill = NA)+
  geom_sf(data = coast_sf, fill = 'grey')+
  theme_bw()

  
# remember that these are for visualization purpose - the actual indices are computed over ROMS depth