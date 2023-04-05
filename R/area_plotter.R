# Plot spatial objects for GOA REOMS indices

pacman::p_load(dplyr, sf, raster, devtools, ggplot2, maps, mapdata)
select <- dplyr::select

# NMFS areas original
nmfs <- st_read("Data/NMFS management area shapefiles/gf95_nmfs.shp")
nmfs <- nmfs %>% filter(NMFS_AREA %in% c(610,620,630,640,650)) %>% # subset to 610-650
  filter(GF95_NMFS1 %in% c(186,194,259,585,870)) # Removes inter-coastal in SE AK

# NMFS areas clipped to isobaths
nmfs_clipped <- st_read("Data/Depth trimmed NMFS shapefiles/NMFS610-650.shp")

# coast
coast <- maps::map(database = 'worldHires', regions = c('USA','Canada'), plot = FALSE, fill=TRUE)

coast_sf <- coast %>% 
  st_as_sf(crs = 4326) %>% 
  st_transform(crs = st_crs(nmfs_clipped)) %>% 
  st_combine() %>%
  st_crop(nmfs_clipped %>% st_bbox())

# view
ggplot()+
  geom_sf(data = nmfs, color = 'red', fill = NA)+
  geom_sf(data = nmfs_clipped, color = 'blue', fill = NA)+
  geom_sf(data = coast_sf, fill = 'grey')
  
