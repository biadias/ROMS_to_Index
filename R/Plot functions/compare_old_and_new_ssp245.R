# ALberto Rovellini
# 6/20/2024
# Script to compare old NEP 10K ssp245 indices (before May 2024 revision from Al Hermann) to indices after revision.
# the main interest here is in temperature and salinity

library(tidyverse)
library(lubridate)

# work on avg files

# read files in

# 300 m
old_ssp245_300 <- read.csv('Data/NEP_10k_revised_indices/nep_avg_wb_ssp245_300_old.csv', header = T)
new_ssp245_300 <- read.csv('Data/NEP_10k_revised_indices/nep_avg_wb_ssp245_300.csv', header = T)

# start from 300 m
old_ssp245_300 <- old_ssp245_300 %>%
  mutate(run = 'old', year = year(date), month = month(date))

new_ssp245_300 <- new_ssp245_300 %>%
  mutate(run = 'new', year = year(date), month = month(date))

ssp245_300 <- rbind(old_ssp245_300, new_ssp245_300) 

# make some plots
p1 <- ssp245_300 %>%
  filter(varname %in% c('temp','salt'), NMFS_AREA == 'All') %>%
  ggplot()+
  geom_line(aes(x = as.Date(date), y = value, color = run), linewidth = 0.85)+
  theme_bw()+
  labs(title = 'Time series by depth from ROMS ssp245 for temp and salt')+
  facet_grid(varname~depthclass, scales = 'free_y')
p1

# see difference
old_ssp245_300_thin <- old_ssp245_300 %>%
  filter(varname %in% c('temp','salt'), NMFS_AREA == 'All') %>%
  dplyr::select(year, month, varname, depthclass, value) %>%
  rename(old_val = value)

new_ssp245_300_thin <- new_ssp245_300 %>%
  filter(varname %in% c('temp','salt'), NMFS_AREA == 'All') %>%
  dplyr::select(year, month, varname, depthclass, value) %>%
  rename(new_val = value)

ssp245_300_thin <- old_ssp245_300_thin %>%
  left_join(new_ssp245_300_thin) %>%
  mutate(delta = new_val - old_val,
         date = ymd(paste(year, month, '15', sep = '-')))

p2 <- ssp245_300_thin %>%
  ggplot()+
  geom_line(aes(x = date, y = delta), linewidth = 1)+
  geom_hline(yintercept = 0, linewidth = 1.2, color = 'red')+
  theme_bw()+
  labs(title = 'Delta between ROMS hindcast for temp and salt (new - old value) by depth')+
  facet_grid(varname~depthclass, scales = 'free_y')
p2
