# ALberto Rovellini
# 6/16/2023
# Script to compare old NEP 10K hindcast indices (before June 2023 revision from Al Hermann) to indices after revision.
# the main interest here is in temperature and salinity

library(tidyverse)
library(lubridate)

# work on avg files

# read files in

# 300 m
old_hind_300 <- read.csv('Data/NEP_10k_indices/nep_avg_hind_300.csv', header = T)
new_hind_300 <- read.csv('Data/NEP_10k_revised_indices/nep_avg_revised_hind_300.csv', header = T)

# start from 300 m
old_hind_300 <- old_hind_300 %>%
  mutate(run = 'old', year = year(date), month = month(date))

new_hind_300 <- new_hind_300 %>%
  mutate(run = 'new', year = year(date), month = month(date))

hind_300 <- rbind(old_hind_300, new_hind_300) %>%
  filter(year >= 2000)

# make some plots
p1 <- hind_300 %>%
  filter(varname %in% c('MZS','MZL'), NMFS_AREA == 'All') %>%
  ggplot()+
  geom_line(aes(x = as.Date(date), y = value, color = run), linewidth = 0.85)+
  theme_bw()+
  labs(title = 'Time series by depth from ROMS hindcast for temp and salt')+
  facet_grid(varname~depthclass, scales = 'free_y')
p1
#ggsave('hindcast_ts_salt_temp.png', p1, width = 12, height = 9)

# see difference
old_hind_300_thin <- old_hind_300 %>%
  filter(varname %in% c('temp','salt'), NMFS_AREA == 'All') %>%
  dplyr::select(year, month, varname, depthclass, value) %>%
  rename(old_val = value)

new_hind_300_thin <- new_hind_300 %>%
  filter(varname %in% c('temp','salt'), NMFS_AREA == 'All') %>%
  dplyr::select(year, month, varname, depthclass, value) %>%
  rename(new_val = value)

hind_300_thin <- old_hind_300_thin %>%
  left_join(new_hind_300_thin) %>%
  mutate(delta = new_val - old_val,
         date = ymd(paste(year, month, '15', sep = '-')))

p2 <- hind_300_thin %>%
  ggplot()+
  geom_line(aes(x = date, y = delta), linewidth = 1)+
  geom_hline(yintercept = 0, linewidth = 1.2, color = 'red')+
  theme_bw()+
  labs(title = 'Delta between ROMS hindcast for temp and salt (new - old value) by depth')+
  facet_grid(varname~depthclass, scales = 'free_y')
p2
#ggsave('hindcast_delta_salt_temp.png', p2, width = 12, height = 9)

