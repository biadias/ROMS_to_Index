library(tidyverse)

old_goa_temp_610_to_630_FebApril <- read.csv('Output/Old_goa_temp_610_to_630_FebApril_300M.csv')
revised_goa_temp_610_to_630_FebApril <- read.csv('Output/goa_temp_610_to_630_FebApril_300M.csv')

old_goa_temp_610_to_630_FebApril <- old_goa_temp_610_to_630_FebApril %>%
  mutate(revised = FALSE)

revised_goa_temp_610_to_630_FebApril <- revised_goa_temp_610_to_630_FebApril %>%
  mutate(revised = TRUE)


goa_temp_610_to_630_FebApril <- rbind(old_goa_temp_610_to_630_FebApril, revised_goa_temp_610_to_630_FebApril)

# Plot

# -- Plot
# - SST
ggplot(goa_temp_610_to_630_FebApril %>% filter(depthclass == "Surface" & hind == "yes"), aes(year, mean_value_dc_610_to_630, colour = revised)) + 
  geom_line() +
  ylab("SST (Celsius)") + 
  xlab("Year") + 
  facet_wrap(~simulation)

# - Bottom temp
ggplot(goa_temp_610_to_630_FebApril %>% filter(depthclass == "Bottom" & hind == "yes"), aes(year, mean_value_dc_610_to_630, colour = revised)) + 
  geom_line() +
  ylab("Bottom temp (Celsius)") + 
  xlab("Year") + 
  facet_wrap(~simulation)
