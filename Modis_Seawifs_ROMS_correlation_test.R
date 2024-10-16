
library(tidyverse)

roms <- read.csv("Output/CGOA_NEP_ROMZ_PP_monthly_630_640_ssp126_300_surface_mt_km2.csv")
roms_phL <- read.csv("CGOA_NEP_ROMZ_PhLBiom_monthly_630_640_ssp126_300_mt_km2.csv")
roms_phS <- read.csv("CGOA_NEP_ROMZ_PhSBiom_monthly_630_640_ssp126_300_mt_km2.csv")
sat <- read.csv("Data/modis_sewifs_chloa.csv")

combined_data <- merge(roms, sat, by = c("year", "month")) %>% 
  select(c("year", "month", "CGOA_prod_Ph", "Mean_month_chloa"))



plot(combined_data$Mean_month_chloa, combined_data$CGOA_prod_Ph, main = "Chlorophyll a vs Primary Production",
     xlab = "Chlorophyll a", ylab = "Primary Production", pch = 19)

qqnorm(combined_data$Mean_month_chloa); qqline(chl_a, col = "red")
qqnorm(combined_data$CGOA_prod_Ph); qqline(pp, col = "red")

cor_coeff_spearman <- cor(combined_data$Mean_month_chloa, combined_data$CGOA_prod_Ph, method = "spearman")


cor_test_spearman <- cor.test(combined_data$Mean_month_chloa, combined_data$CGOA_prod_Ph, method = "spearman")
print(cor_test_spearman)

#There is a significant correlation between Chlorophyll a and primary production.

plot(combined_data$Mean_month_chloa, combined_data$CGOA_prod_Ph,
     main = "Chlorophyll a vs Primary Production",
     xlab = "Mean Monthly Chlorophyll a",
     ylab = "Primary Production",
     pch = 19)
lines(lowess(combined_data$Mean_month_chloa, combined_data$CGOA_prod_Ph), col = "blue")


hist(combined_data$Mean_month_chloa, main = "Histogram of Chlorophyll a",
     xlab = "Mean Monthly Chlorophyll a")
hist(combined_data$CGOA_prod_Ph, main = "Histogram of Primary Production",
     xlab = "Primary Production")

# Q-Q Plots
qqnorm(combined_data$Mean_month_chloa); qqline(combined_data$Mean_month_chloa, col = "red")
qqnorm(combined_data$CGOA_prod_Ph); qqline(combined_data$CGOA_prod_Ph, col = "red")


# Log Transformation
combined_data$log_chloa <- log(combined_data$Mean_month_chloa + 1) # Adding 1 to avoid log(0)
combined_data$log_prod_Ph <- log(combined_data$CGOA_prod_Ph + 1)

library(mgcv)
gam_model <- gam(CGOA_prod_Ph ~ s(Mean_month_chloa), data = combined_data)
plot(gam_model, main = "GAM Smooth of Primary Production vs Chlorophyll a")



# Scatter plot of log-transformed data
plot(combined_data$log_chloa, combined_data$log_prod_Ph,
     main = "Log-Transformed Chlorophyll a vs Primary Production",
     xlab = "Log of Mean Monthly Chlorophyll a",
     ylab = "Log of Primary Production",
     pch = 19)

# Add regression line
abline(lm(combined_data$log_prod_Ph ~ combined_data$log_chloa), col = "blue", lwd = 2)

# Q-Q Plots for log-transformed data
qqnorm(combined_data$log_chloa); qqline(combined_data$log_chloa, col = "red")
qqnorm(combined_data$log_prod_Ph); qqline(combined_data$log_prod_Ph, col = "red")


lm_model <- lm(log_prod_Ph ~ log_chloa, data = combined_data)
summary(lm_model)


# Plot residuals
plot(lm_model$fitted.values, lm_model$residuals,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

qqnorm(lm_model$residuals)
qqline(lm_model$residuals, col = "red")

## Phyto roms and sat
df_list <- list(roms,roms_phL, roms_phS, sat)
combined_df <- df_list %>% reduce(full_join, by=c("year", "month")) %>% 
  select(c("year", "month", "CGOA_prod_Ph","PhL_biom", "PhS_biom", 
           "Mean_month_chloa", "Mean_SD"))

chl_PhL <- combined_df$PhL_biom
chl_sat <- combined_df$Mean_month_chloa

plot(chl_PhL, chl_sat, main = "Chlorophyll a vs Primary Production",
     xlab = "Chlorophyll a roms", ylab = "Chlorophyll a sat", pch = 19)

# Histograms
hist(chl_PhL, main = "Histogram of Chlorophyll a roms", xlab = "Chlorophyll a roms")
hist(chl_sat, main = "Histogram of Chlorophyll a sat", xlab = "Chlorophyll a sat")

# Q-Q Plots
qqnorm(chl_PhL); qqline(chl_PhL, col = "red")
qqnorm(chl_sat); qqline(chl_sat, col = "red")

cor_coeff_spearman <- cor(chl_PhL, chl_sat, method = "spearman")
cor_test_spearman <- cor.test(chl_PhL, chl_sat, method = "spearman")
print(cor_test_spearman)

# Adding a regression line
model <- lm(chl_PhL ~ chl_sat)
abline(model, col = "blue", lwd = 2)


#write.csv(combined_df, "Output/CGOA_ROMZpp_modis_mo_630_640_ssp126_300_surface_mt_km2.csv")

#Transformation to generate the PP forcing function centered in 1 
# Load the dataset (assuming the dataset is in a CSV file)
data <- read.csv("Output/CGOA_ROMZpp_modis_mo_630_640_ssp126_300_surface_mt_km2.csv")

# Calculate the mean of CGOA_prod_Ph
mean_CGOA_prod_Ph <- mean(data$CGOA_prod_Ph, na.rm = TRUE)

# Transform CGOA_prod_Ph to be centered at 1
mean_CGOA_prod_Ph <- 132.5985951	#1990-2015 mean -ROMS hind period


data <- data %>%
  mutate(CGOA_prod_Ph_centered_1 = ((CGOA_prod_Ph - mean_CGOA_prod_Ph) / mean_CGOA_prod_Ph) + 1)


# View the updated dataset
head(data)
write.csv(data, "Output/CGOA_ROMZpp_modis_mo_630_640_ssp126_300_surface_mt_km2_v2.csv")

#Now,lets make a time series object
ts <- ts(data$CGOA_prod_Ph_centered_1,start = c(1980,1),end = c(2099,12),frequency = 12,class = "ts")

#decomposition of ts object into trend,seasonality and error by additive model
decomposed <- decompose(ts,type = "additive")

#plot the components of an additive time series
theme_set(theme_bw())
forecast::autoplot(decomposed)

set.seed(12345)
decomposed_trend <- ts - decomposed$seasonal

decomposed_trend_df <- data.frame(date =as.Date(zoo::as.yearmon(time(decomposed_trend))),
                                  PP=as.matrix(decomposed_trend))

write.csv(decomposed_trend_df, "Output/CGOA_ROMSpp_detrended.csv")
