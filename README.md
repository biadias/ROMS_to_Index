# Overview

This repo overviews the process for making environmental indices from the ROMS NEP data including: data overview, index creation, and delta correction. We don't proved generic functions because each GOA-CLIM modeler may have indices specific to a particular spatial and temporal domain. We therefore provide an example Quarto that folks can copy and derive their desired indices.

# Data overview

Projections of environmental variables biomass are derived from a downscaled Intergovernmental Panel on Climate Change (IPCC) projection model from the coupled model intercomparison project (CMIP): GFDL-ESM2M. The IPCC projection model was downscaled using the North East Pacific (NEP) Regional Ocean Modeling System (ROMS) developed for the Gulf of Alaska.

ROMS simulations are referred to as the following:

- hindcast: representing model spinup (2000 to 2020),

- projection: GFDL-ESM2M downscaled projection (2015 to 2099)

- historical: representing the final year of the spinup forced with observed oceanographic conditions to better represent historical conditions (1980 to 2014).

Variables from ROMS simulations averaged for each layer and grid cell on a monthly time scale were averaged or summed across the following vertical distributions:

- Surface (from 0 to 10 m depth)

- Bottom (from bottom to 10 m from bottom)

- Midwater (from 10 m depth to 10 m from bottom)

The averaged for each NMFS management area to the 1,000 m isobath (areas are calculated from st_area in the sf package):

- 610 (area = 63986698621 $m^2$)

- 620 (area = 69583703140 $m^2$)

- 630 (area = 105918077937 $m^2$)

- 640 (area = 37270389681 $m^2$)

- 650 (area = 43952466109 $m^2$)

- All areas

The "raw" monthly data are located in the "/Data/NEP_10k_indices" folder. They represent the monthly indices either averaged across each depth and spatial domain ("nep_avg_…") or summed across each depth domain and then averaged across the spatial domain ("nep_sum…").


# Delta correction overview

Delta correction is done to correct for shifts in the time series between the hindcast and projection. For normally distributed variables (*Y*):

$$
X^{proj'}_t=\bar{X}^{hind}_{\bar{T}}+
\left(\frac{\sigma^{hind}_{\bar{T}}}{\sigma^{hist}_{\bar{T}}} * 
\left(X^{proj}_t-\bar{X}_{\bar{T}}^{hist} \right)\right)
$$

where $X^{proj'}_t$ is the bias corrected projection index in time-step $t$, $\bar{X}^{hind}_{\bar{T}}$ is the average index value from the hindcast during reference period $T$, $\sigma^{hind}_{\bar{T}}$ is the standard deviation of the index from the hindcast during the reference period, $\sigma^{hist}_{\bar{T}}$ is the standard deviation of the index from the historical run during the reference period, $X^{proj}_t$ is the non bias-corrected projection index, $\bar{X}_{\bar{T}}^{hist}$ is the average value from historical run during the reference period.

For log-normally distributed variables the formula can be adjusted as follows:

$$
X^{proj'}_t=\exp\left(\bar{logX}^{hind}_{\bar{T}}+
\left(\frac{\sigma^{hind}_{\bar{T}}}{\sigma^{hist}_{\bar{T}}} * 
\left(logX^{proj}_t-\bar{logX}_{\bar{T}}^{hist} \right)\right)\right)
$$

where $T$, $\sigma^{hind}_{\bar{T}}$ is the standard deviation of the log index from the hindcast during the reference period, $\sigma^{hist}_{\bar{T}}$ is the standard deviation of the log indext from the histroical run during the reference period.

Because the smallest time step in month, the mean and variance terms should be calculated as follows:

Monthly indices: time-step $t$ should be year $y$ and month $m$ (e.g. $t=y,m$). A function to do the bias correction is in the [R/Delta_correction.r](https://github.com/GOA-CLIM/ROMS_to_Index/blob/main/R/Delta_correction.R)


# Index derivation

Environmental indices can be derived for each modelers following their specific needs. We don't develop all indices here as individual modeler's needs may vary. But, indices integrated across each spatial NMFS management area should be area weighted using the following areas.

A quarto document document exampling some index generation can be found [here](https://github.com/GOA-CLIM/ROMS_to_Index/blob/main/ROMS%20index%20generation.qmd)


# Locations of monthly average output files on loon

*File naming convention:*
The suffix "moave_{yr}_{mo}"  indicates monthly average, year, and month, e.g. /ahr0/hermann/goa-output/monthly_aves/cgoa_moave_2000_01.nc is the monthly average for January of year 2000


*The "master directories" and grid files are here:*

/ahr0/hermann/goa-output/


*The latitude, longitude, and depth locations of each grid point are in these files:*

- zpoints_cgoa.nc

- zpoints_nep.nc


*For variables like temperature and salinity, the relevant variables for latitude, longitude and depth are*
lat_rho, lon_rho, zval_rho

**Note** there are values for u, v, and w variables, which are on slightly different (staggered) grids.


*The major subdirectories of /ahr0/hermann/goa-output/ are as follows:*

- cgoa 3km model hindcast: monthly_aves

- nep 10km model hindcast: monthly_aves_nep_hind

- nep 10km model GFDL "historical" run: monthly_aves_nep_wb_hist

- nep 10km model GFDL ssp126 projection: monthly_aves_nep_wb_ssp126

- nep 10km model GFDL ssp585 projection: monthly_aves_nep_wb_ssp585

# Aclim links
Can be used to develop indices from 10K resolution ROMS models given an input of a spatial field following:
https://kholsman.github.io/ACLIM2/#23_ROMSNPZ_variables
