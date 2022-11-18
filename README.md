# ROMS_to_Index
Code forked from Alberto Rovellini for developing environmental indices from Regional Ocean Modelling System models (ROMS)
Can be used to develop indices from 10K resolution ROMS models given an input of a spatial field following:
https://kholsman.github.io/ACLIM2/#23_ROMSNPZ_variables


**Locations of monthly average output files on loon**

*File naming convention:*
The suffix "moave_{yr}_{mo}"  indicates monthly average, year, and month, e.g. /ahr0/hermann/goa-output/monthly_aves/cgoa_moave_2000_01.nc is the monthly average for January of year 2000


*The "master directories" and grid files are here:*

/ahr0/hermann/goa-output/


*The latitude, longitude, and depth locations of each grid point are in these files:*

zpoints_cgoa.nc

zpoints_nep.nc


*For variables like temperature and salinity, the relevant variables for latitude, longitude and depth are*
lat_rho, lon_rho, zval_rho

**Note** there are values for u, v, and w variables, which are on slightly different (staggered) grids.


*The major subdirectories of /ahr0/hermann/goa-output/ are as follows:*

cgoa 3km model hindcast:
monthly_aves

nep 10km model hindcast:
monthly_aves_nep_hind

nep 10km model GFDL "historical" run:
monthly_aves_nep_wb_hist

nep 10km model GFDL ssp126 projection:
monthly_aves_nep_wb_ssp126

nep 10km model GFDL ssp585 projection:
monthly_aves_nep_wb_ssp585
