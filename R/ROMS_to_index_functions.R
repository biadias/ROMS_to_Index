# Purpose
# Function to interpolate a ROMS variable across depth and take summary statistics across depth and area


#' Interpolate variables with cubic splines over 1 m intervals in the water column.
#'
#' @param romsdepths depths of the ROMS layers
#' @param romsvar variable to interpolate
#'
#' @return
#' @export
interp_foo <- function(romsdepths,romsvar) {
  depths_out <- seq(round(min(romsdepths)),0,by=1) # 1m interpolation, starting from deepest
  interp <- spline(romsdepths,romsvar,xout=depths_out) %>% pluck('y')
  return(tibble(depth=depths_out,val=interp))
}


#' Interpolates ROMS variable across depths
#'
#' @param romsfile path to the ROMS NetCDF file
#' @param variable variable name in ROMS netcdf file (e.g. temp, salt, etc)
#' @param time_step oceanographic time step (e.g. 3377894400)
#' @param this_roms_vars tidync of ROMS netcdf
#' @param this_roms_variables data.frame of variables with associated ROMS dim
#' @param max_depth Min depth of the spatial area (e.g. 0 m)
#' @param max_depth Max depth of the spatial area (e.g. 1000 m)
#' @param average TRUE/FALSE. True = take average over depth/spatial range. False = takes sum over depth/spatial range
#'
#' @description 
# 1. Applies the cubic spline interpolation at each rho point.
# 2. Integrates over the water column as appropriate for each variable (i.e. average for state variables like temperature and salinity, sum for concentrations per m3 to obtain values per m2).
#' @return Return average values for each spatial-area in the GOA model domain.
#' @export
#' 
#' 
interpolate_var <- function(romsfile, variable, time_step, this_roms_vars, this_roms_variables, min_depth = 0, max_depth = -1000, average = TRUE){
  
  print(paste("doing", variable, "for", romsfile, sep = " "))
  
  #get unit
  this_unit <- ncmeta::nc_atts(romsfile, variable) %>% filter(name == 'units') %>% tidyr::unnest(cols = c(value)) %>% pull(value)
  if(length(this_unit)==0){this_unit <- NA} # or else salt will break it
  
  grd <- this_roms_variables %>% filter(name==variable) %>% pluck('grd')
  
  dat <- this_roms_vars %>% activate(grd) %>%
    hyper_tibble(select_var=variable, 
                 xi_rho = between(xi_rho, min_xi_rho, max_xi_rho), 
                 eta_rho = between(eta_rho, min_eta_rho, max_eta_rho),
                 ocean_time = ocean_time == time_step)
  
  # filter based on xi and eta of the points inside NMFS areas and with h <= 1000
  dat <- dat %>%
    mutate(xi_eta = paste(xi_rho, eta_rho, sep = '_')) %>%
    filter(xi_eta %in% xi_eta_set) %>%
    select(-xi_eta)
  
  interp_dat <- dat %>% 
    dplyr::select(xi_rho,eta_rho,!!variable) %>% 
    nest(data=c(!!variable))%>% 
    mutate(evar=purrr::map(data,~.[[1]]))
  
  interp_dat <- interp_dat %>%
    inner_join(romsdepthsdf,by=c('xi_rho','eta_rho')) 
  
  # Interpolate variable across depth layers (ROMS layers to 1 m layers)
  interp_dat <- interp_dat %>% 
    mutate(interp = purrr::map2(romsdepth,evar,interp_foo)) %>% 
    dplyr::select(-data,-evar,-romsdepth)
  
  # join to rho_join set to subset to goa mask only
  interp_dat <- rho_join %>%
    st_set_geometry(NULL) %>%
    left_join(interp_dat,by=c('xi_rho','eta_rho'))
  
  # Unnest for averaging
  interp_dat <- interp_dat %>%
    unnest(interp) %>%
    lazy_dt()
  
  #TODO: spit a warning if for any rho point temp at the surface is lower than at depth - may be sign of depths from ROMS being inverted
  
  # # Double check how many ROMS cells are deeper than 1000 after masking
  # check <- interp_dat %>%
  #   group_by(cellindex) %>%
  #   summarise(maxdepth = max(abs(depth))) %>% as.data.frame()
  # sum(check$maxdepth > abs(max_depth))/nrow(check)
  # hist(check$maxdepth)
  
  # # Double check how many ROMS cells are deeper than 1000 after masking
  # check <- interp_dat %>%
  #   group_by(cellindex) %>%
  #   summarise(maxdepth = max(abs(depth))) %>% 
  #   as.data.frame() %>%
  #   pull(maxdepth)
  #   
  # if(length(check[check>abs(max_depth) | check<abs(min_depth)]) > 0) stop("Some depths of rho points are outside the accepted range")

  # Define depth class and remove areas with depth over max_depth
  interp_dat <- interp_dat %>%
    group_by(cellindex) %>%
    mutate(maxdepth = max(abs(depth)),
           depthclass = case_when(
             abs(depth) <= 10 ~ "Surface", # this means that for very shallow points (h = 10) we only get surface variables
             abs(depth) >= (maxdepth-10) ~ "Bottom",
             abs(depth) > 10 & abs(depth) < (maxdepth-10) ~ "Midwater"
           )) %>%
    ungroup() %>%
    filter(maxdepth <= abs(max_depth) & maxdepth >= abs(min_depth)) # now just a sanity check but should no longer be needed
  
  # Calculate summary statistics over depth and spatial range
  if(average){ # For state variables, take the average over the water column. 
    
    # Summarize by spatial area and depthclass
    goa_dat_area_depth <- interp_dat %>%
      group_by(NMFS_AREA, depthclass) %>%
      summarise(value = mean(val, na.rm = TRUE)) %>% 
      as.data.frame()
    
    # Summarize across spatial areas
    goa_dat_depth <- interp_dat %>%
      group_by(depthclass) %>%
      summarise(value = mean(val, na.rm = TRUE)) %>% 
      mutate(NMFS_AREA = "All") %>%
      as.data.frame()
    
    # Summarize for each spatial area across depthclasses
    goa_dat_area <- interp_dat %>%
      group_by(NMFS_AREA) %>%
      summarise(value = mean(val, na.rm = TRUE)) %>% 
      mutate(depthclass = "All") %>%
      as.data.frame()
    
    # Summarize across spatial areas and depthclasses
    goa_dat <- interp_dat %>%
      summarise(value = mean(val, na.rm = TRUE)) %>% 
      mutate(NMFS_AREA = "All", depthclass = "All") %>%
      as.data.frame()
    
    goa_dat <- do.call("rbind", list(goa_dat_area_depth, goa_dat_depth, goa_dat_area, goa_dat))
    
    goa_dat$summaryStat = "mean"
  } else { # For concentrations, sum over the water column.
    
    # Summarize by spatial area and depthclass
    goa_dat_area_depth <- interp_dat %>%
      group_by(cellindex, depthclass, NMFS_AREA) %>%
      summarise(value_m2 = sum(val, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(NMFS_AREA, depthclass)  %>% 
      summarise(value=mean(value_m2,na.rm=TRUE)) %>% 
      as.data.frame()
    
    # Summarize across spatial areas
    goa_dat_depth <- interp_dat %>%
      group_by(cellindex, depthclass) %>%
      summarise(value_m2 = sum(val, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(depthclass)  %>% 
      summarise(value=mean(value_m2,na.rm=TRUE)) %>% 
      mutate(NMFS_AREA = "All") %>%
      as.data.frame()
    
    # Summarize for each spatial area across depthclasses
    goa_dat_area <- interp_dat %>%
      group_by(cellindex, NMFS_AREA) %>%
      summarise(value_m2 = sum(val, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(NMFS_AREA)  %>% 
      summarise(value=mean(value_m2,na.rm=TRUE)) %>% 
      mutate(depthclass = "All") %>%
      as.data.frame()
    
    # Summarize across spatial areas and depthclasses
    goa_dat <- interp_dat %>%
      group_by(cellindex) %>%
      summarise(value_m2 = sum(val, na.rm = TRUE)) %>%
      ungroup() %>%
      summarise(value=mean(value_m2,na.rm=TRUE)) %>% 
      mutate(NMFS_AREA = "All", depthclass = "All") %>%
      as.data.frame()
    
    goa_dat <- do.call("rbind", list(goa_dat_area_depth, goa_dat_depth, goa_dat_area, goa_dat))
    goa_dat$summaryStat = "mean (per area) sum across depthclass"
  }
  
  goa_dat$varname = variable
  goa_dat$time_step = time_step
  goa_dat$unit = this_unit
  return(goa_dat)
}


#' Derive index per one netcdf file
#'
#' @param this_romsfile 
#' @param average True or False. Wether to take the average across depth x area (TRUE) or take the sum across depths and average those sums across areas (FALSE).
#' @param variables The variables to extract from the netCDF files
#' @param max_depth Min depth of the spatial area (e.g. 0 m)
#' @param max_depth Max depth of the spatial area (e.g. 1000 m)
#'
#' @description Given a netcdf file, the function applies the interpolation function in \code{interpolate_var} to all variables in all time steps within the netcdf
#' @return
#' @export
#'
summarize_netcdf <- function(this_romsfile, variables = c("temp"), average = TRUE, min_depth = 0, max_depth = -1000){
  # read roms netcdf file
  
  this_roms_vars <- tidync(this_romsfile)
  this_roms_variables <- hyper_grids(this_roms_vars) %>% # all available grids in the ROMS ncdf
    pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
    purrr::map_df(function(x){
      this_roms_vars %>% activate(x) %>% hyper_vars() %>% 
        mutate(grd=x)
    })
  
  # get time step
  time_grd <- this_roms_variables %>% filter(name=="ocean_time") %>% pluck('grd')
  roms_time <- this_roms_vars %>% activate(time_grd) %>% hyper_tibble() %>% pull()
  
  # get variable/time step combinations
  var_time_combos <- expand.grid(variables,roms_time) %>% mutate(Var1=as.character(Var1)) %>%
    set_names(c('variable','time_step'))
  
  # read variables and carry out interplation
  goa_vals <- apply(var_time_combos, 1, function(x) interpolate_var(romsfile = this_romsfile,
                                                                    variable = as.character(x[1]), 
                                                                    time_step = roms_time,
                                                                    this_roms_vars, this_roms_variables,
                                                                    min_depth = min_depth, max_depth = max_depth, 
                                                                    average = average))
  goa_vals <- do.call("rbind", goa_vals)
  
  # turn ocean_time to a date - make sure you find the correct Epoch for your model
  epoch <- "1900-01-01" #important, check that this is your correct start for ocean_time or the dates will be messed up
  
  goa_vals <- goa_vals %>% mutate(date=as.POSIXct(time_step, origin = epoch, tz='UTC')) %>%
    arrange(varname,time_step)
  
  return(goa_vals)
}



#' Function to derive summary statistics for multiple netcdfs
#'
#' @param netcdf_files vector of file names and directories for netcdfs
#' @param variables variable names for extraction
#' @param max_depth Min depth of the spatial area (e.g. 0 m). Points with max depth shallower than this in the ROMS data will be removed.
#' @param max_depth Max depth of the spatial area (e.g. 1000 m). Points with max depth deeper than this in the ROMS data will be removed.
#' @param average True or False. Wether to take the average across depth x area (TRUE) or take the sum across depths and average those sums across areas (FALSE).
#' 
#' @description wrapper function for lapply
#'
#' @return a data.frame with index values
#' @export
roms_to_goa <- function(netcdf_files, variables = c("temp"), min_depth = 0, max_depth = -1000, average = TRUE){
  
  vals_list <- lapply(netcdf_files, function(x) summarize_netcdf(this_romsfile = x, variables, min_depth = min_depth, max_depth = max_depth, average))
  return(do.call("rbind", vals_list))
}
