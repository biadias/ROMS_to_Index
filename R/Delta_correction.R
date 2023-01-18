#' Function to run delta correction on data
#'
#' @param hindcast 
#' @param historical 
#' @param projection 
#' @param ref_yrs 
#' @param lognormal 
#' @param smooth 
#'
#' @return
#' @export
#'
#' @examples
delta_correction <- function(
    hindcast = hnd,
    historical = hist,
    projection = fut,
    ref_yrs = 2000:2014,
    lognormal = FALSE,
    smooth = FALSE){
  
  # Convert to log
  if(lognormal){
    projection$value = log(projection$value)
    historical$value = log(historical$value)
    hindcast$value = log(hindcast$value)
  }
  
  # Get mean and var for each month/ROMS across years
  # - historical
  goa_clim_hist <- historical %>% 
    dplyr::filter(year %in% ref_yrs) %>%
    dplyr::group_by(month, varname, depthclass, NMFS_AREA) %>%
    dplyr::summarize(mean_hist = mean( value,na.rm=T),
                     sd_hist = sd( value, na.rm = T))
  
  # - hindcast
  goa_clim_hind <- hindcast %>% 
    dplyr::filter(year %in% ref_yrs) %>%
    dplyr::group_by(month, varname, depthclass, NMFS_AREA) %>%
    dplyr::summarize(mean_hind = mean( value,na.rm=T),
                     sd_hind = sd( value, na.rm = T))
  
  
  # calculate index using delta correction for projection
  # - Merge projection with hist and hind summary statistics 
  projection <- Reduce(function(x, y) merge(x, y, all=TRUE, by = c("month", "varname", "depthclass", "NMFS_AREA")), 
                           list(projection, goa_clim_hist, goa_clim_hind)) %>%
    arrange(varname, year, month)
  
  # - Calculate index
  if(!lognormal){
  projection <- projection %>%
    group_by(month, varname, depthclass, NMFS_AREA) %>%
    mutate(value_dc = mean_hind + (sd_hind/sd_hist * (value - mean_hist))) %>%
    ungroup() %>% arrange(year, month) %>%
    select(NMFS_AREA, depthclass, varname, year, month, date, value, value_dc, unit)
  }
  
  
  if(lognormal){
    projection <- projection %>%
      group_by(month, varname, depthclass, NMFS_AREA) %>%
      mutate(value_dc =  exp(mean_hind + (sd_hind/sd_hist * (log(value) - mean_hist)))) %>%
      ungroup() %>% arrange(year, month) %>%
      select(NMFS_AREA, depthclass, varname, year, month, date, value, value_dc, unit)
    
    # Back transform historical and projection
    historical$value = exp(historical$value)
    projection$value = exp(projection$value)
  }
  
  # - Combine with historical
  historical <- historical %>%
    mutate(value_dc = value) %>% 
    arrange(year,month) %>%
    select(NMFS_AREA, depthclass, varname, year, month, date, value, value_dc, unit)
  
  projection <- rbind(projection, historical) %>% 
    arrange(year,month)
  return(projection)
}