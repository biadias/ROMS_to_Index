# Purpose
# This loads the GOA mask and maps the ROMS coordinates to the NMFS coordinates
roms_vars <- tidync(romsfile_vars)
roms_grid <- tidync(romsfile_grid)


# Get variables. We do not need water velocity, so we can ignore u and v points. 
# Grid info
grid_variables <- hyper_grids(roms_grid) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms_grid %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })


# ROMS variables
roms_variables <- hyper_grids(roms_vars) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms_vars %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })
roms_variables$name


# Find appropriate ROMS ncdf grid for the rho points
latlon_rhogrd <- grid_variables %>% filter(name=="lat_rho") %>% pluck('grd')
# pull the lon/lats
roms_rho <- roms_grid %>% activate(latlon_rhogrd) %>% hyper_tibble() %>% dplyr::select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  mutate(rhoidx=row_number()) # add index


# Add NMFS coordinates to ROMS data.
append_xy_coords <- function(lonlatdat, xyproj=st_crs(mask), lon_col="lon_rho", lat_col="lat_rho"){
  lonlatdat %>% 
    st_as_sf(coords=c(lon_col, lat_col), crs=4326, remove=F) %>%  # convert to spatial object
    st_transform(xyproj) %>%  # convert to NMFS coords
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) # grab x and y coordinates and add them as attributes
}

rhoxy<- append_xy_coords(roms_rho,lon_col="lon_rho",lat_col="lat_rho") %>% mutate(rhoidx=row_number())

# Match ROMS and NMFS grid
# Join rho points to the NMFS mask. Do a spatial join.
rho_join <- mask %>% st_join(rhoxy) %>% na.omit()

# Get range of ROMS coordinates that correspond with GOA geometry, to subset ROMS files and reduce memory chokes
min_xi_rho <- min(rho_join$xi_rho, na.rm = TRUE)
max_xi_rho <- max(rho_join$xi_rho, na.rm = TRUE)
min_eta_rho <- min(rho_join$eta_rho, na.rm = TRUE)
max_eta_rho <- max(rho_join$eta_rho, na.rm = TRUE)


# Set up depths
# Using a custom version of Mike Sumner's angstroms::romshcoords(), because GOA ROMS have grid information in a separate file.
ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}

# Function to interporlate over depth
romshcoords_goa <- function(x, y, grid_type = "rho", slice, ..., S = "Cs_r", depth = "h", simple = FALSE){
  h <- romsdata(x, varname = depth)
  Cs_r <- ncget(y, S)
  v <- values(h)
  if (simple) {
    ## simplistic, early version - probably should be defunct
    out <- set_indextent(brick(array(rep(rev(Cs_r), each = length(v)) * v, 
                                     c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE))
  } else {
    grid_type <- match.arg(tolower(grid_type),c("rho","psi","u","v","w"))
    
    Vtransform <- as.integer(ncget(y,"Vtransform"))
    if (!Vtransform %in% c(1,2)) stop("Vtransform must be 1 or 2")
    
    hc <- ncget(y,"hc")
    
    depth_grid <- if (grid_type=="w") "w" else "rho"
    
    zeta <- if (missing(slice)) 0 else stop("slice not supported yet")##angstroms::romsdata2d(x,"zeta",slice=slice,transpose=FALSE)
    N <- length(ncget(y,"Cs_r"))
    Np <- N+1
    
    h <- ncget(x,"h")
    hmin <- min(h)
    hmax <- max(h)
    
    Lp <- dim(h)[1]
    Mp <- dim(h)[2]
    L <- Lp-1
    M <- Mp-1
    
    z <- array(NA,dim=c(Lp,Mp,if (grid_type=="w") Np else N))
    
    ## Compute vertical stretching function, C(k):
    ##stretch <- stretching(x,depth_grid)
    if (depth_grid=="w") {
      stretch <- list(C=ncget(y,"Cs_w"),s=ncget(y,"s_w"))
    } else {
      stretch <- list(C=ncget(y,"Cs_r"),s=ncget(y,"s_rho"))
    }
    
    ## Average bathymetry and free-surface at requested C-grid type.
    if (grid_type=="rho") {
      hr <- h
      zetar <- zeta
    } else if (grid_type=="psi") {
      hp <- 0.25*(h[1:L,1:M]+h[2:Lp,1:M]+h[1:L,2:Mp]+h[2:Lp,2:Mp])
      zetap <- 0.25*(zeta[1:L,1:M]+zeta[2:Lp,1:M]+zeta[1:L,2:Mp]+zeta[2:Lp,2:Mp])
    } else if (grid_type=="u") {
      hu <- 0.5*(h[1:L,1:Mp]+h[2:Lp,1:Mp])
      zetau <- 0.5*(zeta[1:L,1:Mp]+zeta[2:Lp,1:Mp])
    } else if (grid_type=="v") {
      hv <- 0.5*(h[1:Lp,1:M]+h[1:Lp,2:Mp])
      zetav <- 0.5*(zeta[1:Lp,1:M]+zeta[1:Lp,2:Mp])
    } else if (grid_type=="w") {
      hr <- h
      zetar <- zeta
    } else {
      stop("unsupported grid_type: ",grid_type)
    }
    
    ## Compute depths (m) at requested C-grid location.
    
    if (Vtransform == 1) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hp
          z[,,k] <- z0 + zetap*(1.0 + z0/hp)
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hu
          z[,,k] <- z0 + zetau*(1.0 + z0/hu)
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hv;
          z[,,k] <- z0 + zetav*(1.0 + z0/hv)
        }
      } else if (grid_type=="w") {
        z[,,1] <- -hr
        for (k in seq(from=2,to=Np,by=1)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else if (Vtransform == 2) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zeta+hr)*z0
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hp)/(hc+hp)
          z[,,k] <- zetap+(zetap+hp)*z0
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hu)/(hc+hu)
          z[,,k] <- zetau+(zetau+hu)*z0
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hv)/(hc+hv)
          z[,,k] <- zetav+(zetav+hv)*z0
        }
      } else if (grid_type=="w") {
        for (k in seq_len(Np)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zetar+hr)*z0
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else {
      stop("Vtransform must be 1 or 2")
    }
    ## FIXME all these flips and twirls can be applied more efficiently (or avoided)
    ## though should layers start at the surface and go down or ...
    
    out <- raster::flip(set_indextent(raster::brick(z, transpose = TRUE)), "y")
    ## NO - we want to start at the bottom, so we match romsdata3d
    #out <- raster::subset(out, rev(seq_len(raster::nlayers(out))))
    
  } 
  
  out
}

# Apply to get depths at $\rho$ points. For the purpose of the goa model this will only be done once, 
# Altough depth in ROMS is dynamic over time. We only do it once in Atlantis too - 
# for models like these the bending free surface amounts to rounding error.

# convert ROMS s-coordinates to depth with Mike Sumner's angstroms package
romsdepths <- romshcoords_goa(x = romsfile_grid, y = romsfile_vars, S = "Cs_r", depth = "h")

# using tabularaster to convert to tibble
# and a indexing template with "by_column" filling
romsi <- crossing(xi_rho=1:dim(romsdepths)[2],eta_rho=1:dim(romsdepths)[1]) %>% arrange(-eta_rho) %>% mutate(cellindex=row_number()) # making sure that the join by cellindex below is correct - doing this for consistency with the way tabularaster::as_tibble() unpacks the raster cells 

romsdepthsdf <- tabularaster::as_tibble(romsdepths,dim=F) %>% 
  arrange(cellindex) %>% 
  left_join(romsi,by='cellindex') %>% 
  set_names(c("romsdepth","cellindex","xi_rho","eta_rho")) %>% 
  group_by(cellindex,xi_rho,eta_rho) %>% 
  nest(romsdepth=c(romsdepth)) %>% ungroup() %>% 
  mutate(romsdepth=purrr::map(romsdepth,function(x)x[['romsdepth']])) %>%
  filter(between(xi_rho, min_xi_rho, max_xi_rho) & between(eta_rho, min_eta_rho, max_eta_rho))
