# *------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: .R
# | DATE: 
# | CREATED BY:  Jim Stagge         
# *----------------------------------------------------------------
# | PURPOSE:  
# | 
# |
# *------------------------------------------------------------------

###########################################################################
## Set the Paths
###########################################################################
### Path for Data and Output	
data_path <- "../../data"
output_path <- "../../output"
function_path <- "./functions"

prism_path <- file.path(data_path, "prism")

### Set up output folders
write_output_path <- file.path(output_path, "bear_runoff_model")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Load functions
###########################################################################
require(sf)
require(staggefuncs)
require(waterData)
require(lubridate)
require(airGR)
require(tabularaster)
require(raster)
require(viridis)
require(tidyverse)
require(Evapotranspiration)

select <- dplyr::select


###########################################################################
## Load custom functions
###########################################################################
### Load project specific functions
file.sources = list.files(function_path, pattern="*.R", recursive=TRUE)
sapply(file.path(function_path, file.sources),source)

###########################################################################
## Define Projection
###########################################################################
#I will use USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
#http://spatialreference.org/ref/sr-org/8711/
### This is equivalent to 
#  +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs

### Or  EPSG = 42303
usgs_proj <- 42303

###########################################################################
## Read in Watersheds
###########################################################################
watershed_location <- file.path(data_path, "usgs_hydrography/WBD_16-Utah/NHDPLUS_H_1601_HU4_GDB.gdb")

### Read in HU12 and project to USGS equal area projection
hu_12 <- st_read(watershed_location, layer = "WBDHU12") %>% 
	st_transform(usgs_proj)

head(hu_12)
plot(hu_12)

### Cut the HUC 4 resolution to Bear River watershed
hu_4 <- st_read(watershed_location, layer = "WBDHU4") %>% 
	filter(Name == "Bear") %>%
	st_transform(usgs_proj)

st_crs(hu_4)

### List of HUC regions considered here
huc_list <- c("160101010101", "160101010102", "160101010103", "160101010104", "160101010106")


###########################################################################
## Download Daily PRISM Climate data and find mean over each HUC region
###########################################################################
for (i in seq(1,length(huc_list))) {

	to_read <- file.path(prism_path, paste0("PRISM_ppt_tmin_tmean_tmax_stable_4km_", huc_list[i], ".csv") )
	prism_read <- clim_gauge <- read_csv(to_read, skip=10, col_names = TRUE)

	prism_read <- prism_read %>%
		rename(prcp_mm = `ppt (mm)`) %>%
		rename(tmin_c = `tmin (degrees C)`) %>%
		rename(tavg_c = `tmean (degrees C)`) %>%
		rename(tmax_c = `tmax (degrees C)`) %>%
		rename(date = Date) %>%
		mutate(huc = huc_list[i])
	
	if (i == 1) {
		prism_df <- prism_read
	} else {
		prism_df <- rbind(prism_df, prism_read)
	}

}

###########################################################################
## PET calculation using Oudin
###########################################################################
### Need to separate each one because the function can't handle different Latitudes

for (i in seq(1, dim(hu_char)[1])){

lat_temp <- hu_char$Lat[i]

pet_temp <- prism_df %>%
	filter(huc == hu_char$huc[i]) %>%	
	left_join(hu_char, by="huc") %>%
	select(date, huc, tavg_c, Lat) %>%
	mutate(JD = as.POSIXlt(date)$yday) %>%
	drop_na() %>%
	mutate(pet_oudin = PEdaily_Oudin(JD = JD, Temp = tavg_c,  Lat = lat_temp, LatUnit = "deg")) %>%
	select(date, huc, pet_oudin)
	
	
if (i == 1) {
	pet_df <- pet_temp
} else {
	pet_df <- rbind(pet_temp, pet_df)
}
}

ggplot(pet_df, aes(x=date, y=pet_oudin, colour=huc, group=huc)) + geom_line() +scale_colour_brewer(type="qual") + theme_classic()

### Combine back with prism
prism_df <- prism_df %>%
	left_join(pet_df, by= c("date", "huc"))

### And create a POSIX date column
Sys.setenv(TZ='UTC')
prism_df$date_pos <- as.POSIXct(prism_df$date, format="%Y-%m-%d")


###########################################################################
## Calculate PET using Hargreaves Samani
###########################################################################
### Need to separate each one because the function can't handle different Latitudes

for (i in seq(1, dim(hu_char)[1])){

data("constants")

lat_temp <- hu_char$Lat[i]
constants$lat <- lat_temp
constants$lat_rad <- lat_temp * (pi/180)
constants$Elev <- hu_char$elev_m[i]


pet_temp <- prism_df %>%
	filter(huc == hu_char$huc[i])
	
pet_clim <- list(Date.daily = pet_temp$date, 
	J = zoo(yday(pet_temp$date),pet_temp$date), 
	Tmax = zoo(pet_temp$tmax_c,pet_temp$date), 
	Tmin = zoo(pet_temp$tmin_c,pet_temp$date))

pet_results <- ET.HargreavesSamani(pet_clim, constants, ts="daily", message="no", save.csv="no")
pet_results <- tibble(date = as.Date(time(pet_results$ET.Daily)), pet_har_mm = as.double(pet_results$ET.Daily))

pet_temp <- pet_temp %>%
	left_join(pet_results, by="date") %>%
	mutate(JD = as.POSIXlt(date)$yday) %>%
	select(date, huc, pet_har_mm)
	
if (i == 1) {
	pet_df <- pet_temp
} else {
	pet_df <- rbind(pet_temp, pet_df)
}
}

ggplot(pet_df, aes(x=date, y=pet_har_mm, colour=huc, group=huc)) + geom_line() +scale_colour_brewer(type="qual") + theme_classic()

### Combine back with prism
prism_df <- prism_df %>%
	left_join(pet_df, by= c("date", "huc"))

### And create a POSIX date column
Sys.setenv(TZ='UTC')
prism_df$date_pos <- as.POSIXct(prism_df$date, format="%Y-%m-%d")


###########################################################################
## Read in elevations and calcualte hypso curves
###########################################################################

huc_z <- c(2958, 2996, 3048, 2742, 2567)  ### elevations of centroids in meters


nhd_path <- file.path(data_path, "usgs_hydrography/WBD_16-Utah/1601/nhd_plus_raster/HRNHDPlusRasters1601")

elev_raster <- raster(file.path(nhd_path, "elev_cm.tif"))
elev_raster

### Plot to check the cutting process
plot(elev_raster, col=viridis(100))
plot(hu_12, col=NA, add=TRUE)

crs(elev_raster) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

### Takes a long time to do all, just do upper
hu_upper <- hu_12 %>%
	filter(HUC12 %in% huc_list)

### Cut out the first one, it's downstream of gauge
hu_upper <- hu_upper[c(2,3,4,5,6),]
	
### Create a reference grid with the numbers for each HUC so you don't have to redo each time
gcells <- cellnumbers(elev_raster, hu_upper) %>% 
	mutate(object_ = as.integer(object_))

### Loop through elevation data and extract al values
result_temp <- gcells %>% 
	mutate(clim_var = raster::extract(elev_raster, cell_))

for (i in seq(1, max(result_temp$object_))){
	result_i <- result_temp %>%
		filter(object_ == i)
		
	elev_i <- tibble(group=rep(i, 101)) %>%
		mutate(huc = hu_upper$HUC12[i]) %>%
		mutate(elev_m = quantile(result_i$clim_var, seq(0,1,0.01))) %>%
		mutate(elev_m = elev_m/100) %>%
		mutate(quant = seq(0,1,0.01))
	
	if (i == 1) {
		hypso_quant <- elev_i
	} else {
		hypso_quant <- rbind(hypso_quant, elev_i)
	}

}

ggplot(hypso_quant, aes(x=quant, y=elev_m, colour=huc, group=group)) + geom_line() +scale_colour_brewer(type="qual") + theme_classic()


###########################################################################
## Calculate drainage area
###########################################################################
hu_upper_area <- hu_upper %>% 
	mutate(area = st_area(.), area_m_sqr = as.numeric(area)) %>%
	st_set_geometry(NULL) %>%
	select(HUC12, area_m_sqr)
	
wat_area_msqr <- sum(hu_upper_area$area_m_sqr)

###########################################################################
## Calculate centroid
###########################################################################
### This calculates the centroid, but forces it to be within the polygon
hu_upper_cent <- st_point_on_surface(hu_upper) %>% 
	cbind(st_coordinates(.)) %>%
	st_transform(4269) %>%  ### Convert to NAD83 lat lon
	cbind(st_coordinates(.)) %>%
	st_set_geometry(NULL) %>%  ### Remove geometry
	select(HUC12, X, Y, X.1, Y.1) %>%
	rename("X_m" = "X", "Y_m" = "Y", "Lon" = "X.1", "Lat" = "Y.1")

head(hu_upper_cent)


###########################################################################
## Create merged huc characteristics
###########################################################################
### For now use just the PRISM elevations
hu_upper_z <- tibble(HUC12 = huc_list, elev_m = huc_z)

hu_char <- hu_upper_area %>%
	full_join(hu_upper_cent) %>%
	full_join(hu_upper_z) %>%
	rename("huc" = "HUC12")

hu_char


###########################################################################
###  Could use this to save instead
###########################################################################
save(hu_char, hypso_quant, prism_df, file=file.path(write_output_path, "upper_huc_params.rda"))


