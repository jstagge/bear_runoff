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

soil_path <- "D:/SSURGO"
land_cover_path <- "D:/nlcd" 

#### Download all the PRISM data
#prism_path <- file.path(data_path, "prism")
prism_path <- "D:/prism"
options(prism.path = prism_path)

### Set up output folders
write_output_path <- file.path(output_path, "bear_runoff_model")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(raster)
require(rgdal)
require(viridis)
require(sf)
require(staggefuncs)
require(waterData)
require(prism)
require(fasterize)
require(tabularaster)

select <- dplyr::select


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


###########################################################################
###  Read in USGS Gauge data
###########################################################################
#Site Number: 10011500
#Site Name: BEAR RIVER NEAR UTAH-WYOMING STATE LINE
#Definitely use this

#Site Number: 10016900
#Site Name: BEAR RIVER AT EVANSTON, WY
#1984-05-14 	2019-03-23


#Site Number: 10020100
#Site Name: BEAR RIVER ABOVE RESERVOIR, NEAR WOODRUFF, UT
#1961-10-01 	2019-03-23

flow_border <- importDVs("10011500", code="00060", stat = "00003", sdate = "1935-01-01", edate = "2018-12-31")

my_url <- tellMeURL("10011500", code="00060", stat = "00003", sdate = "1935-01-01", edate = "2018-12-31")
my_url


my_sites <- c("10011500")
my_siteInfo <- siteInfo(my_sites)



###########################################################################
## Download Daily PRISM Climate data and find mean over each HUC region
###########################################################################


#"ppt", "tmean", "tmin", or "tmax
#Jan 1981 - Aug 2018
#get_prism_dailys(type="tmean", minDate = "1981-01-01", maxDate = "2018-08-31", keepZip=F)
get_prism_dailys(type="ppt", minDate = "1981-01-01", maxDate = "1983-08-31", keepZip=F)
get_prism_dailys(type="tmin", minDate = "1981-01-01", maxDate = "1983-08-31", keepZip=F)
get_prism_dailys(type="tmax", minDate = "1981-01-01", maxDate = "1983-08-31", keepZip=F)

### Set up a list of all the data to read in 
prism_files <- ls_prism_data(name=TRUE) %>%
	as.tibble() %>%
	separate(files, into=c("source", "variable", "flag", "res", "date", "bil"), sep="_", remove=FALSE) %>%
	mutate(date = as.Date(date, "%Y%m%d")) %>%
	arrange(date)

prism_file_var <- prism_files %>%
	filter(variable == "ppt" & flag =="stable")
### Maybe do a complete command?
	
prism_file_var <- prism_file_var %>%
	mutate(file_path = paste0(prism_path, "/", prism_file_var$files, "/", prism_file_var$files, ".bil"))


### Now load in the first raster for projection
prism_temp <- raster(prism_file_var$file_path[1])
prism_crs <- crs(prism_temp, asText=TRUE)


### Create a reference grid
### Reproject into the grid for PRISM data
hu12_reproject <- hu_12  %>% 
	st_transform(crs = crs(prism_temp, asText=TRUE))

plot(prism_temp, ext = extent(hu12_reproject))
plot(hu12_reproject, col=NA, add=TRUE)

### Create a reference grid with the numbers for each HUC so you don't have to redo each time
gcells <- cellnumbers(prism_temp, hu12_reproject) %>% 
	mutate(object_ = as.integer(object_))

### Loop through PRISM data and extract mean
for (i in seq(1, dim(prism_file_var)[1])){

prism_temp <- raster(prism_file_var$file_path[i])

result_temp <- gcells %>% 
	mutate(clim_var = raster::extract(prism_temp, cell_)) %>% 
	group_by(object_) %>% 
	summarize_at(vars(clim_var), list(mean=~mean(., na.rm = TRUE), sd=~sd(., na.rm = TRUE), n=length)) %>%
	mutate(date = prism_file_var$date[i])

if (i == 1) {
	result <- result_temp
} else {
	result <- rbind(result, result_temp)
}
}

ggplot(result, aes(x=date, y=mean, colour=object_, group=object_)) + geom_line() + scale_colour_viridis() + theme_classic() 





###########################################################################
## Read in PRISM Climate data
###########################################################################

clim_gauge <- read_csv("PRISM_ppt_tmin_tmean_tmax_stable_4km_19810101_20170101_40.7949_-110.7903.csv", skip=10, col_names = TRUE)

clim_gauge <- clim_gauge %>%
	rename(prcp_mm = `ppt (mm)`) %>%
	rename(tmin_c = `tmin (degrees C)`) %>%
	rename(tavg_c = `tmean (degrees C)`) %>%
	rename(tmax_c = `tmax (degrees C)`)
	
head(clim_gauge)
z_clim <- 2996



##to do
##get dailys for 2016 - ppt, tmean
file.path(data_path, "prism")

options(prism.path = file.path(data_path, "prism"))
get_prism_normals(type = 'tmean', resolution = '4km', mon = 12, keepZip = TRUE)

ls_prism_data(name=TRUE)

new_file<-1#this number corresponds to the row of the file of interest
RS <- prism_stack(ls_prism_data()[12,1]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") ##assign projection info

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns

year<-'30yr_norm_12'
name<-paste0("/Volumes/collnell/CAT/data/prism/monthly_temp/USA/PRISM_tmean_", year,".csv")
write.csv(m.df, name)
writeRaster(RS, name)

head(m.df)



get_prism_dailys(type = 'ppt', minDate = "2013-06-01", maxDate = "2013-06-14", keepZip = TRUE)
ls_prism_data(name=TRUE)

new_file<-c(1) ##change to corresponding file numbers
RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file
to_slice <- grep("_201607",RS[,1],value=T)##search through names


df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
month.df <- melt(df, c("x", "y"))
names(month.df)[1:2] <- c("lon", "lat") #rename columns




#### This works, kind of					   
cent_loc <- hu_results %>%
	select(HUC12, Name, Lon, Lat) %>%  ### Only retain a few columns
	st_as_sf(coords = c("Lon", "Lat"), crs = 4269)  ### Convert to simple spatial point data
	

cent_clim <- extract(RS, cent_loc,  fun=mean, na.rm=TRUE, sp=TRUE)


yup <- matrix(c(hu_results$Lon, hu_results$Lat), 2,dim(hu_results)[1], byrow=TRUE) %>%
	data.frame(.)
names(yup) <- paste0("id_",seq(1,length(hu_results$HUC12)))
stencil <- simplegeom(yup)

stencil


stencil <- simplegeom(data.frame(
              'point1' = c(-89, 46), 
              'point2' = c(-88.6, 45.2)))

stencil <- webgeom('state::New Hampshire')
stencil <- webgeom('state::New Hampshire,Wisconsin,Alabama')

stencil <- webgeom('HUC8::09020306,14060009')

stencil
## An object of class "webgeom":
## url: http://cida.usgs.gov/gdp/geoserver/wfs 
## geom: derivative:wbdhu8_alb_simp 
## attribute: HUC_8 
## values: 09020306, 14060009 
## wfs version: 1.1.0

#see what other HUCs could be used via the query function:
HUCs <- query(stencil, 'values')



query(webgeom(), "geoms")
query(webprocess(), "algorithms")
query(webdata("prism"), "variables")

"ppt" "tmx" "tmn"


fabric <- webdata(url = "dods://apdrc.soest.hawaii.edu/dods/public_data/satellite_prod

https://cida.usgs.gov/thredds/catalog.html?dataset=cida.usgs.gov/prism

https://www.esrl.noaa.gov/psd/data/gridded/data.prismday.html


fabric <- webdata(url = "dods://cida.usgs.gov/thredds/dodsC/gmo/GMO_w_meta.ncml")


prism”, “iclus”, “daymet”, “stageiv”, “topowx”, “solar”, “metobs”


data <- result(job, with.units = TRUE)

http://regclim.coas.oregonstate.edu:8080/thredds/dodsC/regcmdata/NCEP/srm/Monthly/RegCM3_Monthly_srm_NCEP.ncml

fabric <- webdata(url = "http://regclim.coas.oregonstate.edu:8080/thredds/ncss/mwbm_data/prism_mwbm_historical.ncml")

variables(fabric) <- "PRISM_historical_snow"
query(fabric, "times")




PRISM_historical_snow



fabric <- webdata(url = "http://regclim.coas.oregonstate.edu:8080/thredds/ncss/mwbm_data/prism_mwbm_historical.ncml")
times(fabric) <- c('2002-01-01','2010-01-01')
variables(fabric) <- c('PRISM_historical_snow')
fabric

, variables = "PRISM_historical_snow", times = c("1990-01-01", 
    "1995-05-01"))




fabric <- webdata("prism_daily", variables = "tmx", times = c("1895-01-01", 
    "1895-05-01"))
job <- geoknife(stencil = "state::Wisconsin,Virginia", fabric, wait = TRUE)
result(job, with.units = TRUE)

fabric <- webdata(url = "https://cida.usgs.gov/thredds/iso/prism?catalog=https%3A%2F%2Fcida.usgs.gov%2Fthredds%2Fcatalog.html&dataset=cida.usgs.gov%2Fprism")


"Yearly Conterminous U.S. actual evapotranspiration data"
OPeNDAP Service for Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 3 (Continental North America)
"Monthly Conterminous U.S. actual evapotranspiration data
Gridded Observed Meteorological Data: 1950-2010

fabric <- webdata('prism')
fabric
times(fabric) <- c('2002-01-01','2010-01-01')
variables(fabric) <- c('tmx')
fabric

job <- geoknife(stencil, fabric)
#use convienence functions to check on the job:

check(job)
running(job)
error(job)
successful(job)

data <- result(job)
plot(data[,1:2], ylab = variables(fabric))




prism_path <- file.path(data_path, "prism")
options(prism.path = file.path(data_path, "prism"))

get_prism_dailys(type="tmean", minDate = "2013-06-01", maxDate = "2013-06-14", keepZip=F)

file_list <- ls_prism_data(name=TRUE) %>%
	as.tibble() %>%
	separate(files, into=c("source", "variable", "flag", "res", "date"), sep="_", remove=FALSE) %>%
	arrange(date)

ppt_files <- file_list %>%
	filter(variable == "ppt")
### Maybe do a complete command?
	
p <- prism_slice(c(-73.2119,44.4758),ppt_files$files)
print(p)


ppt_file_list <- paste0(prism_path, "/", ppt_files$files, "/", ppt_files$files, ".bil")


## Now we'll load the rasters.
jnorm_rast <- raster(ppt_file_list[1])

j2013_rast <- raster(j2013)
## Now we can do simple subtraction to get the anomaly by subtracting 2014 from the 30 year normal map
anomCalc <- function(x, y) {
  return(x - y)
  }

anom_rast <- overlay(j2013_rast,jnorm_rast,fun = anomCalc)

plot(anom_rast)








get_prism_dailys(type="tmean", minDate = "2013-06-01", maxDate = "2013-06-14", keepZip=F)
p <- prism_slice(c(-73.2119,44.4758),file_list)
print(p)

file_list <- ls_prism_data()[1:5,1]

to_slice = grep("tmean",ls_prism_data(), value = T)
slice <- prism_slice(boulder, "PRISM_tmean_stable_4kmD1_20130603_bil.bil"

get_prism_dailys(type="tmean", minDate = "2013-06-01", maxDate = "2013-06-14", keepZip=FALSE)
p <- prism_slice(c(-73.2119,44.4758),ls_prism_data())
print(p)

to_slice <- grep("_[0-9]{4}[0][1]",ls_prism_data()[,1],value=T)













# First, you need to query for all web data
all_webdata <- query("webdata")
all_titles <- title(all_webdata)
all_abstracts <- abstract(all_webdata)

# Then start sleuthing using `grep`
#keyword_str <- "sea level rise|sea level|sea level|sea rise"
keyword_str <- "precip"
keyword_str <- "daily"
index_t <- grep(keyword_str, all_titles, ignore.case = TRUE)
index_a <- grep(keyword_str, all_abstracts, ignore.case = TRUE)
index_both <- unique(c(index_t, index_a))

# Look at the titles and abstracts of datasets that match your criteria
length(index_both)

all_titles[index_both]

all_abstracts[index_both]

index_both[44]
index_both[16]
index_both[22]

all_titles[index_both[76]]
all_abstracts[76]

all_titles[110]
all_abstracts[110]

all_titles[142]
all_abstracts[142]


"OPeNDAP Service for Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 3 (Continental North America)"

 "University of Idaho Daily Meteorological data for continental US"




# First, you need to query for all web data
all_webdata <- query("webdata")

# Setup the maximum air temp fabric using the URL in all_webdata
#airtemp_title <- "TopoWx: Topoclimatic Daily Air Temperature Dataset for the Conterminous United States"
airtemp_title <- "OPeNDAP Service for Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 3 (Continental North America)"
airtemp_url <-  url(all_webdata[airtemp_title])
airtemp_fabric <- webdata(list(
  url = airtemp_url,
  variables = "tmax",
  times = as.POSIXct(c("2007-07-04 12:01", "2007-07-05 11:59"), tz = "UTC")
))

# Now setup the stencil
texas <- webgeom(geom = "sample:CONUS_states", 
                 attribute = "STATE",
                 values = "Texas")

# Leave the default knife since we want an average over the stencil
# Execute the geoknife job
airtemp_job <- geoknife(stencil = texas, fabric = airtemp_fabric, wait=TRUE)

# Download the data
air_max_data <- result(airtemp_job)
air_max_data




# First, you need to query for all web data
all_webdata <- query("webdata")

# Use the all_webdata object to create the fabric
us_meterology <- webdata(all_webdata["OPeNDAP Service for Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 3 (Continental North America)"])
### 1980 to the present

# Now use query to see what variables are available
metero_vars <- query(us_meterology, "variables")
metero_vars


#### Set up stencil based on HUC12
stencil <- webgeom(url = "http://cida.usgs.gov/nwc/geoserver/WBD/ows")
query(stencil,'geoms')
# [1] "WBD:huc08"    "WBD:huc12"    "WBD:huc12agg" "WBD:huc12all"
geom(stencil) <- 'WBD:huc12'
query(stencil, 'attributes')
# [1] "ogc_fid"    "tnmid"      "metasource" "sourcedata" "sourceorig" "sourcefeat" "loaddate"   "gnis_id"   
# [9] "areaacres"  "areasqkm"   "states"     "huc12"      "hutype"     "humod"      "tohuc"      "noncontrib"
# [17] "noncontr_1" "shape_leng" "shape_area"
attribute(stencil) <- 'huc12'
# See: http://cida.usgs.gov/nwc/#!waterbudget/achuc/020600031101 or other lookup for where to find a HUC.
values(stencil) <- '160102010200'

### Set up fabric
#fabric <- webdata("metobs", variables = "pr", times = c("1990-01-01", "2005-05-01"))
#fabric <- webdata("metobs", variables = "pr", times = c("2002-01-01", "2005-05-01"))
fabric <- webdata(url = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc", variables = "precipitation_amount", times = c("1990-01-01", "2005-05-01"))

job <- geoknife(stencil = stencil, fabric, wait = TRUE)
data.out <- result(job, with.units = TRUE)

plot(data.out[,1:2], ylab = variables(fabric))



fabric = webdata(url = "https://thredds.daac.ornl.gov/thredds-daymet/dodsC/daymet-v3-agg/na.ncml", 
    times = c("1990-01-01", "2010-01-01"), variables = c("prcp", "tmax"))
fabric

fabric = webdata(url = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_pr_1979_CurrentYear_CONUS.nc")

precipitation_amount


http://services.nacse.org/prism/data/public/4km/tmin/20090405

http://services.nacse.org/prism/data/public/releaseDate/ppt/19990101/19991231













library(geoknife)
# Use the WBD feature service available from the National Water Census.
stencil <- webgeom(url = "http://cida.usgs.gov/nwc/geoserver/WBD/ows")
query(stencil,'geoms')
# [1] "WBD:huc08"    "WBD:huc12"    "WBD:huc12agg" "WBD:huc12all"
geom(stencil) <- 'WBD:huc12'
query(stencil, 'attributes')
# [1] "ogc_fid"    "tnmid"      "metasource" "sourcedata" "sourceorig" "sourcefeat" "loaddate"   "gnis_id"   
# [9] "areaacres"  "areasqkm"   "states"     "huc12"      "hutype"     "humod"      "tohuc"      "noncontrib"
# [17] "noncontr_1" "shape_leng" "shape_area"
attribute(stencil) <- 'huc12'
# See: http://cida.usgs.gov/nwc/#!waterbudget/achuc/020600031101 or other lookup for where to find a HUC.
values(stencil) <- '160102010200'
# Run a short time when flash flooding happened.
startDate <- "2016-07-30 22:00:00"
endDate <- "2016-07-31 01:00:00"
fabric <- webdata(url = 'http://cida.usgs.gov/thredds/dodsC/stageiv_combined')

,
                  variables = "Total_precipitation_surface_1_Hour_Accumulation",
                  times = c(as.POSIXct(startDate),
                            as.POSIXct(endDate)))
job <- geoknife(stencil, fabric, wait=TRUE)
data.out <- result(job, with.units=TRUE)
data.out <- result(job, with.units=TRUE)
data.out





job <- geoknife(stencil = "state::Wisconsin,Virginia", fabric, wait = TRUE)
result(job, with.units = TRUE)


daymet, gldas", "nldas", "topowx", "solar", "metobs"

fabric <- webdata('metobs')
fabric
#query(fabric, "variables")

variables(fabric) <- c('pr')
times(fabric) <- as.POSIXct(c("2000-07-04", "2007-07-05"), tz = "UTC")


as.POSIXct(c("2007-07-04 12:01", "2007-07-05 11:59"), tz = "UTC")
fabric

as.POSIXct(c("2000-07-04", "2005-07-05"))

webdata('prism', times=as.POSIXct(c('1990-01-01', '1995-01-01')))


job <- geoknife(stencil, fabric)
#use convienence functions to check on the job:

check(job)
running(job)
error(job)
successful(job)

data <- result(job)
plot(data[,1:2], ylab = variables(fabric))






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

###########################################################################
## Calculate centroids
###########################################################################

### This calculates the centroid, but forces it to be within the polygon
huc12_cent <- st_point_on_surface(hu_12) %>% 
	cbind(st_coordinates(.)) %>%
	st_transform(4269) %>%  ### Convert to NAD83 lat lon
	cbind(st_coordinates(.)) %>%
	st_set_geometry(NULL) %>%  ### Remove geometry
	select(HUC12, X, Y, X.1, Y.1) %>%
	rename("X_m" = "X", "Y_m" = "Y", "Lon" = "X.1", "Lat" = "Y.1")

head(huc12_cent)

### Quick check plot
plot(hu_12["HUC12"], add=TRUE)
points(huc12_cent$X_m, huc12_cent$Y_m,  add=TRUE)


###########################################################################
## Calculate area
###########################################################################
huc12_area <- hu_12 %>% 
	mutate(area = st_area(.), area_m_sqr = as.numeric(area)) %>%
	st_set_geometry(NULL) %>%
	select(HUC12, area_m_sqr)
	
###########################################################################
## Read in Curve Number
###########################################################################
cn_raster <- raster(file.path(write_output_path, "curve_number.tif"))
cn_raster

### Plot to check the cutting process
plot(cn_raster, col=viridis(100))
plot(hu_12, col=NA, add=TRUE)

###########################################################################
## Calculate mean Curve Number for each HUC 12 unit
###########################################################################
### Extract the min, max, and mean CN for each HUC12 region
cn_byhu = raster::extract(x = cn_raster, y = hu_12, df = TRUE) 

### Numbering scheme is how they were originally written
### Create a table to hold the results
hu_frame <- data.frame(HU12=hu_12$HUC12) %>%
	mutate(HU12 = factor(HU12, levels=levels(hu_12$HUC12))) %>%
	mutate(ID = seq(1,dim(.)[1]))

cn_byhu_df <- cn_byhu %>%
	as.tibble() %>%
	left_join(hu_frame) %>%
	mutate(HU12 = as.character(HU12)) %>%
	group_by(HU12) %>% 	
	summarize(min_cn = min(curve_number, na.rm=TRUE), 
		mean_cn = mean(curve_number, na.rm=TRUE),
		max_cn = max(curve_number, na.rm=TRUE))

cn_byhu_df <- hu_12 %>%
	left_join(cn_byhu_df, by=c("HUC12" = "HU12"))

plot(cn_byhu_df["mean_cn"])	

### Write each cell CN in case you have to come back
write.csv(cn_byhu, file.path(write_output_path, "cn_byhu.csv"))


###########################################################################
## Extract the important parts and start to write out
###########################################################################
### Create a base tibble
hu_results <- hu_12 %>%
	st_set_geometry(NULL) %>%  ### Remove geometry
	select(HUC12, ToHUC, HUType, HUMod, Name, States) %>%
	as.tibble()

### Merge with centroid
hu_results <- hu_results %>%
	left_join(huc12_cent)
	
### Merge with drainage area
hu_results <- hu_results %>%
	left_join(huc12_area)

### Merge with Curve Number
hu_results <- hu_results %>%
	left_join(select(cn_byhu_df, HUC12, min_cn,  mean_cn, max_cn))

hu_results <- hu_results %>%
	select(-Shape)
	
head(hu_results)
 
 
### Write each cell CN in case you have to come back
write.csv(hu_results, file.path(write_output_path, "watershed_char.csv"))


###########################################################################
## Read in SSURGO and STATSGO layers for SMA
###########################################################################

ssurgo_sma_df <- st_read(file.path(write_output_path, "ssurgo_sma.shp"))
statsgo_sma_df <- st_read(file.path(write_output_path, "statsgo_sma.shp"))



ssurgo_sma_df %>%
	st_transform(usgs_proj)

ggplot(statsgo_sma_df) +
  geom_sf(aes(fill = ksat_tp)) +
  scale_fill_viridis("Ksat top layer")
  theme_bw()




###########################################################################
## Calculate mean SMA properties for each HUC 12 unit
###########################################################################
### First calculate the intersection of HUC regions and ssurgo map units
huc_ssurgo_inter = st_intersection(hu_12, ssurgo_sma_df)

### Calculate surface area for each new intersected area
huc_ssurgo_inter <- huc_ssurgo_inter %>%
	mutate(area = st_area(.))

### Summarize the area-weighted mean for each HUC unit
### Need to remove the geometry to make this calculation correct
sma_summ_byhuc <- huc_ssurgo_inter %>% 
	st_set_geometry(NULL) %>%
	group_by(HUC12) %>%
	summarize(ksat_avg = weighted.mean(ksat_avg, area, na.rm = TRUE),
		ksat_top = weighted.mean(ksat_top, area, na.rm = TRUE),
		hzdepb_bottom = weighted.mean(hzdepb_bottom, area, na.rm = TRUE) ,
		wsat_avg = weighted.mean(wsat_avg, area, na.rm = TRUE) ,
		wthirdbar_avg = weighted.mean(wthirdbar_avg, area, na.rm = TRUE),
		slope = weighted.mean(slope, area, na.rm = TRUE)
	)

### Merge back with geometry of HUC units
hu_12_summary <- hu_12 %>%
	left_join(sma_summ_byhuc)
	
### Test plot
plot(hu_12_summary["slope"])	
plot(hu_12_summary["wsat_avg"])	
plot(hu_12_summary["ksat_top"])	





###########################################################################
## Read in Land Cover GIS layers
###########################################################################
land_cover <- raster(paste0(land_cover_path, "/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img"))
land_cover

plot(land_cover)

### Create a bounding box based on bear river watershed
### Work with the original non-projected data and project to US grid
bound_box <- as(extent(hu_4), 'SpatialPolygons')

bound_box <- st_as_sf(bound_box) %>%  # convert polygons to 'sf' object
	st_set_crs(usgs_proj) %>% ### Original HUC data is in 4269 projection NAD83 (EPSG:4269) 
	st_buffer(20000) %>%  ### 20,000 meters or about 12 miles
	st_transform(42303)  ### Land cover is in NAD83 / Albers NorthAm

### Crop the data
land_cover <- crop(land_cover, bound_box)

### Project raster to USGS projection
land_cover <- projectRaster(land_cover, crs=st_crs(hu_4)$proj4string)

### Plot to check the cutting process
plot(land_cover)
plot(hu_12, add=TRUE)

### Write the cut raster
writeRaster(land_cover, file.path(write_output_path, "land_cover.tif"))

#land_cover <- raster(file.path(write_output_path, "land_cover.tif"))

###########################################################################
## Calculate CN
###########################################################################
### Rasterize hydrologic soil group

### First have to deal with breaks
### Possible to do as both?  Or carry through as split and apply later

### Reclassify land cover



ssurgo_df_orig <- ssurgo_df

hyd_group_levels <- levels(ssurgo_df$hydgrpdcd)
ssurgo_df <- ssurgo_df %>%
	mutate(hydgrpdcd = factor(ssurgo_df$hydgrpdcd, levels=hyd_group_levels)) %>%
	mutate(hyd_numer = as.numeric(hydgrpdcd))

hyd_group_raster <- fasterize(ssurgo_df, land_cover, field =  "hyd_numer" , fun = "min")



### Come back to this
# mask Raster* with sf object
#mb <- mask(r, b)

nlcd_reclass <- read.csv("nlcd_reclass.csv")
curve_number <- read.csv("curve_number.csv")

curve_number <- curve_number %>% 
	select(-NLCD) %>% 
	gather(soil, curve, -ID) %>%
	rename("land_cover" = "ID") %>%
	mutate(soil_num = as.numeric(factor(soil, levels=hyd_group_levels)))

### Take the average of A and D
ad_cn <- curve_number %>% 
	filter(soil == "A" | soil == "D") %>%
	group_by(land_cover) %>%
	summarize(curve = mean(curve, na.rm=TRUE)) %>%
	mutate(soil = "A/D") %>%
	mutate(soil_num = as.numeric(factor(soil, levels=hyd_group_levels))) %>%
	select(land_cover, soil, curve, soil_num)

### Take the average of B and D
bd_cn <- curve_number %>% 
	filter(soil == "B" | soil == "D") %>%
	group_by(land_cover) %>%
	summarize(curve = mean(curve, na.rm=TRUE)) %>%
	mutate(soil = "B/D") %>%
	mutate(soil_num = as.numeric(factor(soil, levels=hyd_group_levels))) %>%
	select(land_cover, soil, curve, soil_num)

### Take the average of C and D
cd_cn <- curve_number %>% 
	filter(soil == "C" | soil == "D") %>%
	group_by(land_cover) %>%
	summarize(curve = mean(curve, na.rm=TRUE)) %>%
	mutate(soil = "C/D") %>%
	mutate(soil_num = as.numeric(factor(soil, levels=hyd_group_levels))) %>%
	select(land_cover, soil, curve, soil_num)


curve_number <- rbind(curve_number, ad_cn, bd_cn, cd_cn)


cn_lookup <- function(lc_val, soil_val){
	temp_df <- tibble(land_cover = lc_val, soil_num=soil_val)
	temp_df <- temp_df %>% left_join(curve_number, by = c("land_cover", "soil_num"))
	return(temp_df$curve)
}

### Convert land_cover and hydrologic soil group to vectors
land_cover_vec <- c(as.matrix(land_cover))
soil_vec <- c(as.matrix(hyd_group_raster))

### Run the CN lookup
cn_comb <- cn_lookup(land_cover_vec, soil_vec)

### Create Curve Number raster
cn_raster <- matrix(cn_comb, dim(hyd_group_raster))
cn_raster <- raster(cn_raster)
extent(cn_raster) <- extent(hyd_group_raster)
projection(cn_raster) <- projection(hyd_group_raster)

plot(cn_raster)

### Write as a test
writeRaster(cn_raster, file.path(write_output_path, "curve_number.tif"), overwrite=TRUE)

### Extract the min, max, and mean CN for each HUC12 region
zion_srtm_values = raster::extract(x = cn_raster, y = hu_12, df = TRUE) 

#### Not sure if ID is in order of levels or how they are written
### It is how they were originally written
hu_frame <- data.frame(HU12=hu_12$HUC12) %>%
	mutate(HU12 = factor(HU12, levels=levels(hu_12$HUC12))) %>%
	mutate(ID = seq(1,dim(.)[1]))

huh <- zion_srtm_values %>%
	as.tibble() %>%
	left_join(hu_frame) %>%
	mutate(HU12 = as.character(HU12)) %>%
	group_by(HU12) %>% 	
	summarize(min_cn = min(layer, na.rm=TRUE), 
		mean_cn = mean(layer, na.rm=TRUE),
		max_cn = max(layer, na.rm=TRUE))

yup <- hu_12 %>%
	left_join(huh, by=c("HUC12" = "HU12"))

plot(yup["mean_cn"])	
	
################### This worked
huh














cn_lookup(c(11, 21, 41), c("A", "B", "F"))


r3 <- overlay(land_cover, hyd_group_raster, fun=cn_lookup)

### Says it won't work because not vectorized

v = raster(land_cover); v[]=NA
v[] = merge(data.frame(land_cover=land_cover[], soil_num=hyd_group_raster[]), df, all.x=TRUE)$curve


lc_test <- land_cover[] == 82
plot(land_cover[lc_test])

soil_test <- hyd_grou_raster[] == 1
plot(hyd_grou_raster[soil_test])





# Create raster stack 
r <- raster(ncol=50, nrow=50) 
s <- stack(lapply(1:10, function(x) setValues(r, rep(1:100, each=(ncell(r)/100))))) 

# Create lookup table 
m <- as.data.frame(matrix(sample(100, 100*10, TRUE), 100, 10)) 

# create new stack with lookup 
new.stack <- stack(mapply(FUN = function(x, y) { 
  lookup <- matrix(cbind(1:length(y), y), ncol=2) 
  reclassify(x, lookup) 
}, x = as.list(s), y = m)) 

plot(s); plot(new.stack) 












###########################################################################
## Don't need this
###########################################################################
sma_nogeom <- ssurgo_sma_df %>%
	st_set_geometry(NULL) %>%
	select(-AREASYMBOL, -SPATIALVER, -MUSYM)
	
yup <- ssurgo_df %>%
	left_join(sma_nogeom, by="MUKEY")






huc4.ps <- SpatialPolygons2PolySet(huc4.polygons)




### Come back and figue this out later
bear_soilmu_a_merge.shp


ssurgo_df <- readOGR(soil_path, "bear_soilmu_a_merge")



library(raster)
# read data    
p <- shapefile("path/file.shp")
d <- read.csv("path/file.csv")

# merge on common variable, here called 'key'
m <- merge(p, d, by='key')

# perhaps save as shapefile again
shapefile(m, "path/merged.shp")


shapefile = readOGR(dsn = "DIRECTORY WITH SHAPEFILES", layer = "THE ACTUAL SHAPEFILE")
shapefile@data$id = rownames(shapefile@data)
shapefile.points = fortify(shapefile, region = "id")
shapefile.df = join(shapefile.points, shapefile@data, by = "id")
shapefile.df = subset(shapefile.df, select = c(long, lat, group, GEOID))
names(shapefile.df) = c("long", "lat", "group", "GEOID")








###########################################################################
## Set years to read in from ArcGIS
###########################################################################
lc_source <- data.frame(source="Backcast", year=seq(1940, 1990, 10))
lc_source <- rbind(lc_source, data.frame(source="Historical", year=c(1992,2000,2005)))
lc_source <- rbind(lc_source, data.frame(source="NLCD", year=c(2001,2006,2011)))

## Add in all HUC regions
#lc_source$huc <- c(rep("huc_4", 1),rep("huc_5", 11) )
lc_source$huc <- "huc_4"
lc_source <- lc_source %>% 
	complete(nesting(source, year), huc) %>%
	distinct() %>%
	arrange(huc, year)


###########################################################################
## Read in HUC watersheds
###########################################################################

### huc2 = region
### huc4 = subregion
### huc6 = basin
### huc8 = subbasin
### huc10 = watershed
### huc12 = subwatershed

landc_path <- file.path(data_path, "ohio_land_cover/processed/huc12")


### Read in HUC 4 at resolution 12 and remove unnecesary columns
huc_4_12 <- read_csv(file.path(landc_path, "huc_4_12_info.txt"))
huc_4_12 <- huc_4_12 %>%
	select(FID, OBJECTID, NonContrib, NonContr_1, AreaSqKm, AreaAcres, GNIS_ID, Name, States, HUC12, HUType, HUMod, ToHUC)

head(huc_4_12)

### Organize into columns that have each HUC resolution
### Leave the original intact
### If you read in multiple, this is where you rbind
huc_info <- huc_4_12

huc_info <- huc_info %>%
	mutate(hu_2 = substr(HUC12, 1, 2)) %>%
	mutate(hu_4 = substr(HUC12, 1, 4)) %>%
	mutate(hu_6 = substr(HUC12, 1, 6)) %>%
	mutate(hu_8 = substr(HUC12, 1, 8)) %>%
	mutate(hu_10 = substr(HUC12, 1, 10)) %>%
	mutate(hu_12 = substr(HUC12, 1, 12)) %>%
	select(-FID, -OBJECTID, -NonContrib, -NonContr_1, -GNIS_ID, -AreaSqKm, -AreaAcres) %>%
	arrange(HUC12)

head(huc_info)

###########################################################################
## Read in LandCover Codes
###########################################################################
sce_lc_codes <- read.csv(file.path(data_path, "land_cover_fore_sce/landcover_codes.csv"))
#nlcd_lc_codes <- read.csv(file.path(data_path, "nlcd/nlcd_landcover_codes.csv"))
nlcd_lc_codes <- read.csv("D://nlcd/nlcd_landcover_codes.csv")

sce_lc_codes <- sce_lc_codes %>% 
	left_join(select(nlcd_lc_codes, level1, land_cover_level1), by = c("nlcd_level1" = "level1")) %>%
	distinct()




###########################################################################
## Calculate centroids
###########################################################################
st_centroid(x)












require(rgdal)
require(maptools)
require(sp)

# Read Hydrologic Unit Code polygons, convert to PBSMapping PolySet object.

huc4.polygons <- readOGR(".", "HucsLower48Simple")
huc4.ps <- SpatialPolygons2PolySet(huc4.polygons)

# Calculate the centroids 
 
huc4.centroids <- calcCentroid(huc4.ps, rollup=1)

# Combine the centroid coordinates with the HUC4 codes,
# and promote to a Spatial Points Data Frame

huc4.with.cent <- cbind(huc4.centroids, huc4.polygons)
coordinates(huc4.with.cent) <- c("X", "Y")

# Write the updated HUC cetroids to a CSV file.

write.csv(huc4.with.cent, file="Huc4WithCentroids.csv", row.names=FALSE)



#This works
#ssurgo_df <- st_read(file.path(soil_path, "ID709/spatial/soilmu_a_id709.shp"))
ssurgo_df <- st_read(file.path(soil_path, "bear_soilmu_a_merge.shp"))

# convert to SpatialPolygonsDataFrame
#nc_sp <- as(nc, "Spatial")

muaggatt_col <- c("musym"
, "muname"
, "mustatus"
,"slopegraddcp"
,"slopegradwta"
,"brockdepmin"
,"wtdepannmin"
,"wtdepaprjunmin"
,"flodfreqdcd"
,"flodfreqmax"
,"pondfreqprs"
,"aws025wta"
,"aws050wta"
,"aws0100wta"
,"aws0150wta"
,"drclassdcd"
,"drclasswettest"
,"hydgrpdcd"
,"iccdcd"
,"iccdcdpct"
,"niccdcd"
,"niccdcdpct"
,"engdwobdcd"
,"engdwbdcd"
,"engdwbll"
,"engdwbml"
,"engstafdcd"
,"engstafll"
,"engstafml"
,"engsldcd"
,"engsldcp"
,"englrsdcd"
,"engcmssdcd"
,"engcmssmp"
,"urbrecptdcd"
,"urbrecptwta"
,"forpehrtdcp"
,"hydclprs"
,"awmmfpwwta"
,"mukey")

muaggatt <- read_delim(file.path(soil_path, "ID709/tabular/muaggatt.txt"), delim="|", col_names=muaggatt_col)

muaggatt <- muaggatt%>%
	mutate(MUKEY = factor(mukey, levels=levels(ssurgo_df$MUKEY))) %>%
	select(MUKEY, muname, hydgrpdcd)
	
### This is the Hydrologic Soil Group
### hydgrpdcd


ssurgo_df <- ssurgo_df %>%
	full_join(muaggatt, by="MUKEY")

plot_df <- ssurgo_df %>% 
	filter(AREASYMBOL == "ID709")
	
ggplot(plot_df) +
  geom_sf(aes(fill = hydgrpdcd)) +
 # scale_fill_viridis("Hydrologic Group") +
  ggtitle("Hydrologic Group") +
  theme_bw()
  
  
  
  
  
  
  



#### Pull out CN and statsgo stuff


