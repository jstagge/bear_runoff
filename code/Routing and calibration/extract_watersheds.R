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
require(fasterize)

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
## Calculate centroids
###########################################################################
yup <- centroid %>% 
	st_set_geometry(NULL) %>%
	cbind(st_coordinates(centroid)) %>%
	select(HUC12, X, Y)

### This does the true centroid in meters
yup <- st_centroid(hu_12) %>% 
	cbind(st_coordinates(.)) %>%
	st_set_geometry(NULL) %>%
	select(HUC12, X, Y)

### This does the true centroid but first converts back to lat lon
yup2 <- hu_12	%>%
	st_transform(4269) %>%  ### Convert to NAD83 lat lon
	st_centroid() %>% 
	cbind(st_coordinates(.)) %>%	
	st_set_geometry(NULL)  %>%
	select(HUC12, X, Y)
	
### This does the centroid, but forces it to be within the polygon
yup2 <- st_point_on_surface(hu_12) %>% 
	cbind(st_coordinates(.)) %>%
	st_set_geometry(NULL) %>%
	select(HUC12, X, Y)



###########################################################################
## Read in SSURGO GIS layers for Hyd Soil Group
###########################################################################
### Set up the file names to read in
file_names <- c("ID700", "ID709", "ID710", "ID711", "ID712", "ID713", "ID714", "ID715", "ID716", "ID758", "UT601", "UT602", "UT603", "UT604", "UT613", "UT647", "WY041", "WY623", "WY638", "WY663", "WY723")

for (j in seq(1,length(file_names))){

	i <- file_names[j]
	
	### Read in soil
	ssurgo_temp <- st_read(paste0(soil_path, "/",i,"/spatial/soilmu_a_",tolower(i), ".shp")) %>% 
	st_transform(usgs_proj)

	### Read in attributes
	muaggatt <- read_delim(paste0(soil_path, "/",i,"/tabular/muaggatt.txt"), delim="|", col_names=FALSE) %>%
	select(c(1,2, 18, 40)) %>%
	rename(
		musym = X1,
		muname = X2,
		hydgrpdcd = X18,
		mukey = X40) %>%
	mutate(mukey = factor(mukey, levels=levels(ssurgo_temp$MUKEY)))

	ssurgo_temp <- ssurgo_temp %>%
		left_join(muaggatt, by= c("MUKEY" = "mukey")) %>%
		select(-musym)

	if(j == 1){
		ssurgo_df <- ssurgo_temp
	} else {
		ssurgo_df <- rbind(ssurgo_df, ssurgo_temp)
	}
  
}


write_sf(ssurgo_df, file.path(write_output_path, "ssurgo_hyd_group.shp"))

#ssurgo_df <- st_read(file.path(write_output_path, "ssurgo_hyd_group.shp"))

ggplot(ssurgo_temp) +
  geom_sf(aes(fill = hydgrpdcd)) +
 # scale_fill_viridis("Hydrologic Group") +
  ggtitle("Hydrologic Group") +
  theme_bw()
 
 
###########################################################################
## Read in and process SSURGO chorizon data for sma analysis
###########################################################################

### Filter out the regions with no chorizon data
valid_chorizon <- rep(FALSE,length(file_names))

for (j in seq(1,length(file_names))){

	i <- file_names[j]
	chorizon <- read_delim(paste0(soil_path, "/",i,"/tabular/chorizon.txt"), delim="|", col_names=FALSE)
	
	if (dim(chorizon)[1] > 0){valid_chorizon[j] <- TRUE}

}

file_names <- file_names[valid_chorizon]
  
for (j in seq(1,length(file_names))){

	i <- file_names[j]
	
	### Read in soil
	ssurgo_temp <- st_read(paste0(soil_path, "/",i,"/spatial/soilmu_a_",tolower(i), ".shp")) %>% 
	st_transform(usgs_proj)


	### Read in attributes
	#ksat_r Representative saturated hydraulic conductivity
	#hzdepb_r Representative depth from soil surface to bottom of layer
	#wsatiated_r Representative soil porosity
	#wthirdbar_r Representative field capacity

	chorizon <- read_delim(paste0(soil_path, "/",i,"/tabular/chorizon.txt"), delim="|", col_names=FALSE) %>%
	rename(
		hzname = X1,
		cokey = X170, 
		chkey = X171,
		ksat_r = X83,
		hzdepb_r = X10,
		wsatiated_r = X98,
		wthirdbar_r = X92) %>%
	select(cokey, chkey, hzname, ksat_r, hzdepb_r, wsatiated_r,wthirdbar_r) 
	
	### Read in components
	component <- read_delim(paste0(soil_path, "/",i,"/tabular/comp.txt"), delim="|", col_names=FALSE) %>%
	rename(
		mukey = X108,
		cokey = X109, 
		comppct_r = X2, 
		slope_r = X10) %>%
	select(mukey, cokey, comppct_r, slope_r) 
	
	### Calculate the total area of each map unit (mu) covered by each component
	mu_area_weight <- component %>%
		select(-slope_r) %>%
		group_by(mukey) %>%
		summarise(total_area = sum(comppct_r))

	### Calculate percent for each component
	component <- component %>%
		left_join(mu_area_weight) %>%
		mutate(perc_area = comppct_r/total_area) %>%
		select(-total_area, -comppct_r)

	ssurgo_sma <- component %>%
		full_join(chorizon, by="cokey") %>%
		arrange(mukey, cokey, -chkey)

	### Find average of ksat_r, wsatitated_r and wthirdbar_r for each component
	### Find the ksat_r for topmost horizon
	### Find hzdepb_r for bottommost	
	ssurgo_by_co <- ssurgo_sma %>%
		group_by(mukey, cokey) %>%
		summarise(ksat_avg = mean(ksat_r), ksat_top=first(ksat_r), hzdepb_bottom = last(hzdepb_r), wsat_avg = mean(wsatiated_r), wthirdbar_avg = mean(wthirdbar_r) )
	
	### Calculate a weighed mean for each map unit
	ssurgo_sma_summary <- ssurgo_by_co %>%
		left_join(component, by=c("mukey", "cokey"))

	ssurgo_sma_summary <- ssurgo_sma_summary %>% 
		group_by(mukey) %>%
		summarise(ksat_avg = weighted.mean(ksat_avg, perc_area), 
			ksat_top=weighted.mean(ksat_top, perc_area), 
			hzdepb_bottom = weighted.mean(hzdepb_bottom, perc_area), 
			wsat_avg = weighted.mean(wsat_avg, perc_area), 
			wthirdbar_avg = weighted.mean(wthirdbar_avg, perc_area), 
			slope= weighted.mean(slope_r, perc_area))

	ssurgo_sma_summary <- ssurgo_sma_summary %>% 
		mutate(mukey = factor(mukey, levels=levels(ssurgo_temp$MUKEY)))


	ssurgo_temp <- ssurgo_temp %>%
		left_join(ssurgo_sma_summary, by= c("MUKEY" = "mukey"))
	
	if(j == 1){
		ssurgo_sma_df <- ssurgo_temp
	} else {
		ssurgo_sma_df <- rbind(ssurgo_sma_df, ssurgo_temp)
	}
  
}

	
write_sf(ssurgo_sma_df, file.path(write_output_path, "ssurgo_sma.shp"))


ggplot(ssurgo_temp) +
  geom_sf(aes(fill = ksat_avg)) +
  scale_fill_viridis("K sat avg") +
  theme_bw()
  
ggplot(ssurgo_temp) +
  geom_sf(aes(fill = hzdepb_bottom)) +
  scale_fill_viridis("Depth to bottom horizon") +
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


