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
prism_path <- file.path(data_path, "prism")

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
require(airGR)
require(waterData)

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


huc_list <- c("160101010101", "160101010102", "160101010103", "160101010104", "160101010106")

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

### Convert names and select columns
flow_border <- flow_border %>%
	rename("date" = "dates") %>%
	select(-qualcode)


my_url <- tellMeURL("10011500", code="00060", stat = "00003", sdate = "1935-01-01", edate = "2018-12-31")
my_url


my_sites <- c("10011500")
my_siteInfo <- siteInfo(my_sites)



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
		mutate(huc = huc_list[i])
	
	if (i == 1) {
		prism_df <- prism_read
	} else {
		prism_df <- rbind(prism_df, prism_read)
	}

}

huc_z <- c(2958, 2996, 3048, 2742, 2567)  ### elevations of centroids in meters



###########################################################################
## Read in elevations and calcualte hypso curves
###########################################################################
require(tabularaster)

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

### Loop through PRISM data and extract mean
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
## Convert observed gauge to mm per day
###########################################################################

flow_obs <- flow_border %>%
	select(dates, val) %>%
	mutate(flow_ls = val * 28.316846592) %>%  ### convert to l/s
	mutate(flow_cuft_day = val * 86400) %>% ### Convert to cubic feet per day
	mutate(flow_cum_day = flow_cuft_day * .0283168 ) %>% ### convert to m3 per day
	mutate(flow_m_day = flow_cum_day/wat_area_msqr) %>%  ### Divide by watershed area in m2 to get m per day runoff
	mutate(flow_mm_day = flow_m_day * 1000) %>% ### Convert to mm per day
	select(dates, flow_ls, flow_mm_day) %>%
	rename("date" = "dates") %>%
	complete(date = seq.Date(min(date), max(date), by="day")) %>%
	arrange(date)
	

###########################################################################
## Handle Climate data
###########################################################################

### Need to separate each one because the function can't handle different Latitudes

for (i in seq(1, dim(hu_char)[1])){

lat_temp <- hu_char$Lat[i]

pet_temp <- prism_df %>%
	filter(huc == hu_char$huc[i]) %>%	
	rename("date" = "Date") %>%
	left_join(hu_char) %>%
	select(date, huc, tavg_c, Lat) %>%
	mutate(JD = as.POSIXlt(date)$yday) %>%
	drop_na() %>%
	mutate(pet = PEdaily_Oudin(JD = JD, Temp = tavg_c,  Lat = lat_temp, LatUnit = "deg")) %>%
	select(date, huc, pet)
	
	
if (i == 1) {
	pet_df <- pet_temp
} else {
	pet_df <- rbind(pet_temp, pet_df)
}
}

ggplot(pet_df, aes(x=date, y=pet, colour=huc, group=huc)) + geom_line() +scale_colour_brewer(type="qual") + theme_classic()

### Combine back with prism
prism_df <- prism_df %>%
	rename("date" = "Date")

prism_df <- prism_df %>%
	left_join(pet_df, by= c("date", "huc"))

### And create a POSIX date column
Sys.setenv(TZ='UTC')
prism_df$date_pos <- as.POSIXct(prism_df$date, format="%Y-%m-%d")

	
###########################################################################
## Try to set up a run
###########################################################################

model_param <- tibble(huc= hu_char$huc) %>%
	mutate(X1 = 220) %>%
	mutate(X2 = 2.127) %>%
	mutate(X3 = 125.116) %>%
	mutate(X4 = 20) %>%
	mutate(X5  =0.01) %>%
	mutate(X6 = 50)
	
model_param
	
	
for(i in seq(1,length(hu_char$huc))){

huc_i <- hu_char$huc[i]

### Cut Prism Climate data to HUC
prism_i <- prism_df %>%
	filter(huc == huc_i)

### Cut hypso data to HUC
hypso_i <- hypso_quant %>%
	filter(huc == huc_i) %>%
	arrange(quant)

### Cut z for gauge
z_clim <- hu_char$elev_m[i]

### Create input model
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, 
	DatesR = prism_i$date_pos, 
	Precip = prism_i$prcp_mm, 
	PotEvap = prism_i$pet, 
	TempMean = prism_i$tavg_c,
	TempMin = prism_i$tmin_c, 
	TempMax = prism_i$tmax_c, 
	ZInputs = z_clim, 
	HypsoData = hypso_i$elev_m, 
	NLayers = 10)
	
str(InputsModel)

 ### Create index period for run 
Ind_Run <- seq(which(format(prism_i$date_pos, format = "%d/%m/%Y")=="01/01/1981"), 
               which(format(prism_i$date_pos, format = "%d/%m/%Y")=="31/12/2016"))
str(Ind_Run)

### Cut model parameters to HUC 
param_i <- model_param  %>%
	filter(huc == huc_i) %>%
	select(-huc) %>%
	unlist()

## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel, IndPeriod_Run = Ind_Run)
	
OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i)
	
	
	
	
OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i)


InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeige, 
	DatesR = prism_i$date_pos, 
	Precip = prism_i$prcp_mm, 
	PotEvap = prism_i$pet, 
	TempMean = prism_i$tavg_c,
	TempMin = prism_i$tmin_c, 
	TempMax = prism_i$tmax_c, 
	ZInputs = z_clim, 
	HypsoData = hypso_i$elev_m, 
	NLayers = 10)
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeige, InputsModel = InputsModel, IndPeriod_Run = Ind_Run)
SnowOutputsModel <- RunModel_CemaNeige(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i[5:6])


plot(OutputsModel$Precip[1:500], type="l")
lines(OutputsModel$Qsim[1:500], col="red")


lines(OutputsModel$Qsim[1:500], col="red")

]$PliqAndMelt

SnowOutputsModel

names(OutputsModel$CemaNeigeLayers$Layer01)


plot(OutputsModel$CemaNeigeLayers$Layer01$Pliq[1:600], type="l")


 plot(OutputsModel$CemaNeigeLayers$Layer01$Pliq[1:600], type="l")
 lines(OutputsModel$CemaNeigeLayers$Layer02$Pliq[1:600], col="red")
 lines(OutputsModel$CemaNeigeLayers$Layer10$Pliq[1:600], col="green")

 plot(OutputsModel$CemaNeigeLayers$Layer01$Pliq[2:601], type="l")
 lines(SnowOutputsModel$CemaNeigeLayers$Layer01$Pliq[1:600], col="red")
 
 plot(OutputsModel$CemaNeigeLayers$Layer01$PliqAndMelt[2:201], type="l")
 lines(SnowOutputsModel$CemaNeigeLayers$Layer01$PliqAndMelt[1:200], col="red")


 plot(OutputsModel$CemaNeigeLayers$Layer10$Pliq[2:601], type="l")
 lines(SnowOutputsModel$CemaNeigeLayers$Layer10$Pliq[1:600], col="red")
 
 
$PotEvap [numeric] series of input potential evapotranspiration [mm/d]
$Precip


plot(OutputsModel$Precip[2:301], type="l")
lines(SnowOutputsModel$CemaNeigeLayers$Layer05$Pliq[1:300] + SnowOutputsModel$CemaNeigeLayers$Layer05$Psol[1:300], col="red")


	Precip = c(0, SnowOutputsModel$CemaNeigeLayers$Layer05$PliqAndMelt),

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, 
	DatesR = prism_i$date_pos, 
	Precip = c(0, yup2), 
	PotEvap = prism_i$pet, 
	TempMean = prism_i$tavg_c,
	TempMin = prism_i$tmin_c, 
	TempMax = prism_i$tmax_c, 
	ZInputs = z_clim, 
	HypsoData = hypso_i$elev_m, 
	NLayers = 10)
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel, IndPeriod_Run = Ind_Run)
PrecipOutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i[1:4])


plot(OutputsModel$Qsim[2:301], type="l")
lines(PrecipOutputsModel$Qsim[2:301], col="red")

plot(OutputsModel$Qsim[10002:10301], type="l")
lines(PrecipOutputsModel$Qsim[10002:10301], col="red")


CemaNeigeStateEnd[seq_len(2 * NLayers)[seq_len(2 * NLayers)%%2 == 1]]
				
				
seq_len(2 * 10)[seq_len(2 * 10)%%2 == 1]




yup <- data.frame(Layer01 = SnowOutputsModel$CemaNeigeLayers$Layer01$PliqAndMelt, 
Layer02 = SnowOutputsModel$CemaNeigeLayers$Layer02$PliqAndMelt, 
Layer03 = SnowOutputsModel$CemaNeigeLayers$Layer03$PliqAndMelt, 
Layer04 = SnowOutputsModel$CemaNeigeLayers$Layer04$PliqAndMelt, 
Layer05 = SnowOutputsModel$CemaNeigeLayers$Layer05$PliqAndMelt, 
Layer06 = SnowOutputsModel$CemaNeigeLayers$Layer06$PliqAndMelt, 
Layer07 = SnowOutputsModel$CemaNeigeLayers$Layer07$PliqAndMelt, 
Layer08 = SnowOutputsModel$CemaNeigeLayers$Layer08$PliqAndMelt, 
Layer09 = SnowOutputsModel$CemaNeigeLayers$Layer09$PliqAndMelt, 
Layer10 = SnowOutputsModel$CemaNeigeLayers$Layer10$PliqAndMelt)

yup2 <- apply(yup,1,mean, na.rm=TRUE)



#### Close enough






## results preview
plot(OutputsModel)

#OutputsModel$QR [numeric] series of routing store outflow (QR) [mm/d]
#OutputsModel$QD [numeric] series of direct flow from UH2 after exchange (QD) [mm/d]
#OutputsModel$Qsim [numeric] series of simulated discharge [mm/d]

qsim_temp <- tibble(date = as.Date(OutputsModel$DatesR), huc = huc_i, qsim_mm = OutputsModel$Qsim)

sub_area_msqr <- hu_char$area_m_sqr[i]

qsim_temp <- qsim_temp %>%
	mutate(qsim_m_day = qsim_mm / 1000) %>% ### Convert to m per day
	mutate(qsim_cum_day = qsim_m_day * sub_area_msqr) %>%  ### Multiply by watershed area in m2 to get m3 per day runoff
	mutate(flow_cuft_day = qsim_cum_day / .0283168 ) %>% ### convert to ft3 per day
	mutate(flow_cfs = flow_cuft_day / 86400)  ### Convert to cubic feet per day


if (i == 1) {
	qsim_df <- qsim_temp
} else {
	qsim_df <- rbind(qsim_df, qsim_temp)
}

}


yup <- qsim_df %>%
	select(date, huc, flow_cfs) %>%
	spread(huc, flow_cfs)

yup$at_gauge <- apply(yup[,2:6],1,sum)


p <- ggplot(yup, aes(x=date, y=at_gauge)) + geom_line() + theme_classic()
p + geom_line(data = flow_border, aes(x=dates,y=val), colour="red") + coord_cartesian(xlim=c(as.Date("1981-01-01"), as.Date("2018-01-01")))





############ Check the goodness of fit
head(yup)
head(flow_border)




smsemoa(fitness.fun, n.objectives = NULL, n.dim = NULL, minimize = NULL,
lower = NULL, upper = NULL, mu = 100L, ref.point = NULL,
mutator = setup(mutPolynomial, eta = 25, p = 0.2, lower = lower, upper =
upper), recombinator = setup(recSBX, eta = 15, p = 0.7, lower = lower, upper
= upper), terminators = list(stopOnIters(100L)), .



require(ecr)

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}


optim(c(-1.2,1), fr, lower=c(-4, -5), upper=c(10,10), method="L-BFGS-B")


(res <- optim(c(-1.2,1), fr, grr, method = "BFGS"))


optimHess(res$par, fr, grr)
optim(c(-1.2,1), fr, NULL, method = "BFGS", hessian = TRUE)
## These do not converge in the default number of steps
optim(c(-1.2,1), fr, grr, method = "CG")
optim(c(-1.2,1), fr, grr, method = "CG", control = list(type = 2))
optim(c(-1.2,1), fr, grr, method = "L-BFGS-B")











qsim_df %>%
	group_by(date) %>%
	summarize(flow_atgauge = sum(flow_cfs, na.rm=TRUE), count = length(.))




MeanAnSolidPrecip

### Function to determine the percent of each day below freezing
perc_freeze_func <- function(temp_vals) {
	freeze_ecdf <- ecdf(quantile(temp_vals, probs=seq(0.05,0.95, 0.05)))
	perc_freeze <- freeze_ecdf(0)
	return(perc_freeze)
}
### Calculate percent of each day below freezing
prism_i$perc_freeze <- apply(select(prism_i, tmin_c, tavg_c, tmax_c), 1, perc_freeze_func)
prism_i <- prism_i %>%
	mutate(snow_mm = perc_freeze * prcp_mm)
	
	
	



#Drainage area: 172 square miles
wat_area_milesqr <- 172
wat_area_msqr <- wat_area_milesqr * 2589988.11

#convert cfs to mm/watershed/day
#Qls: outlet discharge [l/s]

#clim_join <- clim_gauge %>%
#	rowwise() %>% 
#	mutate(TAVG=mean(c(TMIN,TMAX), na.rm=TRUE)) %>%
#	select(DATE, PRCP, TAVG, TMIN, TMAX)

comb_df <- flow_border %>% 
	dplyr::select(dates, val) %>%
	inner_join(clim_gauge, by=c("dates" = "Date"))

comb_df <- comb_df %>%
	mutate(flow_ls = val * 28.316846592) %>%  ### convert to l/s
	mutate(flow_cuft_day = val * 86400) %>% ### Convert to cubic feet per day
	mutate(flow_cum_day = flow_cuft_day * .0283168 ) %>% ### convert to m3 per day
	mutate(flow_m_day = flow_cum_day/wat_area_msqr) %>%  ### Divide by watershed area in m2 to get m per day runoff
	mutate(flow_mm_day = flow_m_day * 1000) %>% ### Convert to mm per day
	dplyr::select(dates, tmin_c, tavg_c, tmax_c, prcp_mm, flow_ls, flow_mm_day) %>%
	rename("date" = "dates") %>%
	complete(date = seq.Date(min(date), max(date), by="day")) %>%
	arrange(date)


### Calculate PET, had to do this because PET doesn't allow NAs
pet_df <- comb_df %>%
	dplyr::select(date, tavg_c) %>%
	mutate(JD = as.POSIXlt(comb_df$date)$yday) %>%
	drop_na() %>%
	mutate(pet = PEdaily_Oudin(JD = JD, Temp = tavg_c,  Lat = 41.265, LatUnit = "deg")) %>%
	dplyr::select(date, pet)

### Merge back	
comb_df <- comb_df %>%
	full_join(pet_df, by="date")

Sys.setenv(TZ='UTC')
comb_df$date_pos <- as.POSIXct(comb_df$date, format="%Y-%m-%d")


### This is one way to fill gaps
#comb_df$tavg_c_fill <- na.approx(zoo(comb_df$tavg_c, comb_df$date), maxgap=3, na.rm=FALSE)
#comb_df$prcp_mm_fill <- na.approx(zoo(comb_df$prcp_mm, comb_df$date), maxgap=3, na.rm=FALSE)

### Can I do with missing data?
comb_df <- comb_df %>%
	drop_na()

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, DatesR = comb_df$date_pos, Precip = comb_df$prcp_mm, PotEvap = comb_df$pet, TempMean = comb_df$tavg_c,
TempMin = comb_df$tmin_c, TempMax = comb_df$tmax_c, ZInputs = z_clim, HypsoData = hypso_data, NLayers = 10)
str(InputsModel)

 
 
Ind_Run <- seq(which(format(comb_df$date_pos, format = "%d/%m/%Y")=="01/01/1981"), 
               which(format(comb_df$date_pos, format = "%d/%m/%Y")=="31/12/2016"))
str(Ind_Run)


## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel, IndPeriod_Run = Ind_Run)

## calibration criterion: preparation of the InputsCrit object 
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, RunOptions = RunOptions, Qobs = comb_df$flow_mm_day[Ind_Run])
#InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_RMSE, InputsModel = InputsModel, RunOptions = RunOptions, Qobs = comb_df$flow_mm_day[Ind_Run])
str(InputsCrit)


## preparation of CalibOptions object
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_CemaNeigeGR4J, FUN_CALIB = Calibration_Michel)

## calibration
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
InputsCrit = InputsCrit, CalibOptions = CalibOptions,
FUN_MOD = RunModel_CemaNeigeGR4J,
FUN_CRIT = ErrorCrit_NSE)

Param <- OutputsCalib$ParamFinalR
Param

OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)
str(OutputsModel)


## results preview
plot(OutputsModel, Qobs = comb_df$flow_mm_day[Ind_Run])

## efficiency criterion: Kling-Gupta Efficiency and NSE
OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE, InputsModel = InputsModel, RunOptions = RunOptions, comb_df$flow_mm_day[Ind_Run])
OutputsCrit <- ErrorCrit_KGE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)






## loading generalist parameter sets
data(Param_Sets_GR4J)

Param_Sets_GR4J$X4 <- Param_Sets_GR4J$X4u / 5.995 * (wat_area_msqr/1E6)^0.3
Param_Sets_GR4J$X4u <- NULL
Param_Sets_GR4J$X5 <- 0.95
Param_Sets_GR4J$X6 <- 15
Param_Sets_GR4J <- as.matrix(Param_Sets_GR4J)



## preparation of RunOptions object
### Completely made up annual solid precip
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel, IndPeriod_Run = Ind_Run, MeanAnSolidPrecip = c(30,40,50,60,150))


crit_val <- rep(NA, dim(Param_Sets_GR4J)[1])

for (i in seq(1,dim(Param_Sets_GR4J)[1])) {
	OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param_Sets_GR4J[i,])
	OutputsCrit <- ErrorCrit_RMSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)
	
	crit_val[i] <- OutputsCrit$CritValue
}

crit_val
which.min(crit_val)

i <- 14

OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param_Sets_GR4J[i,])
	

## results preview
plot(OutputsModel, Qobs = comb_df$flow_mm_day[Ind_Run])


param_manual <- c(220, 2.127, 125.116, 20, .01, 50)

OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_manual)
	
## results preview
plot(OutputsModel, Qobs = comb_df$flow_mm_day[Ind_Run])





