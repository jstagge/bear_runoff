

 # load waterData package, assuming it has already been installed on the system
require(waterData)
require(airGR)
require(tidyverse)
#require(colorout)


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


### From PRISM
clim_gauge <- read_csv("PRISM_ppt_tmin_tmean_tmax_stable_4km_19810101_20170101_40.7949_-110.7903.csv", skip=10, col_names = TRUE)

clim_gauge <- clim_gauge %>%
	rename(prcp_mm = `ppt (mm)`) %>%
	rename(tmin_c = `tmin (degrees C)`) %>%
	rename(tavg_c = `tmean (degrees C)`) %>%
	rename(tmax_c = `tmax (degrees C)`)
	
head(clim_gauge)
z_clim <- 2996



####################### Read in and process elevation for watershed
########################THIS NEEDS TO BE FIXED
elev_upper_bear <- as.matrix(read.table("elev_10011500.txt", header=FALSE, skip=6, na.strings = "-9999", as.is=TRUE))
plot(raster(elev_upper_bear))


hypso_quant <- quantile(elev_upper_bear, seq(0.01,0.99, 0.01), na.rm=TRUE) 
hypso_data <- c(min(elev_upper_bear, na.rm=TRUE), hypso_quant, max(elev_upper_bear, na.rm=TRUE))
hypso_data <- hypso_data/100 ### Convert to m

plot(hypso_data, type="l")
abline(h=z_clim,col="red")

rm(elev_upper_bear, hypso_quant)


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





