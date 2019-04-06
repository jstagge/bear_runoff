



Pretty significant tear
Pretty big gap between tear
Good indication for doing the surgery

crane, take hour or less
splint  immediately after surgery
followup 2 weeks later to remove stitches, then back in boot

knee walker pick up any time 
until 5 pm




Clear
9.8°C per 1,000 meters

Cloud
6°C per 1,000 meters

3800- 2600

So, almost 8 C difference in one huc




### Cut z for gauge
z_clim <- huc_info$elev_m[i]

### Create input model
extrap <- DataAltiExtrapolation_Valery(DatesR = prism_monthly$date_pos, 
	PrecipScale = TRUE,
	TempMean = prism_monthly$tavg_c,
	Precip = prism_monthly$prcp_mm,
	ZInputs = z_clim,
	HypsoData = hypso_i$elev_m,
	NLayers = 5)


plot(extrap$LayerTempMean[[1]][1:15], type="l")
lines(extrap$LayerTempMean[[5]][1:15], col="red")

plot(extrap$LayerFracSolidPrecip[[1]][1:15], type="l")
lines(extrap$LayerFracSolidPrecip[[5]][1:15], col="red")

plot(extrap$LayerPrecip[[1]][1:15], type="l")
lines(extrap$LayerPrecip[[5]][1:15], col="red")

$LayerFracSolidPrecip[[5]]


$LayerPrecip[[1]]




2567

zlayers
2502.67 2568.59 2599.20 2645.03 2731.83



huc_i <- param_huc[i]

### Cut Prism Climate data to HUC
prism_i <- prism_df %>%
	filter(huc == huc_i)
	
prism_monthly <- prism_i %>%
	mutate(year_month = paste0(year(date), "-",month(date))) %>%
	group_by(year_month) %>%
	summarise(prcp_mm = sum(prcp_mm, na.rm=TRUE), tmin_c = mean(tmin_c, na.rm=TRUE), tavg_c = mean(tavg_c, na.rm=TRUE), tmax_c = mean(tmax_c, na.rm=TRUE), date_pos = min(date_pos)) %>%
	arrange(date_pos)


### Cut hypso data to HUC
hypso_i <- hypso_quant %>%
	filter(huc == huc_i) %>%
	arrange(quant)

### Cut z for gauge
z_clim <- huc_info$elev_m[i]

###########################
###  Snow Module
###########################
snow_module <- function(prism_monthly , hypso , elev_gauge, melt_index = 2.74, n_layers = 5){

### Extrapolate climate data based on elevation
extrap <- DataAltiExtrapolation_Valery(DatesR = prism_monthly$date_pos, 
	PrecipScale = TRUE,
	TempMin = prism_monthly$tmin_c,
	TempMean = prism_monthly$tavg_c,
	TempMax = prism_monthly$tmax_c,
	Precip = prism_monthly$prcp_mm,
	ZInputs = z_clim,
	HypsoData = hypso_i$elev_m,
	NLayers = n_layers)

### Create tibble of fraction solid precipitation
LayerFracSolidPrecip <- extrap$LayerFracSolidPrecip
names(LayerFracSolidPrecip) <- paste0("Layer_", seq(1,length(LayerFracSolidPrecip)))

LayerFracSolidPrecip <- as_tibble(LayerFracSolidPrecip) %>%
	mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerFracSolidPrecip" , -DatesR)
	
### Create tibble of Mean Temp
LayerTempMean <- extrap$LayerTempMean
names(LayerTempMean) <- paste0("Layer_", seq(1,length(LayerTempMean)))

LayerTempMean <- as_tibble(LayerTempMean) %>%
	mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerTempMean" , -DatesR)
	
### Create tibble of Precip
LayerPrecip <- extrap$LayerPrecip
names(LayerPrecip) <- paste0("Layer_", seq(1,length(LayerPrecip)))

LayerPrecip <- as_tibble(LayerPrecip) %>%
	mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerPrecip" , -DatesR)

### Create a object with snow and rain for each of the layers
layer_clim <- LayerTempMean %>%
	full_join(LayerFracSolidPrecip, by=c("DatesR", "Layer") ) %>%
	full_join(LayerPrecip, by=c("DatesR", "Layer") ) 
	
layer_clim <- layer_clim %>%
	mutate(layer_snow = LayerFracSolidPrecip * LayerPrecip) %>%
	mutate(layer_rain = LayerPrecip - layer_snow)
	
### Melt algorithm
for (k in unique(layer_clim$Layer)) {

layer_clim_k <- layer_clim %>%
	filter(Layer == k)

length_months <-  dim(layer_clim_k)[1]
 
snow_pack <- rep(NA, length_months)
snow_pack[1] <- layer_clim_k$layer_snow[1]
melt <- rep(NA, length_months)
melt[1] <- 0

layer_clim_k <- layer_clim_k %>%
	mutate(potential_melt = (1-LayerFracSolidPrecip) * melt_index * 30 * LayerTempMean)### This would use the same melt rate each day that is above freezing  (2.74 mm/degree-day C is good assumption)
	### Still isnt quite right

	for (j in seq(2,length_months)) {
		snow_pack_j <- snow_pack[j-1] + layer_clim_k$layer_snow[j]
		melt_j <- layer_clim_k$potential_melt[j]
		
		if(melt_j > snow_pack_j){
			melt_j <- snow_pack_j
		}
		snow_pack_j <- snow_pack_j - melt_j
		snow_pack[j] <- snow_pack_j
		melt[j] <- melt_j
	}

layer_clim_k$melt <- melt
layer_clim_k$snowpack <- snow_pack

if (k == "Layer_1" ) {
	snow_layer <- layer_clim_k
} else {
	snow_layer <- rbind(snow_layer, layer_clim_k)
}
}

snow_layer <- snow_layer %>% 
	mutate(Layer = factor(Layer, levels=paste0("Layer_", seq(1,n_layers))))
	
snow_layer$layer_rain[snow_layer$layer_rain < 0] <- 0
snow_layer$melt[snow_layer$melt < 0] <- 0

snow_layer <- snow_layer %>%
	mutate(rain_n_melt = layer_rain + melt)
	
### sum across layers and divide by n_layers
snow_total <- snow_layer %>%
	group_by(DatesR) %>%
	summarise(rain = sum(layer_rain)/n_layers, 
		snow = sum(layer_snow)/n_layers, 
		snowpack = sum(snowpack)/n_layers, 
		melt = sum(melt)/n_layers, 
		rain_n_melt = sum(rain_n_melt)/n_layers)
	
return (list(snow_total = snow_total, by_layers = snow_layer))
}



huzzah <- snow_module(prism_monthly = prism_monthly , hypso = hypso_i , elev_gauge = z_clim,  n_layers = 20)




### Spread across columns
yup_doodle <- snow_layer %>%
	select(DatesR, Layer, rain_n_melt) %>%
	spread(Layer, rain_n_melt)




ggplot(filter(huzzah$by_layers, DatesR < as.Date("1985-01-01")), aes(x=DatesR, y=snowpack, colour=Layer)) + geom_line() + scale_colour_viridis(discrete=TRUE)

ggplot(filter(snow_layer, DatesR < as.Date("1985-01-01")), aes(x=DatesR, y=rain_n_melt, colour=Layer)) + geom_line() + scale_colour_viridis(discrete=TRUE)



layer_clim <- layer_clim %>%
	left_join(layer_clim_k, by=c("DatesR", "Layer"))


DataAltiExtrapolation_Valery(DatesR, Precip, PrecipScale = TRUE,
TempMean, TempMin = NULL, TempMax = NULL,
ZInputs, HypsoData, NLayers, verbose = TRUE)



DatesR [POSIXt] vector of dates
Precip [numeric] time series of daily total precipitation (catchment average) [mm/d]
PrecipScale (optional) [boolean] indicating if the mean of the precipitation interpolated on
the elevation layers must be kept or not, required to create CemaNeige module
inputs, default = TRUE (the mean of the precipitation is kept to the original value)
TempMean [numeric] time series of daily mean air temperature [°C]
TempMin (optional) [numeric] time series of daily min air temperature [°C]
TempMax (optional) [numeric] time series of daily max air temperature [°C]
ZInputs [numeric] real giving the mean elevation of the Precip and Temp series (before
extrapolation) [m]
HypsoData [numeric] vector of 101 reals: min, q01 to q99 and max of catchment elevation
distribution [m]
NLayers [numeric] integer giving the number of elevation layers requested [-]



