
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
	ZInputs = elev_gauge,
	HypsoData = hypso$elev_m,
	NLayers = n_layers)

### Create tibble of fraction solid precipitation
LayerFracSolidPrecip <- extrap$LayerFracSolidPrecip
names(LayerFracSolidPrecip) <- paste0("Layer_", seq(1,length(LayerFracSolidPrecip)))

LayerFracSolidPrecip <- as_tibble(LayerFracSolidPrecip) %>%
	dplyr::mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerFracSolidPrecip" , -DatesR)
	
### Create tibble of Mean Temp
LayerTempMean <- extrap$LayerTempMean
names(LayerTempMean) <- paste0("Layer_", seq(1,length(LayerTempMean)))

LayerTempMean <- as_tibble(LayerTempMean) %>%
	dplyr::mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerTempMean" , -DatesR)
	
### Create tibble of Precip
LayerPrecip <- extrap$LayerPrecip
names(LayerPrecip) <- paste0("Layer_", seq(1,length(LayerPrecip)))

LayerPrecip <- as_tibble(LayerPrecip) %>%
	dplyr::mutate(DatesR = prism_monthly$date_pos) %>%
	gather("Layer", "LayerPrecip" , -DatesR)

### Create a object with snow and rain for each of the layers
layer_clim <- LayerTempMean %>%
	full_join(LayerFracSolidPrecip, by=c("DatesR", "Layer") ) %>%
	full_join(LayerPrecip, by=c("DatesR", "Layer") ) 
	
layer_clim <- layer_clim %>%
	dplyr::mutate(layer_snow = LayerFracSolidPrecip * LayerPrecip) %>%
	dplyr::mutate(layer_rain = LayerPrecip - layer_snow)
	
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
	dplyr::mutate(potential_melt = (1-LayerFracSolidPrecip) * melt_index * 30 * LayerTempMean)### This would use the same melt rate each day that is above freezing  (2.74 mm/degree-day C is good assumption)
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
	dplyr::mutate(Layer = factor(Layer, levels=paste0("Layer_", seq(1,n_layers))))
	
snow_layer$layer_rain[snow_layer$layer_rain < 0] <- 0
snow_layer$melt[snow_layer$melt < 0] <- 0

snow_layer <- snow_layer %>%
	dplyr::mutate(rain_n_melt = layer_rain + melt)
	
### sum across layers and divide by n_layers
snow_total <- snow_layer %>%
	group_by(DatesR) %>%
	summarise(rain = sum(layer_rain)/n_layers, 
		snow = sum(layer_snow)/n_layers, 
		snowpack = sum(snowpack)/n_layers, 
		melt = sum(melt)/n_layers, 
		rain_n_melt = sum(rain_n_melt)/n_layers)

### Add proper dates back in
date_df <- prism_df %>%
	select(year, month, date, date_pos) %>%
	distinct()

snow_layer <- date_df %>%
	right_join(snow_layer, by=c("date_pos" = "DatesR")) %>%
	arrange(date)

snow_total <- date_df %>%
	right_join(snow_total, by=c("date_pos" = "DatesR")) %>%
	arrange(date)

return (list(snow_total = snow_total, by_layers = snow_layer))
}

