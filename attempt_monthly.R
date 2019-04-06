
flow_monthly <- flow_border %>%
	mutate(year_month = paste0(year(date), "-",month(date))) %>%
	group_by(year_month) %>%
	summarise(val = mean(val, na.rm=TRUE), date = min(date)) %>%
	arrange(date)


fds <- left_join(pred_df, flow_monthly)



prism_monthly <- prism_i %>%
	mutate(year_month = paste0(year(date), "-",month(date))) %>%
	group_by(year_month) %>%
	summarise(prcp_mm = sum(prcp_mm, na.rm=TRUE), tavg_c = mean(tavg_c, na.rm=TRUE), pet=sum(pet, na.rm=TRUE), date_pos = min(date_pos)) %>%
	arrange(date_pos)

##prism_monthly <- prism_monthly %>%
##	mutate(prcp_mm = prcp_mm * 30, pet = pet * 30)
  
   
  
### Create input model
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M, 
	DatesR = prism_monthly$date_pos, 
	Precip = prism_monthly$prcp_mm, 
	PotEvap = prism_monthly$pet)

 ### Create index period for run 
Ind_Run <- seq(which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/1982"), 
               which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/2017"))
			   
Ind_Warm <- seq(which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/1981"), 
               which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/12/1981"))

### Cut model parameters to HUC 
param_i <- model_param  %>%
	filter(huc == huc_i) %>%
	select(-huc) %>%
	select(X1, X2) %>%
	unlist() 
	
## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M, InputsModel = InputsModel, IndPeriod_Run = Ind_Run, IndPeriod_WarmUp = Ind_Warm)

## Calculate Outputs
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i)










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


###########################
###  Model for the upper basin - needs to be a function to run optimization
###########################

upper_basin_model <- function(x, huc_info) {

param_c <- x

param_names <- c("X1", "X2")#, "melt_ind")
param_huc <- huc_info$huc
huc_n <- length(param_huc)


model_param <- matrix(param_c, length(param_huc), length(param_names))
rownames(model_param) <- param_huc
colnames(model_param) <- param_names

model_param <- as_tibble(model_param) %>%
	mutate(huc = param_huc)
	
for(i in seq(1,huc_n)){

huc_i <- param_huc[i]

### Cut Prism Climate data to HUC
prism_i <- prism_df %>%
	filter(huc == huc_i)
	
prism_monthly <- prism_i %>%
	mutate(year_month = paste0(year(date), "-",month(date))) %>%
	group_by(year_month) %>%
	summarise(prcp_mm = sum(prcp_mm, na.rm=TRUE), tmin_c = mean(tmin_c, na.rm=TRUE), tavg_c = mean(tavg_c, na.rm=TRUE), tmax_c = mean(tmax_c, na.rm=TRUE), pet=sum(pet, na.rm=TRUE), date_pos = min(date_pos)) %>%
	arrange(date_pos)

### Cut hypso data to HUC
hypso_i <- hypso_quant %>%
	filter(huc == huc_i) %>%
	arrange(quant)

### Cut z for gauge
z_clim <- huc_info$elev_m[i]

### Add in Snow Module here
snow_results <- snow_module(prism_monthly = prism_monthly , hypso = hypso_i , elev_gauge = z_clim,  melt_index = 2.74, n_layers = 10)

### Join back with prism data
prism_monthly <- prism_monthly %>%
	left_join(snow_results$snow_total, by=c("date_pos" = "DatesR"))

### Create input model
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M, 
	DatesR = prism_monthly$date_pos, 
	Precip = prism_monthly$rain_n_melt, 
	PotEvap = prism_monthly$pet)

### Create index period for run 
Ind_Run <- seq(which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/1982"), 
               which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/2017"))
			   
Ind_Warm <- seq(which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/01/1981"), 
               which(format(prism_monthly$date_pos, format = "%d/%m/%Y")=="01/12/1981"))

### Cut model parameters to HUC 
param_i <- model_param  %>%
	filter(huc == huc_i) %>%
	select(-huc) %>%
	select(X1, X2) %>%
	unlist() 
	
## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M, InputsModel = InputsModel, IndPeriod_Run = Ind_Run, IndPeriod_WarmUp = Ind_Warm)

## Calculate Outputs
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i)

### Extract area of subbasin
sub_area_msqr <- hu_char$area_m_sqr[i]

### Calculate flow and convert units
qsim_temp <- tibble(date = as.Date(OutputsModel$DatesR), huc = huc_i, qsim_mm = OutputsModel$Qsim) %>%
	mutate(qsim_m_month = qsim_mm / 1000) %>% ### Convert to m per month
	mutate(qsim_cum_month = qsim_m_month * sub_area_msqr) %>%  ### Multiply by watershed area in m2 to get m3 per month runoff
	mutate(flow_cuft_month = qsim_cum_month / .0283168 ) %>% ### convert to ft3 per month
	mutate(flow_cuft_day = flow_cuft_month / 30) %>% ### Convert to ft3 per day
	mutate(flow_cfs = flow_cuft_day / 86400)  ### Convert to cubic feet per day

if (i == 1) {
	qsim_df <- qsim_temp
} else {
	qsim_df <- rbind(qsim_df, qsim_temp)
}

}

### Spread this over HUC regions
pred_df <- qsim_df %>%
	select(date, huc, flow_cfs) %>%
	spread(huc, flow_cfs)

pred_df$at_gauge <- apply(pred_df[,2:6],1,sum)

comb_df <- pred_df %>% 
	select(date, at_gauge) %>%
	left_join(flow_monthly, by="date") %>%
#	select(-staid) %>%
	rename("pred" = "at_gauge", "obs" = "val")
	
comb_df <- comb_df %>%
	mutate(resid = obs - pred)
	
#gof_run <- gof(comb_df$pred, comb_df$obs, na.rm=TRUE)
gof_run <- gof(comb_df$pred, comb_df$obs, na.rm=TRUE)

return(- gof_run[[2]])

}

#upper_basin_model(x=param_c, huc_info = huc_char)


#yup <- optim(x=x, upper_basin_model)#, huc_info=hu_char)

require(GA)
#require(snow)
lower_bounds <- c(600, 0.5)#, 1)
lower_bounds <- rep(lower_bounds, each=5)
upper_bounds <- c(2000, 1.2)#, 100)
upper_bounds <- rep(upper_bounds, each=5)

yup <- ga(type = "real-valued", fitness = upper_basin_model, lower = lower_bounds, upper = upper_bounds, huc_info=hu_char, maxiter = 50) #, parallel=TRUE)


plot(yup)
yup2 <- summary(yup)
x <- yup2$solution[1,]

param_c <- x

param_names <- c("X1", "X2")#, "melt_ind")
param_huc <- huc_info$huc
huc_n <- length(param_huc)


model_param <- matrix(param_c, length(param_huc), length(param_names))
rownames(model_param) <- param_huc
colnames(model_param) <- param_names

Best = 44.510

 x
          x1           x2           x3           x4           x5           x6           x7           x8           x9          x10 
 758.8100691  695.1356203  839.5295322  662.9173023 1165.0901269    0.7253586    1.0044219    0.5580676    0.8486673    1.1407137 
> 

> model_param
                    X1        X2
160101010103  758.8101 0.7253586
160101010104  695.1356 1.0044219
160101010101  839.5295 0.5580676
160101010102  662.9173 0.8486673
160101010106 1165.0901 1.1407137




comb_df <- comb_df %>%
	mutate(month = month(date))
	
plot_df <- comb_df %>%
	select(date, obs, pred) %>%
	gather("variable", "flow_cfs", -date)

p <- ggplot(plot_df, aes(x=date, y=flow_cfs, colour=variable)) %>%
	+ geom_line(size=0.3) %>%
	+ theme_classic(9) %>%
	+ scale_x_date(name = "Date", date_breaks = "5 years", date_labels = "%Y") %>%
	+ scale_y_continuous(name="Monthly Mean Flow (cfs)") %>%
	+ scale_colour_manual(name = "", values = c("Black", "Red"), labels=c("Observed", "Predicted")) %>%
	+ theme(legend.position="bottom")

p

ggsave("fit_ts.png", p, width=10.5, height=4, dpi=320)

p <- ggplot(comb_df, aes(x=obs, y=pred)) %>%
	+ geom_abline(intercept = 0, slope = 1, color="grey50", size=0.5) %>%
	+ geom_point() %>%
	+ theme_classic(9) %>%
	+ scale_x_continuous(name = "Observed (cfs)") %>%
	+ scale_y_continuous(name="Predicted (cfs)") %>%
	+ coord_fixed(ratio = 1)

p


ggsave("fit_pred_obs.png", p, width=4.5, height=4.5, dpi=320)

p <- ggplot(comb_df, aes(x=obs, y=pred)) %>%
	+ geom_abline(intercept = 0, slope = 1, color="grey50", size=0.5) %>%
	+ geom_point() %>%
	+ theme_bw(9) %>%
	+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) %>%
	+ scale_x_continuous(name = "Observed (cfs)") %>%
	+ scale_y_continuous(name="Predicted (cfs)") %>%
	+ facet_wrap(~month, scales = "free")


ggsave("fit_pred_obs_bymonth.png", p, width=10, height=10, dpi=320)


p <- ggplot(comb_df, aes(x=obs, y=resid)) %>%
	+ geom_abline(intercept = 0, slope = 0, color="grey50", size=0.5) %>%
	+ geom_point() %>%
	+ theme_classic(9) %>%
	+ scale_x_continuous(name = "Observed (cfs)") %>%
	+ scale_y_continuous(name="Residuals (cfs)") %>%
	+ coord_fixed(ratio = 1)

p

ggsave("fit_resid_obs.png", p, width=4.5, height=4.5, dpi=320)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	