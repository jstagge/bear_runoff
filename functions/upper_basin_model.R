
###########################################################################
## Upper Basin Model Function
###########################################################################

upper_basin_model <- function(x, huc_info, clim_data, hypso) {

param_c <- x

param_names <- c("X1", "X2")#, "melt_ind")
param_huc <- huc_info$huc
huc_n <- length(param_huc)


model_param <- matrix(param_c, length(param_huc), length(param_names))
rownames(model_param) <- param_huc
colnames(model_param) <- param_names

model_param <- as_tibble(model_param) %>%
	dplyr::mutate(huc = param_huc)
	
for(i in seq(1,huc_n)){

huc_i <- param_huc[i]

### Cut Prism Climate data to HUC
prism_i <- clim_data %>%
	filter(huc == huc_i)
	
prism_monthly <- prism_i %>%
	dplyr::mutate(year_month = paste0(year(date), "-",month(date))) %>%
	dplyr::mutate(month = month(date), year = year(date)) %>%
	dplyr::mutate(date_pos = as.POSIXct(date, tz = "America/Denver")) %>%
	group_by(year, month) %>%
	summarise(date = min(date), prcp_mm = sum(prcp_mm, na.rm=TRUE), tmin_c = mean(tmin_c, na.rm=TRUE), tavg_c = mean(tavg_c, na.rm=TRUE), tmax_c = mean(tmax_c, na.rm=TRUE), pet=sum(pet, na.rm=TRUE), date_pos = min(date_pos)) %>%
	arrange(date)

### Cut hypso data to HUC
hypso_i <- hypso %>%
	filter(huc == huc_i) %>%
	arrange(quant)

### Cut z for gauge
z_clim <- huc_info$elev_m[i]

### Add in Snow Module here
snow_results <- snow_module(prism_monthly = prism_monthly , hypso = hypso_i , elev_gauge = z_clim,  melt_index = 2.74, n_layers = 10)

### Join back with prism data
snow_total <- snow_results$snow_total
snow_total <- snow_total %>%
	select(-year, -month, -date_pos)

prism_monthly <- prism_monthly %>%
	left_join(snow_total, by="date")

### Create input model
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M, 
	DatesR = prism_monthly$date_pos, 
	Precip = prism_monthly$rain_n_melt, 
	PotEvap = prism_monthly$pet)

### Create index period for run 
Ind_Warm <- seq(which(prism_monthly$date == "1981-01-01"), 
               which(prism_monthly$date == "1981-12-01"))

Ind_Run <- seq(which(prism_monthly$date == "1982-01-01"), 
               which(prism_monthly$date == "2016-12-01"))

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
#as.Date(OutputsModel$DatesR)
qsim_temp <- tibble(date = prism_monthly$date[Ind_Run] , huc = huc_i, qsim_mm = OutputsModel$Qsim) %>%
	dplyr::mutate(qsim_m_month = qsim_mm / 1000) %>% ### Convert to m per month
	dplyr::mutate(qsim_cum_month = qsim_m_month * sub_area_msqr) %>%  ### Multiply by watershed area in m2 to get m3 per month runoff
	dplyr::mutate(flow_cuft_month = qsim_cum_month / .0283168 ) %>% ### convert to ft3 per month
	dplyr::mutate(flow_cuft_day = flow_cuft_month / 30) %>% ### Convert to ft3 per day
	dplyr::mutate(flow_cfs = flow_cuft_day / 86400)  ### Convert to cubic feet per day

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

return(pred_df)

}






###########################################################################
## Optim function
###########################################################################
#huc_info, clim_data, hypso, obs

optim_upper_basin <- function(x, ...) {
require(hydroGOF)

huc_info <- huc_info
clim_data <- clim_data
hypso <- hypso
obs <- obs

### Run upper basin model
pred_df <- upper_basin_model(x=x, huc_info = huc_info, clim_data = clim_data, hypso = hypso)

### Combine to get the gauge 
pred_df$at_gauge <- apply(pred_df[,2:6],1,sum)

comb_df <- pred_df %>% 
	select(date, at_gauge) %>%
	left_join(obs, by="date") %>%
#	select(-staid) %>%
	rename("pred" = "at_gauge", "obs" = "val")
	
comb_df <- comb_df  %>%
	dplyr::mutate(month = month(date), year = year(date))  %>%
	dplyr::mutate(resid = obs - pred) %>%
	arrange(date)

	
#gof_run <- gof(comb_df$pred, comb_df$obs, na.rm=TRUE)
gof_run <- gof(comb_df$pred, comb_df$obs, na.rm=TRUE)

### NSE
#gof_run[[9]]

### Find annual min and max, test the RMSE for each
extremes <- comb_df %>%
	dplyr::mutate(year = year(date)) %>%
	group_by(year) %>%
	summarise(obs_min = min(obs), pred_min = min(pred), obs_max = max(obs), pred_max = max(pred))

gof_min <- gof(extremes$pred_min, extremes$obs_min, na.rm=TRUE)
gof_max <- gof(extremes$pred_max, extremes$obs_max, na.rm=TRUE)

return_results <- c(gof_run[[9]], gof_min[[4]], gof_max[[4]])
names(return_results) <- c("NSE", "RMSE_ann_min", "RMSE_ann_max")

return(return_results)
}