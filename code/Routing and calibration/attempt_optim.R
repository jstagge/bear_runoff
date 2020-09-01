require(hydroGOF)

mutate <- dplyr::mutate



param_c <- c(rep(220,huc_n), rep(2.127,huc_n), rep(125.116,huc_n), rep(20,huc_n), rep(0.01,huc_n), rep(50,huc_n))



###########################################################################
## Try to set up a run
###########################################################################


upper_basin_model <- function(x, huc_info) {

param_c <- x

param_names <- paste0("X", seq(1,6))
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

### Cut hypso data to HUC
hypso_i <- hypso_quant %>%
	filter(huc == huc_i) %>%
	arrange(quant)

### Cut z for gauge
z_clim <- huc_info$elev_m[i]

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

 ### Create index period for run 
Ind_Run <- seq(which(format(prism_i$date_pos, format = "%d/%m/%Y")=="01/01/1981"), 
               which(format(prism_i$date_pos, format = "%d/%m/%Y")=="31/12/2016"))

### Cut model parameters to HUC 
param_i <- model_param  %>%
	filter(huc == huc_i) %>%
	select(-huc) %>%
	unlist()

## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel, IndPeriod_Run = Ind_Run)

## Calculate Outputs
OutputsModel <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = param_i)

### Extract area of subbasin
sub_area_msqr <- hu_char$area_m_sqr[i]

### Calculate flow and convert units
qsim_temp <- tibble(date = as.Date(OutputsModel$DatesR), huc = huc_i, qsim_mm = OutputsModel$Qsim) %>%
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

### Spread this over HUC regions
pred_df <- qsim_df %>%
	select(date, huc, flow_cfs) %>%
	spread(huc, flow_cfs)

pred_df$at_gauge <- apply(pred_df[,2:6],1,sum)

comb_df <- pred_df %>% 
	select(date, at_gauge) %>%
	left_join(flow_border) %>%
	select(-staid) %>%
	rename("pred" = "at_gauge", "obs" = "val")
	
comb_df <- comb_df %>%
	mutate(resid = obs - pred)
	
### Calculate monthly mean
comb_df <- comb_df %>% 
	mutate(year_month = paste0(year(date), "-",month(date)))
	
summ_df <- comb_df %>%	
	group_by(year_month) %>%
	summarise(obs_mean = mean(obs, na.rm=TRUE), pred_mean = mean(pred, na.rm=TRUE))
  
#gof_run <- gof(comb_df$pred, comb_df$obs, na.rm=TRUE)
gof_run <- gof(summ_df$pred_mean, summ_df$obs_mean, na.rm=TRUE)

return(- gof_run[[2]])

}

#upper_basin_model(x=param_c, huc_info = huc_char)


#yup <- optim(x=param_c, upper_basin_model)#, huc_info=hu_char)

require(GA)
require(snow)
lower_bounds <- c(0.0001, -20, 0.0001, 0.5, 0, 1 )
lower_bounds <- rep(lower_bounds, each=5)
upper_bounds <- c(3000, 20, 2000, 15, 1, 50)
upper_bounds <- rep(upper_bounds, each=5)

yup <- ga(type = "real-valued", fitness = upper_basin_model, lower = lower_bounds, upper = upper_bounds, huc_info=hu_char, maxiter = 30) #, parallel=TRUE)




require(ecr)

res = ecr(fitness.fun = upper_basin_model, representation = "float",
  n.dim = 30, survival.strategy = "plus")
  
  
  ,
  lower = lower, upper = upper,
  mu = MU, lambda = LAMBDA,
  mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
  terminators = list(stopOnIters(MAX.ITER)))





### Net Needed yet
### Route gauges
pred_df <- pred_df %>%
	mutate()
	
inflow <- pred_df %>% 
	select( "160101010101") %>% 
	unlist()
inflow <- inflow / 35.3147

begin_date <- as.character(min(pred_df$date))
end_date <- as.character(max(pred_df$date))

reachRouting(inflow,"muskingumcunge",list(bedWith=10,sideSlope=2,channelSlope=0.0156, manningRoughness=0.04, riverLength=15.1),list(start=begin_date,end=end_date,by=86400))








p <- ggplot(pred_df, aes(x=date, y=at_gauge)) + geom_line() + theme_classic()
p + geom_line(data = flow_border, aes(x=date,y=val), colour="red") + coord_cartesian(xlim=c(as.Date("1981-01-01"), as.Date("2018-01-01")))

p + geom_line(data = flow_border, aes(x=date,y=val), colour="red") + coord_cartesian(xlim=c(as.Date("1981-01-01"), as.Date("1985-01-01")))









inflow<-c(100,500,1500,2500,5000,11000,22000,28000,28500,26000,
22000,17500,14000,10000,7000,4500,2500,1500,1000,500,100)


inflow <- pred_df %>% select( "160101010101") %>% unlist()


routingMethod<-c("muskingum","muskingumcunge")

routingParams<-list(k=3,x=0.2,bedWith=50,sideSlope=2,channelSlope=0.0001,manningRoughness=0.01,riverLength=100)

simulation<-list(start='2000-01-01',end='2000-01-04',by=3600)''



inflow<-pred_df %>% select( "160101010101") %>% unlist()

inflow <- inflow / 35.3147

begin_date <- min(pred_df$date)
end_date <- max(pred_df$date)


routingMethod <-c("muskingum","muskingumcunge")

#k and x for muskingum
#bedwith [m] side slope[m/m] channel slope[m/m] manningrough and length [km] for muskingumcunge

routingParams<-list(k=3,x=0.2,bedWith=20,sideSlope=2,channelSlope=0.003, manningRoughness=0.04, riverLength=25)

simulation<-list(start="1981-01-01",end="2016-12-31",by=86400)

route1 <- reachRouting(inflow,routingMethod[1],routingParams,simulation)

route2 <- reachRouting(inflow,"muskingumcunge",list(k=3,x=0.2,bedWith=20,sideSlope=2,channelSlope=0.003, manningRoughness=0.04, riverLength=25),list(start="1981-01-01",end="2016-12-31",by=86400))

inflow <- inflow * 35.3147
outflow <- route2$operation$"23 Km" * 35.3147

#convert cms to cfs then back

plot(inflow[1:100], type="l")
lines(outflow[1:100], col="red")

plot(inflow[7000:8000], type="l")
lines(outflow[7000:8000], col="red")


inflow <- pred_df %>% 
	select( "160101010101") %>% 
	unlist()
inflow <- inflow / 35.3147

begin_date <- as.character(min(pred_df$date))
end_date <- as.character(max(pred_df$date))

reachRouting(inflow,"muskingumcunge",list(bedWith=10,sideSlope=2,channelSlope=0.0156, manningRoughness=0.04, riverLength=15.1),list(start=begin_date,end=end_date,by=86400))

Reach	Length (m)	Up elev (cm)	Down elev (cm)	Width (m)
Hayden Fork	15067	282975	259426	10
Stillwater Fork	18048	292606	259426	10
East Fork	17567	283559	253382	10
Bear 1	6413	259426	253382	10
Bear 2	6924	253382	243003	10
West Fork	18680	282406	243003	10





createBasin
createReach

createDiversion
createJunction(name, downstream, inflow, delayInflow)
