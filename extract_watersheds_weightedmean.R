


st_area()

x_and_y = st_intersection(hu_12, ssurgo_sma_df)


plot(hu_12["HUC12"])
plot(x_and_y["HUC12"])



nz_avheight = aggregate(x = ssurgo_sma_df["ksat_avg"], hu_12["HUC12"], FUN = mean, na.rm=TRUE)



x_and_y <- x_and_y %>%
	mutate(area = st_area(.))
	
subw <- x_and_y %>%
	filter(HUC12 == "160102010200")
	
weighted.mean(subw$ksat_avg, subw$area, na.rm=TRUE)

trier <- st_join(hu_12, nz_avheight, join=st_equals, left=TRUE)

trier %>%
 filter(HUC12 == "160102010200")


summ_byhuc <- x_and_y %>% 
	st_set_geometry(NULL) %>%
	group_by(HUC12) %>%
	summarize(mean_ksat = weighted.mean(ksat_avg, area, na.rm = TRUE))

trier <- hu_12 %>%
	left_join(summ_byhuc)
	
	
plot(trier["mean_ksat"])	
	
	