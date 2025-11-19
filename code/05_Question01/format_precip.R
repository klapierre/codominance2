#Interannual Precipitation addendum --------------------------------------
#Calculation of the Coefficient of Variation for Annual Precipitation at Sites
#Precipitation Data from MSWEP, modified code from Ingrid Slette, Feb 19, 2025
sites <- FinalAridity #

# make that file a SpatVector
site <- sites %>% vect(geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# list all of the monthly mswep precip data files
# change this to location to which you downloaded these files
r_paths <- list.files("MSWEP_Daily",
                      full.names = TRUE) %>%
  sort()

# make that a SpatRaster
r <- rast(r_paths)

# extract daily precip data for each site
ppt_daily <- terra::extract(r, site, bind = TRUE)

df <- as.data.frame(ppt_daily)

names(df) <- c("OID_", "site_code", "site_proj_comm", "MAP", "MAT", "GDiv", "ANPP", "NDep", "HumanFootprint", "Aridity", "lumpMode", "LumpNames", paste0("precip_", time(r)))

out <- pivot_longer(df, -c("OID_", "site_code", "site_proj_comm", "MAP", "MAT", "GDiv", "ANPP", "NDep", "HumanFootprint", "Aridity", "lumpMode", "LumpNames"), names_to = "date",
                    values_to = "precip") %>%
  mutate(date = str_replace(date, "^precip_", ""))

# create a new column for the year
out$year <- substr(out$date, 1, 4)

# create a new column for the month
out$month <- substr(out$date, 6, 7)

# create a new column for the day
out$day <- substr(out$date, 9, 10)

#calculate cv for annual precipitation at each site
interannual_precip <- filter(out, year != 2020 & year != 2025) %>% group_by(OID_,year) %>%
  summarise(annual_total = sum(precip, na.rm = T)) %>%
  summarise(
    mean_annual = mean(annual_total, na.rm = TRUE),
    sd_annual = sd(annual_total, na.rm = TRUE),
    cv_Precip = sd_annual / mean_annual)

Interannual.precip <- merge(FinalAridity, interannual_precip[,c(1,4)], by = "OID_" )
df_IAP <- as.data.frame(Interannual.precip)
write.csv(Interannual.precip, file = "C:/Users/msgrabda/Downloads/IAP.csv")
