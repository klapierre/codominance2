
################################################################################
##  05_Q1_site_modes.R: Calculate mode of codominance numbers in control plots and 
##  compare to environmental data.
##
##  Authors: Ashley LaRoque, Jordan Winter, Elise Grabda, Alyssa Young, Rachael Brenneman (modified K. Komatsu)
##  Date created: 
################################################################################

source("code/01_library.R")
source("code/02_functions.R")


# Read data ---------------------------------------------------------------
 
numCodomPlotYear <- readRDS("data/numCodomPlotYear.rds") %>% 
  separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 
                            'plot_id', 'treatment', 'calendar_year'), 
           sep='::', remove=F) %>% 
  left_join(readRDS("data/expInfo.rds"))

 
# Calculate mode across years for all plots ----------------------------------------------------------
 
# create a dataframe of codominant group for plots with a single timepoint (because can't calculate mode of singleton)
singletonCodomPlotYear <- numCodomPlotYear %>% 
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   mutate(length=length(plot_id)) %>% 
   ungroup() %>% 
   filter(length==1) %>% 
   rename(plot_codom=num_group) %>% 
   dplyr::select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, plot_codom) 
 
# calculate mode across years for all plots
modePlotTrue <- numCodomPlotYear %>%  
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   reframe(plot_codom = DescTools::Mode(num_group)) %>% # mode function must be capital here 
   ungroup() %>% 
   filter(!is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom = round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
   rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlotYear %>% 
   select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
   unique() %>% 
   full_join(modePlotTrue) %>% 
   filter(is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom = round(mean(num_group), digits=0), .groups='drop')
 
# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length = length(community_type)) %>% 
  ungroup() %>% 
  filter(length == 1) %>% 
  rename(mode_site = plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, mode_site) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
   filter(trt_type == 'control') %>%
   group_by(database, site_code, project_name, community_type) %>% # mode generated from these
   reframe(mode_site = DescTools::Mode(plot_codom)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type) %>% 
   summarise(mode_site = round(mean(mode_site), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_site)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  filter(trt_type == 'control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_codom) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(mode_site)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site = round(mean(plot_codom), digits=0), .groups='drop')

#############################################
## envData either needs to be changed or removed 
# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lumpMode = ifelse(mode_site == 3, 2, mode_site)) 

# how many sites?
siteCount <- modeSite %>%
  mutate(lat = round(Latitude, 1),
         long = round(Longitude, 1)) %>% 
  select(lat, long) %>% 
  unique()



# Fig: Distribution of values ---------------------------------------------


######################################################################
# within this chunk of '#', uses the df_aridity csv file 
## compared to the histogram made with 'df_hist', this hist has a much higher count of even

# Ashley's upload of aridity data: already with mode (mono, co, tri, even) and lumped mode (mono, co, even)
df_aridity <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Grad School/Terui Lab/Codominance/FinalAridity.csv") %>% 
  mutate(LumpNames = factor(LumpNames, c("Monodominated", "Codominated", "Even")))


df_arid <- df_aridity %>% 
  dplyr::select(lumpMode, LumpNames, lumpMode, MAP, MAT, GDiv, ANPP, HumanFootprint, NDep, Aridity)

df_h <- df_arid %>% 
  pivot_longer(cols = c("MAP", "MAT", "GDiv", "ANPP", "HumanFootprint", "NDep", "Aridity"), 
               names_to = "variable", values_to = "value") %>% 
  mutate(variable1 = case_when(variable == "MAP" ~ "MAP",
                               variable == "MAT" ~ "MAT",
                               variable == "GDiv" ~ "Gamma Diversity",
                               variable == "ANPP" ~ "ANPP",
                               variable == "HumanFootprint" ~ "Human Footprint Index",
                               variable == "NDep" ~ "N Deposition",
                               variable == "Aridity" ~ "Aridity")) 

# plot histogram 
ggplot(df_h, aes(value)) + # df_hist comes from formatted df above 
  geom_histogram(aes(fill = LumpNames), bins = 50) +
  facet_wrap(~ variable1,
             scales = "free") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(name = "",
                    values = c("#007BA7", "#A63922", "#D8B573")) +
  labs(x = "Value",
       y = "Count")


###################################################################

## format data
# df_wide <- modeSite %>% 
#   dplyr::select(mode_site, lumpMode, MAP, MAT, gamma_rich, anpp, HumanDisturbance, N_Deposition)
# 
# df_hist <- df_wide %>% 
#   pivot_longer(cols = c("MAP", "MAT", "gamma_rich", "anpp", "HumanDisturbance", "N_Deposition"), 
#                names_to = "variable", values_to = "value") %>% 
#   mutate(lumpMode = ifelse(lumpMode == 1, "Monodominated", 
#                            ifelse(lumpMode == 2, "Codominated", "Even")),
#          mode.levels = factor(lumpMode, levels = c("Monodominated", "Codominated", "Even")),
#          variable1 = case_when(variable == "MAP" ~ "MAP",
#                                variable == "MAT" ~ "MAT",
#                                variable == "gamma_rich" ~ "Gamma Diversity",
#                                variable == "anpp" ~ "ANPP",
#                                variable == "HumanDisturbance" ~ "Human Footprint Index",
#                                variable == "N_Deposition" ~ "N Deposition")) 
# 
# mode.labs <- c("monodominated", "codominated", "even")
# names(mode.labs) <- c("Monodominated", "Codominated", "Even")
# 
# # plot histogram 
# ggplot(df_hist, aes(value)) + # df_hist comes from formatted df above
#   geom_histogram(aes(fill = mode.levels), bins = 50) +
#   facet_wrap(~ variable1,
#              scales = "free") +
#   theme_bw() +
#   theme(legend.position = "top") +
#   scale_fill_manual(name = "",
#                     labels = c("Monodominated", "Codominated", "Even"),
#                     values = c("#007BA7", "#A63922", "#D8B573")) +
#   labs(x = "Value",
#        y = "Count")


# saveRDS(modePlot, file = "data/modePlot.rds")
# saveRDS(modeSite, file = "data/modeSite.rds")
# write.csv(modeSite, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\1_first author\\codominance\\data\\modeSite_20250624.csv', row.names=F)

# Interannual Precipitation addendum --------------------------------------
##Calculation of the Coefficient of Variation for Annual Precipitation at Sites
##Precipitation Data from MSWEP, modified code from Ingrid Slette, Feb 19, 2025
# sites <- FinalAridity #
# 
# # make that file a SpatVector
# site <- sites %>% vect(geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
# 
# # list all of the monthly mswep precip data files
# # change this to location to which you downloaded these files
# r_paths <- list.files("MSWEP_Daily",
#                       full.names = TRUE) %>% 
#   sort()
# 
# # make that a SpatRaster
# r <- rast(r_paths)
# 
# # extract daily precip data for each site
# ppt_daily <- terra::extract(r, site, bind = TRUE)
# 
# df <- as.data.frame(ppt_daily)
# 
# names(df) <- c("OID_", "site_code", "site_proj_comm", "MAP", "MAT", "GDiv", "ANPP", "NDep", "HumanFootprint", "Aridity", "lumpMode", "LumpNames", paste0("precip_", time(r)))
# 
# out <- pivot_longer(df, -c("OID_", "site_code", "site_proj_comm", "MAP", "MAT", "GDiv", "ANPP", "NDep", "HumanFootprint", "Aridity", "lumpMode", "LumpNames"), names_to = "date",
#                     values_to = "precip") %>% 
#   mutate(date = str_replace(date, "^precip_", ""))
# 
# # create a new column for the year
# out$year <- substr(out$date, 1, 4)
# 
# # create a new column for the month
# out$month <- substr(out$date, 6, 7)
# 
# # create a new column for the day
# out$day <- substr(out$date, 9, 10)
# 
# #calculate cv for annual precipitation at each site
# interannual_precip <- filter(out, year != 2020 & year != 2025) %>% group_by(OID_,year) %>% 
#   summarise(annual_total = sum(precip, na.rm = T)) %>% 
#   summarise(
#     mean_annual = mean(annual_total, na.rm = TRUE),
#     sd_annual = sd(annual_total, na.rm = TRUE),
#     cv_Precip = sd_annual / mean_annual)
# 
# Interannual.precip <- merge(FinalAridity, interannual_precip[,c(1,4)], by = "OID_" )
# df_IAP <- as.data.frame(Interannual.precip)
# write.csv(Interannual.precip, file = "C:/Users/msgrabda/Downloads/IAP.csv")
# 

# Fig: Correlation matrix of all the response variables ------------------------

# shows that correlations are statistically significant but weak-moderately correlated with one another
chart.Correlation(df_IAP[, c("MAP", 
                              "MAT", 
                              "GDiv", 
                              "ANPP", 
                              "HumanFootprint", 
                              "NDep", 
                              "Aridity",
                              "cv_Precip")],
                  method = "pearson",
                  histogram = TRUE,
                  cex = 10)

# Multinomial Analysis ----------------------------------------------------------

#multinom function comes from nnet
#mutinomial model using monodominance as the reference state
multinom.baseline1 <- multinom(factor(df_aridity$lumpMode,
                               levels = c(1,2,4)) ~ Aridity + 
                                                    cv_Precip +
                                                    ANPP + 
                                                    GDiv * 
                                                    NDep * 
                                                    HumanFootprint ,
                              data = df_IAP)

#checking p-value manually
(coef <- summary(multinom.baseline1)$coefficients)

(stderr <- summary(multinom.baseline1)$standard.errors)

z <- coef/stderr
(p_values <- 2 * (1 - pnorm(abs(z))))



#PR2() from 'pscl' package provides many R2 estimates
pR2(multinom.baseline1) #McFaddenR2 = 0.337

#Putting into more convenient table format, P-values and OR's match model output 
results <- tidy(multinom.baseline1, conf.int = TRUE, exponentiate = TRUE) #tidy() is a 'broom' function

results_table <- results %>%
  mutate(OR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
    p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))) %>%
  select(y.level, term, OR_CI,p_value) %>%
  rename( "Outcome Category" = y.level,
    "Predictor" = term,
    "Odds Ratio (95% CI) vs. Monodominated" = OR_CI,
    "P-value" = p_value)


# Format: for figure of multinomial model predictions----------------------------------------------------------

# Generate mean (later used in model predicted data)
df_om <- df_aridity %>% 
  na.omit() %>% 
  mutate(Aridity_mean = mean(Aridity),
         GDiv_mean = mean(GDiv),
         HumanFootprint_mean = mean(HumanFootprint),
         NDep_mean = mean(NDep),
         ANPP_mean = mean(ANPP)) 

# Variables of interest
var <- c("Aridity", "GDiv", "HumanFootprint", "NDep", "ANPP") 

# Generate sequence of values looped for each variable 
df_seq <- foreach(v = var, .combine = bind_cols) %do% {
  
  df_v <- df_om %>% 
    dplyr::select(v)
  
  out <- as_tibble(seq(from = min(df_v), 
             to = max(df_v), 
             length.out = 100))
  
}

# Set column names 
colnames(df_seq) <- var

# Add identifier column to link dataframes 
df_seq <- df_seq %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100))

# Combine dataframes- variable means and sequence 
df_r <- df_om %>% 
  dplyr::select(ends_with("_mean")) %>% 
  slice(1:100) %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100)) %>% 
  full_join(df_seq, by = "seq") 

# Mean of variables of interest
var2 <- c("Aridity_mean", "GDiv_mean", "HumanFootprint_mean", "NDep_mean", "ANPP_mean") 

# Predict data using model, sequence, and mean
df_predicted <- foreach(v = var, 
                        v2 = var2,
                        .combine = bind_cols) %do% {
 # Select variable sequence
   df_v <- df_r %>% 
     dplyr::select(v, seq)
 # Select means for all other variables   
   df_m <- df_r %>% 
     dplyr::select(ends_with("_mean"), seq) %>% 
     dplyr::select(!v2)
 # Combine data sets to predict from   
   df_s <- df_v %>% 
     full_join(df_m, by = "seq") %>% 
     dplyr::select(!seq) %>% 
     rename_with(~str_remove(., '_mean'))
 # Predict data
   df_p <- cbind(df_s,
                 data.frame(predict(multinom.baseline1,
                                    newdata = df_s,
                                    type = "probs")))

   
 # Format for figure interpretation   
   df_c <- df_p %>% 
     dplyr::select(v, starts_with("X")) %>% 
     pivot_longer(cols = starts_with("X"), 
                  names_to = "Codom", 
                  values_to = "Probability")
                        }



# figure out what is going on here and if it can be included in fig 2
#ci <- as.data.frame(confint(multinom.baseline1, level = 0.95))



# Clarify column names 
df_combined <- df_predicted %>% 
  select(Codom...2, Aridity, GDiv, HumanFootprint, NDep, ANPP,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         Aridity = Aridity,
         "Gamma Diversity" = GDiv,
         "Human Footprint Index" = HumanFootprint, 
         "N Deposition" = NDep,
         ANPP = ANPP,
         Prob_Aridity = Probability...3,
         Prob_gamma = Probability...6,
         Prob_Human = Probability...9,
         Prob_N = Probability...12,
         Prob_ANPP = Probability...15)

# Assign names
prob <- c("Prob_Aridity", "Prob_gamma", "Prob_ANPP", "Prob_Human", "Prob_N")
named_var <- c("Aridity", "Gamma Diversity", "ANPP", "Human Footprint Index", "N Deposition")



# Fig: Multinomial model predictions --------------------------------------

axis_limits <- list(
  "Aridity" = list(limits = c(0, 40000), breaks = seq(0, 40000, by = 10000)),
  "Gamma Diversity" = list(limits = c(0, 160), breaks = seq(0, 160, by = 50)),
  "ANPP" = list(limits = c(0, 1100), breaks = seq(0, 1100, by = 500)),
  "Human Footprint Index" = list(limits = c(0, 45), breaks = seq(0, 45, by = 10)),
  "N Deposition" = list(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)))

# Pre-allocate list
output <- list()

# Generate figures using loop across variables and their probabilities 
output <- foreach(i = seq_along(named_var), .combine = 'c') %do% {
  v <- named_var[i]
  p <- prob[i]
  
  df_combo <- df_combined %>% 
    dplyr::select(Codom, !!sym(v), !!sym(p)) %>% 
    rename(v1 = !!sym(v), p1 = !!sym(p)) %>% 
    mutate(Codom = case_when(Codom == "X1" ~ "Monodominated",
                             Codom == "X2" ~ "Codominated",
                             Codom == "X4" ~ "Even"),
           Codom = factor(Codom, levels = c("Monodominated", "Codominated", "Even")))
  
  lims <- axis_limits[[v]]$limits
  brks <- axis_limits[[v]]$breaks
  
  fig <- ggplot(df_combo, aes(x = v1, y = p1, color = Codom)) +
    geom_point(size = 4) +
    labs(y = "Probability", x = v) +
    scale_x_continuous(limits = lims, breaks = brks) +
    ylim(0.0, 0.8) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 35),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          legend.position = "right",
          legend.text = element_text(size = 30),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA)) +
    scale_color_manual(name = "",
                       #labels = c("Monodominated", "Codominated", "Even"),
                       values = c("#007BA7", "#A63922", "#D8B573")) +
    guides(color = guide_legend(override.aes = list(size = 7)))
  
  legend <- get_legend(fig)
  fig2 <- fig + theme(legend.position = "none")
  fig_q1 <- ggExtra::ggMarginal(fig2, type = 'density', margins = 'y',
                                size = 5, groupColour = TRUE, groupFill = TRUE)
  
  list(fig_q1)
}
fig

# Overlay the legend on the top-right of the first plot
plot_with_legend <- grobTree(
  output[[1]],  # already a gtable/grob from ggMarginal
  grobTree(legend, vp = viewport(x = 1.16, y = 1.39, just = c("right", "top"))))

# Replace only the first plot with overlaid legend
output[[1]] <- plot_with_legend


out_hist <- foreach(h = named_var) %do% {
  df_o <- df_h %>%
    filter(variable1 == h) %>% 
    dplyr::select(variable1, LumpNames, value) %>% 
    na.omit()
  
  lims <- axis_limits[[h]]$limits
  brks <- axis_limits[[h]]$breaks
  
  fig_h <- ggplot(df_o, aes(value)) +
    geom_histogram(aes(fill = LumpNames), bins = 50, position = "identity", alpha = 0.7) +
    scale_x_continuous(limits = lims, breaks = brks) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 35),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.margin = unit(c(4, 3.4, -0.5, 0.5), "cm"), #adjusts white space around histograms
          plot.background = element_rect(fill = "white", colour = "white")) +
    scale_fill_manual(name = "",
                      values = c("#007BA7", "#A63922", "#D8B573")) +
    labs(x = "", y = "Count")
}


# Arrange all plots
final_plot <- grid.arrange(out_hist[[1]], out_hist[[2]], out_hist[[3]],  
             output[[1]], output[[2]], output[[3]], 
             out_hist[[4]], out_hist[[5]], #out_hist[[6]], 
             output[[4]], output[[5]], #output[[6]], 
             nrow = 4, ncol = 3, 
             heights = c(2.5, 3, 2.5, 3)) # adjusts height of plots

#ggsave("Fig2_model.png", final_plot, width = 28.5, height = 18, dpi = 400)
# save figure as png
ggsave("Fig2_model_2.0.png", final_plot, width = 28.5, height = 18, dpi = 400)



# Fig: Distribution of codominance across sites ----------------------------------------------------------

mapData <- modeSite %>% 
  mutate(codom_category = ifelse(lumpMode == 1, "Mono", 
                          ifelse(lumpMode == 2, "Codom", "Even"))) %>% 
  filter(!is.na(N_Deposition))

mapData$codom_category <- factor(mapData$codom_category, levels = c("Mono", "Codom", "Even"))

autumnalPalette <- c("#007BA7", "#A63922", "#D8B573")

cc_bins <- mapData %>% 
  group_by(codom_category) %>% 
  summarise(count = n()) 
cc_binsdf <- as.data.frame(cc_bins)  

Fig1Bars <- cc_binsdf %>% 
  ggplot(aes(x = factor(codom_category, levels = c("Mono", "Codom", "Even"))))+
  geom_bar(aes(y=count),stat = "identity", fill = autumnalPalette)+
  geom_text(aes(y=count,label = count),vjust = -0.1, size = 20)+
  labs(x="",y= "Count")+
  theme_classic2()+
  theme(element_text(size = 20),
        axis.text.x = element_text(size = 40, angle = 0, vjust = 0.5, hjust =0.5),
        axis.title.y = element_text(size = 40),
        axis.text.y  = element_text(size = 40),
        plot.margin = margin(1,.4,.4,.4,"cm"))





# # Determining differences in codominance across ecoregions ----------------------------------------------------------
# 
# modeSiteEcoregion <- modeSite %>% 
#   mutate(lump_mode_cat=ifelse(lumpMode==1, 'Monodominated',
#                        ifelse(lumpMode==2, 'Codominated', 'Even'))) 
# 
# 
# ecoregionModel <- xtabs(~ Formation + lump_mode_cat, data = modeSiteEcoregion)
# 
# print(chisq <- chisq.test(ecoregionModel))
# # X-squared = 36.811, df = 28, p-value = 0.1231
# 
# mosaicplot(ecoregionModel, shade = TRUE, las=2,
#            main = "modeSiteCodom")
# 
# summaryModeSiteCodom <- as.data.frame(ecoregionModel) %>% 
#   group_by(lump_mode_cat ) %>%
#   mutate(percent=Freq/sum(Freq)) %>% 
#   ungroup()
# 
# ggplot(summaryModeSiteCodom, aes(x=lump_mode_cat , y=Formation)) +
#   geom_tile(aes(fill=percent)) +
#   geom_text(aes(label=Freq), size=6, color='grey40') +
#   # geom_text(aes(label=round(100*percent, digits=0))) +
#   scale_fill_gradient(low='#F8FBF8', high='#031B88') +
#   scale_y_discrete(limits=rev) +
#   xlab('Control') + ylab('Treatment') + labs(fill='Column\nPercentage')
