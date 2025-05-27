
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
   separate(exp_unit, into=c('site_code', 'project_name', 'community_type', 'plot_id', 'treatment', 'calendar_year'), sep='::', remove=F) %>% 
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
   summarise(plot_codom=round(mean(plot_codom), digits=0), .groups='drop') %>% # calculate mean for ties
   rbind(singletonCodomPlotYear)

# for plots with singleton ties for modes, calculate mean and round to nearest integer
multipleMode <- numCodomPlotYear %>% 
   select(database, site_code, project_name, community_type, plot_id, trt_type, treatment, num_group, calendar_year) %>% 
   unique() %>% 
   full_join(modePlotTrue) %>% 
   filter(is.na(plot_codom)) %>%
   group_by(database, site_code, project_name, community_type, plot_id, trt_type, treatment) %>% 
   summarise(plot_codom=round(mean(num_group), digits=0), .groups='drop')
 
# bind dataframes for averaged ties and true modes at plot level
modePlot <- rbind(modePlotTrue, multipleMode)


# Calculate mode across all control plots for each experiment ------------------

# create a dataframe of codominant group for experiments with a single plot (because can't calculate mode of singleton)
singletonCodomPlot <- modePlot %>% 
  filter(trt_type=='control') %>% 
  group_by(database, site_code, project_name, community_type) %>% 
  mutate(length=length(community_type)) %>% 
  ungroup() %>% 
  filter(length==1) %>% 
  rename(mode_site=plot_codom) %>% 
  dplyr::select(database, site_code, project_name, community_type, mode_site) 

# calculate mode across plots for each experiment, dropping those with ties
modeSiteTrue <- modePlot %>%
   filter(trt_type=='control') %>%
   group_by(database, site_code, project_name, community_type) %>% # mode generated from these
   reframe(mode_site = DescTools::Mode(plot_codom)) %>%  
   ungroup() %>% 
   group_by(database, site_code, project_name, community_type) %>% 
   summarise(mode_site=round(mean(mode_site), digits=0), .groups='drop') %>% # calculate mean for ties
   filter(!is.na(mode_site)) %>% 
   rbind(singletonCodomPlot)

# for plots with ties for modes, calculate mean and round to nearest integer
multipleModeProj <- modePlot %>% 
  filter(trt_type=='control') %>% 
  select(database, site_code, project_name, community_type, plot_id, plot_codom) %>% 
  full_join(modeSiteTrue) %>% 
  filter(is.na(mode_site)) %>%
  group_by(database, site_code, project_name, community_type) %>% 
  summarise(mode_site=round(mean(plot_codom), digits=0), .groups='drop')

# bind dataframes for average ties and true modes at site level
modeSite <- rbind(modeSiteTrue, multipleModeProj) %>% 
  left_join(readRDS("data/envData.rds")) %>% 
  mutate(lumpMode = ifelse(mode_site == 3, 2, mode_site)) 


# Histogram- count of codom level per variable ----------------------------

df_hist <- modeSite %>% 
  select(mode_site, lumpMode, MAP, MAT, gamma_rich, anpp, HumanDisturbance, N_Deposition) %>% 
  pivot_longer(cols = c("MAP", "MAT", "gamma_rich", "anpp", "HumanDisturbance", "N_Deposition"), 
               names_to = "variable", values_to = "value") %>% 
  mutate(lumpMode = ifelse(lumpMode == 1, "Monodominated", 
                           ifelse(lumpMode == 2, "Codominated", "Even")),
         mode.levels = factor(lumpMode, levels = c("Monodominated", "Codominated", "Even")),
         variable1 = case_when(variable == "MAP" ~ "MAP",
                               variable == "MAT" ~ "MAT",
                               variable == "gamma_rich" ~ "Gamma Diversity",
                               variable == "anpp" ~ "ANPP",
                               variable == "HumanDisturbance" ~ "Human Disturbance",
                               variable == "N_Deposition" ~ "N Deposition")) 

mode.labs <- c("monodominated", "codominated", "even")
names(mode.labs) <- c("Monodominated", "Codominated", "Even")

# bar counts are all the same per variable ...?
# this is at the site level 
ggplot(df_hist, aes(mode.levels)) +
  geom_bar(aes(fill = mode.levels)) +
  facet_wrap(~ variable1,
             labeller = labeller(mode.levels = mode.labs)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(name = "",
                    labels = c("Monodominated", "Codominated", "Even"),
                    values = c("#007BA7", "#A63922", "#D8B573")) +
  labs(x = "Codominance Level",
       y = "Count")

# histogram of variable values
ggplot(df_hist, aes(value)) + # df_hist comes from formatted df above 
  geom_histogram(aes(fill = mode.levels), bins = 50) +
  facet_wrap(~ variable1,
             scales = "free") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(name = "",
                    labels = c("Monodominated", "Codominated", "Even"),
                    values = c("#007BA7", "#A63922", "#D8B573")) +
  labs(x = "Value",
       y = "Count")

# saveRDS(modeSite, file = "data/modeSite.rds")


# Visualizing each predictor with boxplots ----------------------------------------------------------
ggplot(modeSite, aes(x=as.factor(mode_site), y=abs(Latitude)))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=MAP))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=MAT))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=gamma_rich))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=anpp))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=HumanDisturbance))+
  geom_boxplot()+
  coord_flip()

ggplot(modeSite, aes(x=as.factor(mode_site), y=N_Deposition))+
  geom_boxplot()+
  coord_flip()

factors <-modeSite[,c(6:13)]
chart.Correlation(factors, method = "spearman")


# Multinomial Analysis ----------------------------------------------------------

multinom.baseline1 <- multinom(factor(modeSite$lumpMode,
                     levels = c(1,2,4)) ~ MAP+gamma_rich+HumanDisturbance+N_Deposition+MAP*MAT*anpp+gamma_rich*N_Deposition*HumanDisturbance ,
              data=modeSite)
a <- multinom(factor(modeSite$lumpMode,
                        levels = c(1,2,4)) ~ MAP+gamma_rich+HumanDisturbance+N_Deposition+MAP*MAT*anpp+gamma_rich*N_Deposition*HumanDisturbance ,
                 data=modeSite)
stepAIC(a, direction = "backward")

coef <- summary(a)$coefficients
coef

stderr <- summary(a)$standard.errors
stderr

z <- coef/stderr
p_values <- 2 * (1 - pnorm(abs(z)))
p_values

exp(coef(a))

head(round(fitted(a),2))


# Code to generate figure of multinomial model----------------------------------------------------------

# Generate mean (later used in model predicted data)
df_om <- modeSite %>% 
  na.omit() %>% 
  mutate(MAP_mean = mean(MAP),
         MAT_mean = mean(MAT),
         gamma_rich_mean = mean(gamma_rich),
         HumanDisturbance_mean = mean(HumanDisturbance),
         N_Deposition_mean = mean(N_Deposition),
         anpp_mean = mean(anpp)) 

# Variables of interest
var <- c("MAP", "MAT", "gamma_rich", "HumanDisturbance", "N_Deposition", "anpp") 

# Generate sequence of values looped for each variable 
df_seq <- foreach(v = var, .combine = bind_cols) %do% {
  
  df_v <- df_om %>% 
    select(v)
  
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
  select(ends_with("_mean")) %>% 
  slice(1:100) %>% 
  add_column(seq = seq(from = 1, to = 100, length.out = 100)) %>% 
  full_join(df_seq, by = "seq") 

# Mean of variables of interest
var2 <- c("MAP_mean", "MAT_mean", "gamma_rich_mean", "HumanDisturbance_mean", "N_Deposition_mean", "anpp_mean") 

# Predict data using model, sequence, and mean
df_predicted <- foreach(v = var, 
                        v2 = var2,
                        .combine = bind_cols) %do% {
 # Select variable sequence
   df_v <- df_r %>% 
    select(v, seq)
 # Select means for all other variables   
   df_m <- df_r %>% 
     select(ends_with("_mean"), seq) %>% 
     select(!v2)
 # Combine data sets to predict from   
   df_s <- df_v %>% 
     full_join(df_m, by = "seq") %>% 
     select(!seq) %>% 
     rename_with(~str_remove(., '_mean'))
 # Predict data
   df_p <- cbind(df_s,
                 data.frame(predict(multinom.baseline1,
                                    newdata = df_s,
                                    type = "probs")))
 # Format for figure interpretation   
   df_c <- df_p %>% 
     select(v, starts_with("X")) %>% 
     pivot_longer(cols = starts_with("X"), 
                  names_to = "Codom", 
                  values_to = "Probability")
}

# Clarify column names 
df_combined <- df_predicted %>% 
  select(Codom...2, MAP, MAT, gamma_rich, HumanDisturbance, N_Deposition, anpp,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         "Gamma Diversity" = gamma_rich,
         "Human Disturbance" = HumanDisturbance, 
         "N Deposition" = N_Deposition,
         ANPP = anpp,
         Prob_MAP = Probability...3,
         Prob_MAT = Probability...6,
         Prob_gamma = Probability...9,
         Prob_Human = Probability...12,
         Prob_N = Probability...15,
         Prob_anpp = Probability...18)

# Assign names
prob <- c("Prob_MAP", "Prob_MAT", "Prob_gamma", "Prob_anpp", "Prob_Human", "Prob_N")
named_var <- c("MAP", "MAT", "Gamma Diversity", "ANPP", "Human Disturbance", "N Deposition")

axis_limits <- list(
  "MAP" = list(limits = c(0, 2500), breaks = seq(0, 2500, by = 500)),
  "MAT" = list(limits = c(-10, 30), breaks = seq(-10, 30, by = 10)),
  "Gamma Diversity" = list(limits = c(0, 250), breaks = seq(0, 250, by = 50)),
  "ANPP" = list(limits = c(0, 1500), breaks = seq(0, 1500, by = 500)),
  "Human Disturbance" = list(limits = c(0, 50), breaks = seq(0, 50, by = 10)),
  "N Deposition" = list(limits = c(0, 2000), breaks = seq(0, 2000, by = 500))
)

# Pre-allocate list
output <- list()

# Generate figures using loop across variables and their probabilities 
output <- foreach(i = seq_along(named_var), .combine = 'c') %do% {
  v <- named_var[i]
  p <- prob[i]
  
  df_combo <- df_combined %>% 
    select(Codom, !!sym(v), !!sym(p)) %>% 
    rename(v1 = !!sym(v), p1 = !!sym(p))
  
  lims <- axis_limits[[v]]$limits
  brks <- axis_limits[[v]]$breaks
  
  fig <- ggplot(df_combo, aes(x = v1, y = p1, color = Codom)) +
    geom_point(size = 2) +
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
                       labels = c("Monodominated", "Codominated", "Even"),
                       values = c("#007BA7", "#A63922", "#D8B573"))+
    guides(color = guide_legend(override.aes = list(size = 7)))
  
  legend <- get_legend(fig)
  fig2 <- fig + theme(legend.position = "none")
  fig_q1 <- ggExtra::ggMarginal(fig2, type = 'density', margins = 'y',
                                size = 5, groupColour = TRUE, groupFill = TRUE)
  
  list(fig_q1)
}


# Overlay the legend on the top-right of the first plot
plot_with_legend <- grobTree(
  output[[1]],  # already a gtable/grob from ggMarginal
  grobTree(legend, vp = viewport(x = 1.16, y = 1.39, just = c("right", "top"))))

# Replace only the first plot with overlaid legend
output[[1]] <- plot_with_legend


out_hist <- foreach(h = named_var) %do% {
  df_o <- df_hist %>%
    filter(variable1 == h) %>% 
    select(variable1, mode.levels, value) %>% 
    na.omit()
  
  lims <- axis_limits[[h]]$limits
  brks <- axis_limits[[h]]$breaks
  
  fig_h <- ggplot(df_o, aes(value)) +
    geom_histogram(aes(fill = mode.levels), bins = 50, position = "identity", alpha = 0.7) +
    scale_x_continuous(limits = lims, breaks = brks) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 35),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.margin = unit(c(2, 3.4, 0, 0.5), "cm"),
          plot.background = element_rect(fill = "white", colour = "white")) +
    scale_fill_manual(name = "",
                      labels = c("Monodominated", "Codominated", "Even"),
                      values = c("#007BA7", "#A63922", "#D8B573")) +
    labs(x = "", y = "")
}

# Show all plots
grid.arrange(grobs = output, nrow = 2, ncol = 3)


# Generate for loop for variable value distribution
out_hist <- foreach(h = named_var) %do% {
  
  # Prepare dataframe
  df_o <- df_hist %>%
    select(variable1, mode.levels, value) %>% 
    filter(variable1 == h) %>% 
    na.omit()
  
  lims <- axis_limits[[h]]$limits
  brks <- axis_limits[[h]]$breaks
  
  fig_h <- ggplot(df_o, aes(value)) +
    geom_histogram(aes(fill = mode.levels), bins = 50, position = "identity", alpha = 0.7) +
    scale_x_continuous(limits = lims, breaks = brks) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 35),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.margin = unit(c(2, 3.4, 0, 0.5), "cm"),
          plot.background = element_rect(fill = "white", colour = "white")) +
    scale_fill_manual(name = "",
                      labels = c("Monodominated", "Codominated", "Even"),
                      values = c("#007BA7", "#A63922", "#D8B573")) +
    labs(x = "", y = "")
}

# Arrange all plots
final_plot <- grid.arrange(out_hist[[1]], out_hist[[2]], out_hist[[3]],  
                           output[[1]], output[[2]], output[[3]], 
                           out_hist[[4]], out_hist[[5]], out_hist[[6]], 
                           output[[4]], output[[5]], output[[6]], 
                           nrow = 4, ncol = 3, heights = c(2, 3, 2, 3))




# Visualizing distribution of codominance across sites ----------------------------------------------------------

mapData <- modeSite %>% 
  mutate(codom_category = ifelse(mode_site == 1, "Monodominated", 
                          ifelse(mode_site == 2, "Codominated", "Even"))) %>% 
  filter(!is.na(N_Deposition))

mapData$codom_category <- factor(mapData$codom_category, levels = c("Monodominated", "Codominated", "Tridominated", "Even"))

autumnalPalette <- c("#007BA7", "#A63922", "#D8B573", 'grey')

mapData %>% 
  ggplot(aes(x="", y=mode_site, fill = codom_category)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y",start=0)

ggplot(mapData, aes(x = codom_category, fill = codom_category)) +
  geom_histogram(stat = "count") +
  stat_count(binwidth = 1, 
             geom = 'text', 
             color = 'white', 
             aes(label = after_stat(count)),
             position = position_stack(vjust = 0.5))+
  xlab("Site Dominance")+
  ylab("Count")+
  scale_fill_manual(values = autumnalPalette) +
  theme(legend.position = "none")

