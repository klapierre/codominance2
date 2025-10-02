# figure code to plot MAP & MAT from multinomial model 

# Format: for figure of multinomial model predictions----------------------------------------------------------

# Generate mean (later used in model predicted data)
df_om <- df_iap %>% 
  na.omit() %>% 
  mutate(#Aridity_mean = mean(Aridity),
    GDiv_mean = mean(GDiv),
    HumanFootprint_mean = mean(HumanFootprint),
    NDep_mean = mean(NDep),
    ANPP_mean = mean(ANPP),
    #cv_Precip_mean = mean(cv_Precip),
    MAP_mean = mean(MAP),
    MAT_mean = mean(MAT)) 

# Variables of interest
var <- c("GDiv", "HumanFootprint", "NDep", "ANPP", "MAP", "MAT") 

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
var2 <- c( "GDiv_mean", "HumanFootprint_mean", "NDep_mean", "ANPP_mean", "MAP_mean", "MAT_mean") 

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
  select(Codom...2, GDiv, HumanFootprint, NDep, ANPP, MAP, MAT,
         starts_with("Probability")) %>% 
  rename(Codom = Codom...2,
         "Gamma Diversity" = GDiv,
         "Human Footprint Index" = HumanFootprint, 
         "N Deposition" = NDep,
         ANPP = ANPP,
         MAP = MAP,
         MAT = MAT,
         Prob_gamma = Probability...3,
         Prob_Human = Probability...6,
         Prob_N = Probability...9,
         Prob_ANPP = Probability...12,
         Prob_MAP = Probability...15,
         Prob_MAT = Probability...18)

# Assign names
prob <- c( "Prob_gamma", "Prob_ANPP", "Prob_Human", "Prob_N", "Prob_MAP", "Prob_MAT")
named_var <- c( "Gamma Diversity", "ANPP", "Human Footprint Index", "N Deposition", "MAP", "MAT")




# Fig: Multinomial model predictions --------------------------------------
df_iap <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Grad School/Terui Lab/Codominance/IAP.csv") %>% 
  mutate(LumpNames = factor(LumpNames, c("Monodominated", "Codominated", "Even")))

df_sel <- df_iap %>% 
  dplyr::select(lumpMode, LumpNames, lumpMode, MAP, MAT, GDiv, ANPP, HumanFootprint, NDep, Aridity, cv_Precip)

df_h <- df_sel %>% 
  pivot_longer(cols = c("MAP", "MAT", "GDiv", "ANPP", "HumanFootprint", "NDep", "Aridity", "cv_Precip"), 
               names_to = "variable", values_to = "value") %>% 
  mutate(variable1 = case_when(variable == "MAP" ~ "MAP",
                               variable == "MAT" ~ "MAT",
                               variable == "GDiv" ~ "Gamma Diversity",
                               variable == "ANPP" ~ "ANPP",
                               variable == "HumanFootprint" ~ "Human Footprint Index",
                               variable == "NDep" ~ "N Deposition",
                               variable == "Aridity" ~ "Aridity",
                               variable == "cv_Precip" ~ "Precip")) 


axis_limits <- list(
  "Aridity" = list(limits = c(0, 40000), breaks = seq(0, 40000, by = 10000)),
  "Gamma Diversity" = list(limits = c(0, 160), breaks = seq(0, 160, by = 50)),
  "ANPP" = list(limits = c(0, 1100), breaks = seq(0, 1100, by = 500)),
  "Human Footprint Index" = list(limits = c(0, 45), breaks = seq(0, 45, by = 10)),
  "N Deposition" = list(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)),
  "Precip" = list(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)),
  "MAP" = list(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)),
  "MAT" = list(limits = c(-15, 25), breaks = seq(-15, 25, by = 10)))

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
    ylim(0.0, 1) +
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
                           out_hist[[4]], out_hist[[5]], out_hist[[6]], 
                           output[[4]], output[[5]], output[[6]], 
                           nrow = 4, ncol = 3, 
                           heights = c(2.5, 3, 2.5, 3)) # adjusts height of plots

#ggsave("Fig2_model.png", final_plot, width = 28.5, height = 18, dpi = 400)
# save figure as png
ggsave("Fig2_model_2.0.png", final_plot, width = 28.5, height = 18, dpi = 400)
