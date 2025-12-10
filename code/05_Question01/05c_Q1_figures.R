# figures for multinomial model with aridity and precip

# Fig: Distribution of values ---------------------------------------------

## must mutate and factor levels so each group is treated as such 
#change this data file to represent the cutoff of interest
df_iap <- readRDS("data/modeSite.rds") %>% 
  rename(ANPP=anpp,
         GDiv=gamma_rich,
         NDep=NDeposition) %>% 
  mutate(LumpNames=factor(case_when(lumpMode==1 ~ 'Monodominated',
                                        lumpMode==2 ~ 'Codominated', 
                                        lumpMode==4 ~ 'Even'),
                             levels=c('Monodominated', 'Codominated', 'Even')),
         lumpMode=factor(lumpMode, levels=c(1,2,4)))

df_arid <- df_iap %>% 
  dplyr::select(lumpMode, LumpNames, MAP, MAT, GDiv, ANPP, HumanFootprint, NDep, Aridity, cv_Precip)

df_h <- df_arid %>% 
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

# plot histogram 
ggplot(df_h, aes(value)) + # df_hist comes from formatted df above 
  geom_histogram(aes(fill = LumpNames), bins = 50) +
  facet_wrap(~ variable1,
             scales = "free") +
  theme_bw() +
  theme(legend.position = "top")+
  scale_fill_manual(name = c("Monodominated","Codominated","Even"),
                    values = c( "#007BA7","#A63922", "#D8B573")) +
  labs(x = "Value",
       y = "Count")


# Fig: Correlation matrix of all the response variables ------------------------

# shows that correlations are statistically significant but weak-moderately correlated with one another
(fig_corr <- chart.Correlation(df_iap[, c("MAP", 
                                          "MAT", 
                                          "GDiv", 
                                          "ANPP", 
                                          "HumanFootprint", 
                                          "NDep", 
                                          "Aridity",
                                          "cv_Precip")],
                               method = "pearson",
                               histogram = TRUE))


# Fig: Multinomial model predictions --------------------------------------

axis_limits <- list(
  "MAP" = list(limits = c(0, 2600), breaks = seq(0, 2600, by = 1000)),
  "MAT" = list(limits = c(-12, 30), breaks = seq(-12, 30, by = 10)),
  "Gamma Diversity" = list(limits = c(0, 250), breaks = seq(0, 250, by = 50)),
  "ANPP" = list(limits = c(0, 1100), breaks = seq(0, 1100, by = 500)),
  "Human Footprint Index" = list(limits = c(0, 45), breaks = seq(0, 45, by = 10)),
  "N Deposition" = list(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)),
  "Interannual Precip" = list(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)))

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
    geom_smooth(size = 2) +
    geom_point(size = 0.0003)+
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
          legend.text = element_text(size = 24),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA)) +
    scale_color_manual(name = "",
                       labels = c("Monodominated", "Codominated", "Even"),
                       values = c("#007BA7", "#A63922", "#D8B573")) +
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

final_plot <- grid.arrange(out_hist[[1]], out_hist[[2]], out_hist[[3]], nullGrob(),  
                           output[[1]], output[[2]], output[[3]], nullGrob(),
                           out_hist[[4]], out_hist[[5]], out_hist[[6]], out_hist[[7]], 
                           output[[4]], output[[5]], output[[6]], output[[7]],
                           nrow = 4, ncol = 4, 
                           heights = c(2.5, 3, 2.5, 3)) # adjusts height of plots

#ggsave("Fig2_model.png", final_plot, width = 28.5, height = 18, dpi = 400)
# save figure as png
ggsave("Fig2.png", final_plot, width = 28.5, height = 18, dpi = 400)



# Fig: Distribution of codominance across sites ----------------------------------------------------------

mapData <- df_iap %>% 
  mutate(codom_category = ifelse(lumpMode == 1, "Mono", 
                                 ifelse(lumpMode == 2, "Codom", "Even"))) %>% 
  filter(!is.na(NDep))

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




