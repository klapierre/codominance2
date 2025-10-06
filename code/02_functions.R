# source code for functions


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

# #mode
# mode <- function(codes){
#   which.max(tabulate(codes))
# }
# 
# # mode (used to calculate mode for codominance groups across years)
# Mode <- function(x, na.rm = FALSE) {
#   if(na.rm){
#     x = x[!is.na(x)]
#   }
# 
#   ux <- unique(x)
#   return(ux[which.max(tabulate(match(x, ux)))])
# }

get_dist <- function(data, pool, md) {
  
  n_pair <- n_distinct(data$pair_id)
  df_dist <- foreach(i = 1:n_pair,
                     .combine = bind_rows) %do% {
                       
                       ## select a given set of co-dominants
                       df_pair <- data %>% 
                         filter(pair_id == i)
                       
                       ## pairwise co-dominants
                       codom <- pull(df_pair, species)
                       m_combo <- combn(codom, 2)
                       
                       ## calculate trait distance for each pair
                       d <- NULL
                       for (j in seq_len(ncol(m_combo))) {
                         v_codom <- m_combo[, j]
                         key <- sapply(v_codom,
                                       function(x) which(x == pool))
                         d[j] <- md[key[1], key[2]]
                       }
                       
                       ## distance distribution for the pool
                       vd <- md[lower.tri(md)]
                       mu_d <- mean(vd)
                       sd_d <- sd(vd)
                       ses <- (d - mu_d) / sd_d
                       
                       ## return pairwise distance and proportion > observed
                       distinct(df_pair, 
                                site_proj_comm) %>% 
                         bind_cols(tibble(sp1 = m_combo[1, ],
                                          sp2 = m_combo[2, ],
                                          dist = d,
                                          ses = ses,
                                          p = sapply(d, function(x) mean(x > vd)),
                                          n_pool = length(vd),
                                          n_obs = sapply(d, function(x) sum(x > vd)))
                         ) %>% 
                         return()
                     }
  
  return(df_dist)
}