# Trait-Shifts 
# goal: get CWM for each plot by trait combo 

library(data.table)
library(ggplot2)
library(gridExtra)

dat <- fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wTRY.csv", 
             select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean"))

dat_splits <- split(dat, list(dat$Plot), drop = TRUE)

CWM_function <- function(df){

  yrs = length(levels(as.factor(df$Year)))
  
  CWM_inner_function <- function(trait_df){
    trait_df$weighted_mean <- (trait_df$rel_100 * trait_df$mean)/100
    out <- data.frame(Plot = trait_df$Plot[1], 
                      Year = trait_df$Year[1], 
                      TraitNameAbr = trait_df$TraitNameAbr[1],
                      CWM = sum(trait_df$weighted_mean, na.rm = TRUE))
    return(out)
  }
  
  if(yrs == 1){
    df_trait_splits <- split(df, list(df$TraitNameAbr))
    out_all <- do.call(rbind, lapply(df_trait_splits, CWM_inner_function))
    return(out_all)
  }
  
  if(yrs > 1){
    df_yr_splits <- split(df, list(df$Year, df$TraitNameAbr))
    out_all <- do.call(rbind, lapply(df_yr_splits, CWM_inner_function))
    return(out_all)
  }
}


# this takes 10 minutes 
plot_CWM_out <- lapply(dat_splits, CWM_function)
plot_CWM_out_df <- do.call(rbind, plot_CWM_out)

write.csv(plot_CWM_out_df, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM.csv")
