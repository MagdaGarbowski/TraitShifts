# project: Trait-Shifts 
# objective: merge datasets prior to analyses
# author: Magda Garbowski 
# date: November 22, 2023

library(data.table)

# ----------------------------------- load data -------------------------------------

SPCIS_CWM <- as.data.frame(fread("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM.csv", 
                           select = c("Plot", "Year", "Trait", "CWM")))

SPCIS_traits <- as.data.frame(fread("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_wtraits.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "Trait", "Value")))

SPCIS_plot_info <- as.data.frame(fread("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_trait_coverage.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

SPCIS_ecoregion_IV <- as.data.frame(fread("/Users/magdagarbowski/TraitShifts/Data/SPCIS_ecoreg_03182023.csv", 
                                          select = c("Plot", "US_L4CODE", "US_L4NAME", "Lat", "Long")))

# -------------------------------------- functions ---------------------------------------

# return plots with both N and I 
inv_nat_present <- function(df){
  inv_levels_length <- length(levels(as.factor(df$NativeStatus)))
  if(inv_levels_length == 2){
    return(df)
  }
}

# weighted averages of N and I  
CWM_function_N_I <- function(df){
  df_splits <- split(df, df$Trait)
  splits_function <- function(df_out){
    df_n <- df_out[df_out$NativeStatus == "N",]
    df_i <- df_out[df_out$NativeStatus == "I",]
    
    cwm_function_2 <- function(df){
      df$Value <- as.numeric(df$Value)
      tot_cov = sum(df$PctCov_100)
      df$rel_100_ss = df$PctCov_100/tot_cov * 100
      df$weighted_mean = (df$Value * df$rel_100_ss)/100
      cwm = sum(df$weighted_mean)
    }
    
    cwm_i_n <- as.data.frame(do.call(cbind, lapply(list(df_i, df_n), cwm_function_2)))
    colnames(cwm_i_n) <- c("CWM_I", "CWM_N")
    
    cwm_i_n$Plot <- df_out$Plot[1]
    cwm_i_n$Year <- df_out$Year[1]
    cwm_i_n$Trait <- df_out$Trait[1]
    cwm_i_n$inv_relcov <- sum(df_out[df_out$NativeStatus == "I",]$rel_100)
    cwm_i_n$nat_relcov <- sum(df_out[df_out$NativeStatus == "N",]$rel_100)
    return(cwm_i_n)
  }
  out <- do.call(rbind, lapply(df_splits, splits_function))
  return(out)
}

# ----------------------------------- merge datasets -------------------------------------

# merge the datasets and subset plots with >80% coverage by trait 
CWM_plot_info <- merge(SPCIS_CWM, SPCIS_plot_info, by = c("Plot", "Year", "Trait"))
CWM_plot_info_80 <- CWM_plot_info[CWM_plot_info$trait_cov > 80,]
CWM_plot_info_80 <-  CWM_plot_info_80[CWM_plot_info_80$EcoRegionLevelI %in% 
                                        c("EASTERN TEMPERATE FORESTS", "GREAT PLAINS", "NORTHERN FORESTS",
                                          "MEDITERRANEAN CALIFORNIA", "NORTH AMERICAN DESERTS", "NORTHWESTERN FORESTED MOUNTAINS"),]

# drop alaska plots 
CWM_plot_info_80 <- CWM_plot_info_80[!CWM_plot_info_80$Lat > 60,]

# ----------------------------------- abundance gradient datasets -------------------------------------
# (1) full
# (2) >2% introduced species abundance 

# full
CWM_plot_info_80 <- merge(CWM_plot_info_80, SPCIS_ecoregion_IV, by = "Plot", all.x = TRUE)
CWM_plot_info_80 <- CWM_plot_info_80[!CWM_plot_info_80$Lat > 60,]

# 2% > invasion 
CWM_plot_info_80_2percent <- CWM_plot_info_80[CWM_plot_info_80$I_relcov > 1.999,]

# drop alaska plots 
CWM_plot_info_80_2percent <- CWM_plot_info_80_2percent[!CWM_plot_info_80_2percent$Lat > 60,]

# ------------------------------- co-occurring natives and introduced dataset --------------------------

# merge datasets traits dataset with 80% coverage dataset
dat_80_cov <- merge(SPCIS_traits, CWM_plot_info_80, by.x = c("Plot", "Trait", "Year"), 
                    by.y = c("Plot", "Trait", "Year"))

dat_80_splits <- split(dat_80_cov, list(dat_80_cov$Plot), drop = TRUE) 

N_I_plots <- lapply(dat_80_splits, inv_nat_present)
N_I_plots[sapply(N_I_plots, is.null)] <- NULL

inv_nat_cwm <- do.call(rbind, lapply(N_I_plots, CWM_function_N_I))

# drop plots that do not have weighted avgs for both natives and invasive 
inv_nat_cwm <- inv_nat_cwm[!(inv_nat_cwm$inv_relcov == 0 | inv_nat_cwm$nat_relcov == 0),]

# merge with ecoregion level 1 and 4  
eco_r <- unique(CWM_plot_info_80[c("Plot", "Year", "EcoRegionLevelI", "US_L4NAME", "Lat", "Long")])
inv_nat_cwm_er <- merge(inv_nat_cwm, eco_r, by = c("Plot", "Year"))

inv_nat_cwm_er <- inv_nat_cwm_er[!inv_nat_cwm_er$Lat > 60,]

write.csv(CWM_plot_info_80, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80.csv")
write.csv(CWM_plot_info_80_2percent, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80_2percent.csv")
write.csv(inv_nat_cwm_er, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_cooccurring_weighted_means.csv")

write.csv(CWM_plot_info_80, "/Users/magdagarbowski/TraitShifts/Figshare_datasets/SPCIS_CWM_80.csv")
write.csv(CWM_plot_info_80_2percent, "/Users/magdagarbowski/TraitShifts/Figshare_datasets/SPCIS_CWM_80_2percent.csv")
write.csv(inv_nat_cwm_er, "/Users/magdagarbowski/TraitShifts/Figshare_datasets/SPCIS_cooccurring_weighted_means.csv")
