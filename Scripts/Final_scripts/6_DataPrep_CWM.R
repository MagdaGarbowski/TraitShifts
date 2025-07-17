# project: Trait-Shifts 
# objective: get community weighted means by plot
# author: Magda Garbowski 
# date: November 22, 2023

library(data.table)

# ------------------------------- data ---------------------------------------------

traits <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Submitted_datasets/SPCIS_wtraits.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100",
                                      "rel_100", "Trait", "Units", "Value")))

coverage <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Submitted_datasets/SPCIS_trait_coverage.csv"))

# ------------------------------- functions -----------------------------------------

# relative cover for only species with trait data 
rel_cov_function <- function(df){
  df_splits <- split(df, df$Trait)
  rel_cov_inner <- function(df_inner){
    tot_cover_std <- sum(df_inner$rel_100)
    df_inner$rel_100_std <- (df_inner$rel_100/tot_cover_std) * 100
    return(df_inner)
  }
  out <- do.call(rbind, lapply(df_splits, rel_cov_inner))
  return(out)
}

# CWM function using relative cover of only species with trait data
CWM_function <- function(df){
  yrs = length(levels(as.factor(df$Year)))
  CWM_inner_function <- function(trait_df){
    trait_df$weighted_mean <- (trait_df$rel_100_std * trait_df$Value)/100
    out <- data.frame(Plot = trait_df$Plot[1], 
                      Year = trait_df$Year[1], 
                      Trait = trait_df$Trait[1],
                      CWM = sum(trait_df$weighted_mean, na.rm = TRUE))
    return(out)
  }
  
  if(yrs == 1){
    df_trait_splits <- split(df, list(df$Trait), drop = TRUE)
    out_all <- do.call(rbind, lapply(df_trait_splits, CWM_inner_function))
    return(out_all)
  }
  
  if(yrs > 1){
    df_yr_splits <- split(df, list(df$Year, df$Trait), drop = TRUE)
    out_all <- do.call(rbind, lapply(df_yr_splits, CWM_inner_function))
    return(out_all)
  }
}

# ----------------------------- community weighted means  --------------------------------------

# select out plots with >80 coverage for traits excluding "NI"  
coverage_80 <- coverage[coverage$trait_cov > 80,]
coverage_80_ss <- coverage_80[c("Plot", "Trait", "Year", "trait_cov")]

# merge datasets to select plots that have >80 coverage for specific traits 
dat_80_cov <- merge(traits, coverage_80_ss, by.x = c("Plot", "Trait", "Year"), 
                    by.y = c("Plot", "Trait", "Year"), all.y = TRUE)

# recalculate relative cover with just species that have trait values
dat_80_splits <- split(dat_80_cov, list(dat_80_cov$Plot, dat_80_cov$Year), drop = TRUE)
spcis_data_relcov_ls <- lapply(dat_80_splits, rel_cov_function)
spcis_data_relcov <- do.call(rbind, spcis_data_relcov_ls)

# relativized CWM 
dat_splits <- split(spcis_data_relcov, list(spcis_data_relcov$Plot), drop = TRUE)
CWM_out <- lapply(dat_splits, CWM_function)
CWM_out_df <- do.call(rbind, CWM_out)

# -------------------------------- write csv ---------------------------------------------------

write.csv(CWM_out_df, "/Users/MagdaGarbowski 1/TraitShifts/Submitted_datasets/SPCIS_CWM.csv")
