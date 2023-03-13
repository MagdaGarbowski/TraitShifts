# Trait-Shifts 
# goal: get CWM for each plot by trait combo 

library(data.table)

# ------------------------------- data ---------------------------------------------

dat <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wtraits.csv", 
             select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean")))

coverage <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv"))

# ------------------------------- functions ---------------------------------------------

# relative cover function for only species with trait data 
rel_cov_function <- function(df){
  df_splits <- split(df, df$TraitNameAbr)
  rel_cov_inner <- function(df_inner){
    tot_cover_std <- sum(df_inner$rel_100)
    df_inner$rel_100_std <- (df_inner$rel_100/tot_cover_std) * 100
    return(df_inner)
  }
  out <- do.call(rbind, lapply(df_splits, rel_cov_inner))
}

# CWM function using relative cover (only species with trait values)
CWM_function <- function(df){
  yrs = length(levels(as.factor(df$Year)))
  
  CWM_inner_function <- function(trait_df){
    trait_df$weighted_mean <- (trait_df$rel_100_std * trait_df$mean)/100
    out <- data.frame(Plot = trait_df$Plot[1], 
                      Year = trait_df$Year[1], 
                      TraitNameAbr = trait_df$TraitNameAbr[1],
                      CWM = sum(trait_df$weighted_mean, na.rm = TRUE))
    return(out)
  }
  
  if(yrs == 1){
    df_trait_splits <- split(df, list(df$TraitNameAbr), drop = TRUE)
    out_all <- do.call(rbind, lapply(df_trait_splits, CWM_inner_function))
    return(out_all)
  }
  
  if(yrs > 1){
    df_yr_splits <- split(df, list(df$Year, df$TraitNameAbr), drop = TRUE)
    out_all <- do.call(rbind, lapply(df_yr_splits, CWM_inner_function))
    return(out_all)
  }
}


# proportions function for categorical traits 

prop_function <- function(df, var){
  yrs = length(levels(as.factor(df$Year)))
  inner_function <- function(dat){
    tot_cov = sum(dat$rel_100)
    var_df = dat[dat[[var]] == 1,]
    dat_out <- data.frame(
      Plot = dat$Plot[1],
      Year = dat$Year[1],
      CWM = sum(var_df$rel_100)/tot_cov,
      TraitNameAbr = var)
    
    return(dat_out)
  }
  if(yrs == 1){
    out <- inner_function(df)
    return(out)
  }
  if(yrs > 1){
    df_yr_splits <- split(df, list(df$Year), drop = TRUE)
    out <- do.call(rbind, lapply(df_yr_splits, inner_function))
    return(out)
  }
}

# ----------------------------------------------------------------------------------------
# ---------------------------- continuous traits  ---------------------------------------
# ----------------------------------------------------------------------------------------

# ------------------------------- prep data  ---------------------------------------------

# select out plots with >80 coverage for traits 
coverage_80 <- coverage[coverage$trait_cov > 80,]
coverage_80_ss <- coverage_80[c("Plot", "Trait", "Year", "trait_cov")]

# merge datasets to select plots from dat that have >80 for specific traits 
dat_80_cov <- merge(dat, coverage_80_ss, by.x = c("Plot", "TraitNameAbr", "Year"), 
                    by.y = c("Plot", "Trait", "Year"), all.y = TRUE)

# ------------  recalculate relative cover with just species that have trait values -------- 

dat_80_splits <- split(dat_80_cov, list(dat_80_cov$Plot, dat_80_cov$Year), drop = TRUE)
spcis_data_relcov_ls <- lapply(dat_80_splits, rel_cov_function)
spcis_data_relcov <- do.call(rbind, spcis_data_relcov_ls)

# drop categorical traits 
dat_cont <- spcis_data_relcov[!spcis_data_relcov$TraitNameAbr %in% c("Duration","Growth.Habit", "Mycorrhizal.type"),]
dat_cont$mean <- as.numeric(dat_cont$mean)

# -------------------- get relativised CWM for cont traits  -----------------------------------

dat_splits <- split(dat_cont, list(dat_cont$Plot), drop = TRUE)

# this takes 10 minutes 
CWM_out <- lapply(dat_splits, CWM_function)
CWM_out_df <- do.call(rbind, CWM_out)

# ----------------------------------------------------------------------------------------
# ---------------------------- categorical traits  ---------------------------------------
# ----------------------------------------------------------------------------------------

# ------------------------------- growth habit ---------------------------------------------

GH_coverage <- coverage[coverage$Trait == "Growth.Habit",]
GH_coverage_80 <- GH_coverage[GH_coverage$trait_cov >= 80,]

dat_GH <- dat[dat$TraitNameAbr == "Growth.Habit",]

dat_GH$Woodiness <- ifelse(dat_GH$mean %in% c("Shrub", "Tree"), 1,
                       ifelse(dat_GH$mean %in% c("Forb", "Graminoid", "Vine"), 0, NA))

dat_GH_80 <- dat_GH[dat_GH$Plot %in% c(GH_coverage_80$Plot),]

dat_GH_80 <- dat_GH_80[!is.na(dat_GH_80$mean),]

dat_GH_80_splits <- split(dat_GH_80, list(dat_GH_80$Plot))

woody_prop_df <- do.call(rbind, lapply(dat_GH_80_splits, prop_function, "Woodiness"))

# ------------------------------- duration ---------------------------------------------

Duration_coverage <- coverage[coverage$Trait == "Duration",]
Duration_coverage_80 <- Duration_coverage[Duration_coverage$trait_cov >= 80,]

dat_Duration <- dat[dat$TraitNameAbr == "Duration",]

dat_Duration$Annual <- ifelse(dat_Duration$mean %in% c("Annual"), 1,
                              ifelse(dat_Duration$mean %in% c("Biennial", "Perennial"), 0, NA))

dat_Duration_80 <- dat_Duration[dat_Duration$Plot %in% c(Duration_coverage_80$Plot),]

dat_Duration_80 <- dat_Duration_80[!is.na(dat_Duration_80$mean),]

dat_Duration_80_splits <- split(dat_Duration_80, list(dat_Duration_80$Plot))

annual_prop_df <- do.call(rbind, lapply(dat_Duration_80_splits, prop_function, "Annual"))

# ------------------------------------- myco ?? -----------------------------------------------

Myco_coverage <- coverage[coverage$Trait == "Mycorrhizal.type",]
Myco_coverage_80 <- Myco_coverage[Myco_coverage$trait_cov >= 80,]

dat_Myco <- dat[dat$TraitNameAbr == "Mycorrhizal.type",]

# -------------------------------- write csv ---------------------------------------------------
CWM_out_all <- do.call(rbind, list(CWM_out_df, woody_prop_df, annual_prop_df))

write.csv(CWM_out_all, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM.csv")
