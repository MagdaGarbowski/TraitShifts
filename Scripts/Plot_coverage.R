# Trait-Shifts 
# goal: get dataset with plot-level trait coverage and plot level I, IN, N coverage 

SPCIS_TRY <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_sumstats.csv")
SPCIS_other <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/other_db_sumstats.csv")
SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022.csv")

# ---------------------- combine other with TRY  --------------------------------------
# get full species list 
SPCIS_all <- rbind(SPCIS_TRY, SPCIS_other)
sps_list <- data.frame(species_matched = unique(SPCIS_all$sps_try_match))

# TRY into wide format 
SPCIS_TRY <- SPCIS_TRY[c("sps_try_match", "TraitNameAbr", "mean", "sd")]
SPCIS_TRY_wide <- reshape(SPCIS_TRY, idvar = "sps_try_match", timevar = "TraitNameAbr", direction = "wide")
SPCIS_TRY_merged <- merge(sps_list, SPCIS_TRY_wide, by.x = "species_matched", by.y = "sps_try_match", all.x = TRUE)

# FUNGAL 
fungal <- SPCIS_other[SPCIS_other$TraitNameAbr == "Mycorrhizal.type",][c("sps_try_match", "mean")]
colnames(fungal)[2] <- "Mycorrhizal.type"

SPCIS_TRY_fungal_merged <- merge(SPCIS_TRY_merged, fungal, by.x = "species_matched", by.y = "sps_try_match", all.x = TRUE)

# GRoot 

groot <- SPCIS_other[!SPCIS_other$TraitNameAbr %in% c("Mycorrhizal.type", "root_depth"),][c("sps_try_match", "mean", "TraitNameAbr")]


# ------------------------------- functions ---------------------------------------------
# relative cover of species plot/year
rel_cov_function <- function(df){
  tot_cover <- sum(df$PctCov)
  df$rel_100 <- (df$PctCov/tot_cover) * 100
  return(df)
}

# relative cover of N, IN, I plot/year 
nat_int_function <- function(df){
  out <- data.frame(Plot = df$Plot[1],
                    Zone = df$Zone[1],
                    Year = df$Year[1],
                    N_relcov = sum(df[df$NativeStatus == "N",][c("rel_100")], na.rm = TRUE),
                    I_relcov = sum(df[df$NativeStatus == "I",][c("rel_100")], na.rm = TRUE),
                    NI_relcov = sum(df[df$NativeStatus == "NI",][c("rel_100")], na.rm = TRUE),
                    UNKstat_relcov = sum(df[is.na(df$NativeStatus),][c("rel_100")]))
  return(out)
}

# trait coverage by plot/year
# for plots that have been sampled in more than one year, these will need to first be split by year
coverage_function <- function(df){
  yrs = length(levels(as.factor(df$Year)))
  df_function <- function(df){
    out <- data.frame(Plot = df$Plot[1],
                      Zone = df$Zone[1],
                      Year = df$Year[1],
                      Trait = df$TraitNameAbr[1],
                      trait_cov = sum(df$rel_100))
  }
  if(yrs == 1){
    out <- df_function(df)
    return(out)
  }
  if(yrs > 1){
    yr_splits <- split(df, df$Year)
    out_yrs <- do.call(rbind, lapply(yr_splits, df_function))
    return(out_yrs)
  }
}

# -------------------------- simplify datasets and run through functions  ----------------------------

# plot data 
plot_data <- unique(SPCIS[c("Dataset", "Site", "Plot", "Original.Plot", "Long", "Lat", "Zone", "EcoRegionLevelI", "Original.Long", "Original.Lat",
                            "FuzzedCoord", "Year", "Month", "PlotArea.m2", "SamplingMethod", "Resampled")])

# subset SPCIS data
spcis_data_ss <- SPCIS[c("Plot", "Zone", "Year", "AcceptedTaxonName", "Original.TaxonName", "NativeStatus", "PctCov", "PctCov_100")]

# get relative cover of species by plot 
spcis_data_plot_splits <- split(spcis_data_ss, list(spcis_data_ss$Plot, spcis_data_ss$Year), drop = TRUE)
spcis_data_relcov_ls <- lapply(spcis_data_plot_splits, rel_cov_function)
spcis_data_relcov <- do.call(rbind, spcis_data_relcov_ls)

#  get relcov of native vs. introduced 
nat_int_df <- do.call(rbind, lapply(spcis_data_relcov_ls, nat_int_function))

# merge SPCIS data with SPCIS_TRY data 
spcis_data_w_try <- merge(spcis_data_relcov, SPCIS_TRY[c("sps_try_match", "TraitNameAbr", "mean")], 
                          by.x = "AcceptedTaxonName", by.y = "sps_try_match", all.x = TRUE)

# split dataset by plot and trait 
spcis_data_w_try_splits <- split(spcis_data_w_try, 
                                 list(spcis_data_w_try$Plot, 
                                      spcis_data_w_try$TraitNameAbr), drop = TRUE)

# get trait coverage by plot 
coverage_out_ls <- lapply(spcis_data_w_try_splits, coverage_function)
coverage_out <- do.call(rbind, coverage_out_ls)

# merge datasets 
coverage_out_2 <- merge(coverage_out, nat_int_df, by = c("Plot", "Zone", "Year"), all.x = TRUE)
coverage_out_3 <- merge(coverage_out_2, plot_data, by = c("Plot", "Zone", "Year"), all.x = TRUE)

write.csv(coverage_out_3, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv")
write.csv(spcis_data_w_try, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wTRY.csv")
