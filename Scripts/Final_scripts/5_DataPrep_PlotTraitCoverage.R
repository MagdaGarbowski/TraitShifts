# project: Trait-Shifts 
# objective: trait coverage by plot
# author: Magda Garbowski 
# date: April 27 

library(data.table)
# --------------------------------- data ---------------------------------------------

SPCIS_traits <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/TraitShifts_TraitTable.csv")
SPCIS <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_10272022.csv")
basal_area_plots <- as.data.frame(fread("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_BA_plots_08102023.csv"))
woodiness <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/Woodiness_by_sps.csv")

# --------------------------------- functions ---------------------------------------------

# combine abundances of subspecies to species level 
sum_sp_fun <- function(df){
  out_df <- aggregate(df[,c(7,8)], df[,c(1:4,6)], sum)
  out_df_na <- unique(rbind(df[is.na(df$AcceptedTaxonName),], df[is.na(df$NativeStatus),]))
  out_df_all <- rbind(out_df, out_df_na[c(1:4,6:8)])
  return(out_df_all)
}

# relative cover of species plot/year
rel_cov_function <- function(df){
  tot_cover <- sum(df$PctCov)
  df$rel_100 <- (df$PctCov/tot_cover) * 100
  return(df)
}

# relative cover check 
rel_cov_check_function <-function(df){
  out <- data.frame(Plot = df$Plot[1],
                    Zone = df$Zone[1],
                    Year = df$Year[1],
                    tot_cov = sum(df$rel_100))
  return(out)
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
coverage_function <- function(df){
  df$rel_100 <- ifelse(is.na(df$Value), 0, df$rel_100)
  yrs = length(levels(as.factor(df$Year)))
  df_function <- function(df){
    out <- data.frame(Plot = df$Plot[1],
                      Zone = df$Zone[1],
                      Year = df$Year[1],
                      Trait = df$Trait[1],
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

max_year_function <- function(df){
  max_year = max(df$Year)
  df_out <- df[df$Year == max_year,]
  return(df_out)
}

# --------------------------------- trait coverage by plot ---------------------------------------------

# get unique plots
plot_data <- unique(SPCIS[c("Dataset", "Site", "Plot", "Original.Plot", "Long", "Lat", "Zone", "EcoRegionLevelI", "Original.Long", "Original.Lat",
                            "FuzzedCoord", "Year", "Month", "PlotArea.m2", "SamplingMethod", "Resampled")])

# subset SPCIS data
spcis_data_ss <- SPCIS[c("Plot", "Zone", "Year", "AcceptedTaxonName", "Original.TaxonName", "NativeStatus", "PctCov", "PctCov_100")]

# combing like abundances 
spcis_data_plot_splits <- split(spcis_data_ss, list(spcis_data_ss$Plot, spcis_data_ss$Year), drop = TRUE)
spcis_data_plot_splits_2 <- lapply(spcis_data_plot_splits, sum_sp_fun)

# relative cover of species by plot 
spcis_data_relcov_ls <- lapply(spcis_data_plot_splits_2, rel_cov_function)
spcis_data_relcov <- do.call(rbind, spcis_data_relcov_ls)

# check relative abundances 
rel_ab_check <- do.call(rbind, lapply(spcis_data_relcov_ls, rel_cov_check_function))

# get relative cover of native and introduced species
nat_int_df <- do.call(rbind, lapply(spcis_data_relcov_ls, nat_int_function))

# merge spcis data with trait data 
spcis_data_w_traits <- merge(spcis_data_relcov, SPCIS_traits,
                          by.x = "AcceptedTaxonName", by.y = "Taxon", all.x = TRUE)

# drop genus-level traits 
spcis_data_w_traits_2 <- spcis_data_w_traits[spcis_data_w_traits$NativeStatus %in% c("N", "I"),]

# split dataset by plot and trait 
spcis_data_w_traits_splits <- split(spcis_data_w_traits_2, 
                                 list(spcis_data_w_traits_2$Plot, 
                                      spcis_data_w_traits_2$Trait), drop = TRUE)

# get trait coverage by plot 
coverage_out_ls <- lapply(spcis_data_w_traits_splits, coverage_function)
coverage_out <- do.call(rbind, coverage_out_ls)

# merge datasets 
coverage_out_2 <- merge(coverage_out, nat_int_df, by = c("Plot", "Zone", "Year"), all.x = TRUE)
coverage_out_3 <- merge(coverage_out_2, plot_data, by = c("Plot", "Zone", "Year"), all.x = TRUE)

# drop duplicates 
coverage_out_3  <- coverage_out_3[!duplicated(coverage_out_3),]

# drop basal area plots 
coverage_out_3 <- coverage_out_3[!coverage_out_3$Plot %in% c(basal_area_plots$Plot),]

# get max year 
coverage_out_3_splits <- split(coverage_out_3, coverage_out_3$Plot)
coverage_out_4 <- do.call(rbind, lapply(coverage_out_3_splits, max_year_function))

# ---------------------- trait coverage by plot with woody vs. herbacious ------------------------------
spcis_data_w_traits_woodiness <- merge(spcis_data_w_traits, 
                                       woodiness[c("AcceptedTaxonName", "Woody.Herbacious")], 
                                       by = "AcceptedTaxonName")

# drop species with NA for woody.herbacious 
spcis_data_w_traits_woodiness_2 <- spcis_data_w_traits_woodiness[!is.na(spcis_data_w_traits_woodiness$Woody.Herbacious),]

# split dataset by plot and trait 
spcis_data_w_traits_woodiness_2_splits <- split(spcis_data_w_traits_woodiness_2, 
                                    list(spcis_data_w_traits_woodiness_2$Plot, 
                                         spcis_data_w_traits_woodiness_2$Trait), drop = TRUE)

# get trait coverage by plot 
coverage_out_w_woodiness_ls <- lapply(spcis_data_w_traits_woodiness_2_splits, coverage_function)
coverage_out_w_woodiness <- do.call(rbind, coverage_out_ls)

# merge datasets 
coverage_out_w_woodiness_2 <- merge(coverage_out_w_woodiness, nat_int_df, by = c("Plot", "Zone", "Year"), all.x = TRUE)
coverage_out_w_woodiness_3 <- merge(coverage_out_w_woodiness_2, plot_data, by = c("Plot", "Zone", "Year"), all.x = TRUE)

# drop duplicates  
coverage_out_w_woodiness_3  <- coverage_out_w_woodiness_3[!duplicated(coverage_out_w_woodiness_3),]

# drop basal area plots 
coverage_out_w_woodiness_3 <- coverage_out_w_woodiness_3[!coverage_out_w_woodiness_3$Plot %in% c(basal_area_plots$Plot),]

# # get max year 
coverage_out_w_woodiness_3_splits <- split(coverage_out_w_woodiness_3, coverage_out_3$Plot)
coverage_out_w_woodiness_4 <- do.call(rbind, lapply(coverage_out_w_woodiness_3_splits, max_year_function))


# --------------------------------------- output datasets ---------------------------------------------------

# output datasets 
write.csv(coverage_out_4, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_trait_coverage.csv")
write.csv(spcis_data_w_traits_2, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_wtraits.csv")
write.csv(coverage_out_w_woodiness_4, "/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_trait_coverage_w_woodiness.csv")
