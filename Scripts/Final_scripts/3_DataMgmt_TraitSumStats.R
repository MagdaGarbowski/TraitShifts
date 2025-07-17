# Trait shifts 
# goal: get species-level and genus-level averages for traits 
# goal: combine trait databases 


# Notes to remember what is being done here: 
# averages for traits: first take study-level (dataset in TRY) average, then take average of these averages 
# I am doing this to "match" data from aggregate sources (i.e., literature, GRoot)
# for max height and max rooting depth: 
# ALL values from TRY are used, aggregate values from the literature or GRoot are used and the 97.5 quantile is calculated
# This may be biasing values towards TRY since so many more come from that database 
# n studies for these values reflect number of observations, not studies 

TRY_SPCIS <- read.csv ("/Users/magdagarbowski/TraitShifts/Generated_Data/SPCIS_10272022_TRY_22398.csv")
fungal_SPCIS <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/fungalroot_SPCIS.csv")
groot_SPCIS <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/groot_SPCIS.csv")
rootdepth_SPCIS <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/rootdepth_SPCIS.csv")
duration_habit_SPCIS <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/duration_growthhabit_SPCIS.csv")
literature_data_figures <- read.csv("/Users/magdagarbowski/TraitShifts/Data/TShifts_abundantsps_values_from_figures.csv")
literature_data_tables <- read.csv("/Users/magdagarbowski/TraitShifts/Data/TShifts_abundantsps_values_from_tables.csv")

SPCIS_names <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/SPCIS_TNRS.csv")

# ----------------------------- individual edits of TNRS  ----------------------------------------

SPCIS_names$Accepted_species <- ifelse(SPCIS_names$Name_submitted == "Eriogonum nudum", "Eriogonum nudum", 
                                           ifelse(SPCIS_names$Name_submitted == "Viburnum nudum", "Viburnum nudum", SPCIS_names$Accepted_species))

SPCIS_names$Accepted_name <- ifelse(SPCIS_names$Name_submitted == "Eriogonum nudum", "Eriogonum nudum", 
                                        ifelse(SPCIS_names$Name_submitted == "Viburnum nudum", "Viburnum nudum", SPCIS_names$Accepted_name))

SPCIS_names$Name_matched <- ifelse(SPCIS_names$Name_submitted == "Eriogonum nudum", "Eriogonum nudum", 
                                    ifelse(SPCIS_names$Name_submitted == "Viburnum nudum", "Viburnum nudum", SPCIS_names$Name_matched))

# ------------------------------- functions ---------------------------------------------
sum_stats_study <- function(df){
  StdValue = mean(df$StdValue, na.rm = TRUE)
  out <- cbind(data.frame(df[1,c("SpeciesName", "AccSpeciesName", "TraitID", "TraitName", "UnitName", "DatasetID", "Dataset", "spcis_try_match", "TraitNameAbr")]),StdValue)
  return(out)
}

sum_stats_study_2 <- function(df, cols){
  mean = mean(df$mean, na.rm = TRUE)
  n_obs = nrow(df)
  out <- cbind(data.frame(df[1, cols]), mean, n_obs)
  return(out)
}

sum_stats_sps <- function(df, species_col, trait_col, value){
  out <-  data.frame(sps_try_match = df[[species_col]][1],
                     TraitNameAbr = df[[trait_col]][1], 
                     mean = mean(df[[value]], na.rm = TRUE), 
                     sd = sd(df[[value]], na.rm = TRUE), 
                     n_studies = nrow(df))
  return(out)
}

sum_stats_genus <- function(df, species_col, trait_col, value, n_obs){
  out <-  data.frame(sps_try_match = df[[species_col]][1],
                     TraitNameAbr = df[[trait_col]][1], 
                     mean = mean(df[[value]], na.rm = TRUE), 
                     sd = sd(df[[value]], na.rm = TRUE), 
                     n_studies = sum(df[[n_obs]], na.rm = TRUE),
                     n_species = nrow(df))
  return(out)
}

groot_weighted_avg <- function(df){
  total_studies <- sum(df$n_studies)
  df$weighted_mean <- df$mean * (df$n_studies/total_studies)
  overall_mean <- sum(df$weighted_mean)
  out <- data.frame(sps_try_match = df$sps_try_match[1],
                    TraitNameAbr = df$TraitNameAbr[1],
                    mean = overall_mean,
                    n_studies = total_studies,
                    n_species = 1,
                    sd = NA)
  return(out)
}

max_height_function <- function(df, species_col, value){
  out <- data.frame(sps_try_match = df[[species_col]][1],
                    TraitNameAbr = "max_heightveg_m", 
                    mean = max(df[[value]], na.rm = TRUE), 
                    sd = NA,
                    n_studies = nrow(df),
                    n_species = 1)
  return(out)
}

max_root_function <- function(df, species_col, value){
  out <- data.frame(sps_try_match = df[[species_col]][1],
                    TraitNameAbr = "max_rooting_depth_m", 
                    mean = max(df[[value]], na.rm = TRUE), 
                    sd = NA,
                    n_studies = nrow(df),
                    n_species = 1)
  return(out)
}

max_height_function_2 <- function(df, species_col){
  quantile_975 <- quantile(df[["StdValue"]], probs = c(0.975))
  df_out <- data.frame(sps_try_match = df[[species_col]][1],
                       TraitNameAbr = "max_975_heightveg_m", 
                       mean = quantile_975, 
                       sd = NA,
                       n_studies = nrow(df),
                       n_species = 1)
  return(df_out)
}

max_root_function_2 <- function(df, species_col){
  quantile_975 <- quantile(df[["mean"]], probs = c(0.975))
  df_out <- data.frame(sps_try_match = df[[species_col]][1],
                       TraitNameAbr = "max_975_rootdepth_m", 
                       mean = quantile_975, 
                       sd = NA,
                       n_studies = nrow(df),
                       n_species = 1)
  return(df_out)
}

mean_root_depth_function <- function(df){
  out <- data.frame(sps_try_match = df$sps_try_match[1],
                    TraitNameAbr = "mean_rooting_depth_m", 
                    mean = mean(df$mean, na.rm = TRUE),
                    sd = sd(df$mean, na.rm = TRUE), 
                    n_studies = nrow(df),
                    n_species = 1)
  return(out)
}

# -------------------------- literature data  ---------------------------------------
literature_data_tables <- literature_data_tables[c("species", "TraitNameAbr", "mean", "reference")]
literature_data_figures <- literature_data_figures[c("species", "TraitNameAbr", "mean", "reference")]

literature_data_all <- rbind(literature_data_tables, literature_data_figures)
literature_data_all_splits <- split(literature_data_all, list(literature_data_all$TraitNameAbr, 
                                                              literature_data_all$species,
                                                              literature_data_all$reference), drop = TRUE)

literature_study <- do.call(rbind, lapply(literature_data_all_splits, sum_stats_study_2, c("species", "TraitNameAbr", "reference")))

#### These values need to be integrated with TRY, GRoot, Rooting_depth, and then overall means recalculated
#### Subset literature study by "appropriate" traits, bind, get averages, combine datasets 

TRY_lit <- literature_study[literature_study$TraitNameAbr %in% c("heightveg_m","LDMC_g/g","leafN_mg/g","leafP_mg/g","seedmass_mg", "SLA_mm2/mg","SSD_g/cm3"),]
root_lit <- literature_study[literature_study$TraitNameAbr %in% c("Mean_Root_diameter","Root_N_concentration", "Root_tissue_density", "Specific_root_length"),]
depth_lit <- literature_study[literature_study$TraitNameAbr %in% c("max_rooting_depth_m"),]

# get averages for root_lit to match groot 
root_lit_splits <- split(root_lit, list(root_lit$species, root_lit$TraitNameAbr), drop = TRUE)
root_lit_avg <- do.call(rbind, lapply(root_lit_splits, sum_stats_sps, "species", "TraitNameAbr", "mean"))

# ---------------------------------------------- TRY -------------------------------------------------------
# make columns match other datasets 
colnames(TRY_lit) <- c("spcis_try_match", "TraitNameAbr","Dataset", "StdValue")

# select out continuous traits 
TRY_SPCIS_continous <- TRY_SPCIS[TRY_SPCIS$TraitID %in% c("3106", "26", "3116", "47", "4", "3114", "3117", "14", "3115", "15",
                                                          "3113", "3110", "3111", "3107", "3112", "3108", "3109", "3086"),]

# drop rows with na in StdValue column 
TRY_SPCIS_continous <- TRY_SPCIS_continous[!is.na(TRY_SPCIS_continous$StdValue),]

# shorten trait names 
TRY_SPCIS_continous$TraitNameAbr <- ifelse(TRY_SPCIS_continous$TraitID %in% c(  "3115", "3116", "3117", "3086"), "SLA_mm2/mg",
                                          ifelse(TRY_SPCIS_continous$TraitID == "3106", "heightveg_m", 
                                                 ifelse(TRY_SPCIS_continous$TraitID == "3107", "heightgen_m", 
                                                        ifelse(TRY_SPCIS_continous$TraitID %in% c("3108", "3109", "3110", "3111", "3112","3113","3114"), "leafarea_mm2", 
                                                               ifelse(TRY_SPCIS_continous$TraitID == "26", "seedmass_mg", 
                                                                      ifelse(TRY_SPCIS_continous$TraitID == "14", "leafN_mg/g",
                                                                             ifelse(TRY_SPCIS_continous$TraitID == "15", "leafP_mg/g",
                                                                                    ifelse(TRY_SPCIS_continous$TraitID == "47", "LDMC_g/g", 
                                                                                           ifelse(TRY_SPCIS_continous$TraitID == "4", "SSD_g/cm3", TRY_SPCIS_continous$TraitName)))))))))

# First get study level averages 
TRY_sps_study_trait_splits <- split(TRY_SPCIS_continous, 
                                    list(TRY_SPCIS_continous$spcis_try_match, 
                                         TRY_SPCIS_continous$TraitNameAbr, 
                                         TRY_SPCIS_continous$DatasetID), drop = TRUE)

TRY_SPCIS_continous_study <- do.call(rbind, lapply(TRY_sps_study_trait_splits, sum_stats_study))
TRY_SPCIS_continous_study <- TRY_SPCIS_continous_study[c("spcis_try_match", "TraitNameAbr","Dataset", "StdValue")]

# combine with literature data
TRY_SPCIS_continous_study_w_lit <- rbind(TRY_lit[c("spcis_try_match", "TraitNameAbr","Dataset", "StdValue")], TRY_SPCIS_continous_study)

# TRY - species-level summary statistics 
# split dataset for summary stats 
TRY_sps_trait_splits <- split(TRY_SPCIS_continous_study_w_lit, 
                              list(TRY_SPCIS_continous_study_w_lit$spcis_try_match, 
                                   TRY_SPCIS_continous_study_w_lit$TraitNameAbr), drop = TRUE)

TRY_SPCIS_trait_avgs <- do.call(rbind, lapply(TRY_sps_trait_splits, sum_stats_sps, "spcis_try_match", "TraitNameAbr", "StdValue"))
TRY_SPCIS_trait_avgs$n_species <- 1

# ---------------------------------------------- TRY & literature - MAX height  -------------------------------------------------------

TRY_lit <- rbind(TRY_lit[c("spcis_try_match", "TraitNameAbr","Dataset", "StdValue")], 
                 TRY_SPCIS_continous[c("spcis_try_match", "TraitNameAbr","Dataset", "StdValue")])

TRY_lit_veg_height <- TRY_lit[TRY_lit$TraitNameAbr == "heightveg_m",]

TRY_lit_veg_height_splits <- split(TRY_lit_veg_height,list(TRY_lit_veg_height$spcis_try_match), drop = TRUE)

TRY_lit_veg_height_max <- do.call(rbind, lapply(TRY_lit_veg_height_splits, max_height_function_2, "spcis_try_match"))

# bind max height with species avgs
TRY_SPCIS_trait_avgs <- rbind(TRY_SPCIS_trait_avgs, TRY_lit_veg_height_max)

# TRY - genus-level summary statistics 
TRY_genus_sp <- as.data.frame(t(sapply(strsplit(TRY_SPCIS_trait_avgs$sps_try_match, " "), "[")))
TRY_genus_sp <- t(sapply(strsplit(TRY_SPCIS_trait_avgs$sps_try_match, " "), function(x) rbind(x[1:2])))

colnames(TRY_genus_sp) <- c("Genus", "Species")
TRY_SPCIS_trait_avgs_genus <- cbind(TRY_SPCIS_trait_avgs, TRY_genus_sp)
TRY_genus_df <- TRY_SPCIS_trait_avgs_genus[c("Genus","TraitNameAbr", "mean", "n_studies")]

# split genus-level dataset for summary stats 
TRY_genus_df_splits <- split(TRY_genus_df, list(TRY_genus_df$Genus, TRY_genus_df$TraitNameAbr), drop = TRUE)
TRY_SPCIS_genera_avgs <- do.call(rbind, lapply(TRY_genus_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_studies"))

# full SPCIS_TRY dataset 
TRY_all <- rbind(TRY_SPCIS_trait_avgs, TRY_SPCIS_genera_avgs)
TRY_all$X <- NULL
TRY_all$source <- "TRY or literature"
TRY_all[c("mean", "sd")] <- apply(TRY_all[c("mean", "sd")], 2, function(x) round(x, 4))

# ---------------------------------------------- fungal root  ---------------------------------------------
# merge with SPCIS at genus-level since no species-level data are available

SPCIS_names_ss <- SPCIS_names[c("Name_submitted", "Name_matched", "Genus_matched")]
fungal_SPCIS_matched <- merge(fungal_SPCIS, SPCIS_names_ss, by.x = "Genus", by.y = "Genus_matched", all.y = TRUE)
fungal_SPCIS_matched <- fungal_SPCIS_matched[c("Mycorrhizal.type", "Name_matched")]

# myco type will be in "mean" column to match other datasets 
colnames(fungal_SPCIS_matched) <- c("mean", "sps_try_match")
fungal_SPCIS_matched$TraitNameAbr <- "Mycorrhizal.type"

# add columns for rbind of datasets 
fungal_SPCIS_matched$sd <- NA
fungal_SPCIS_matched$n_studies <- NA
fungal_SPCIS_matched$n_species <- NA
fungal_all <- fungal_SPCIS_matched
fungal_all$source <- "FungalRoot"

# --------------------------------------- duration, growth habit  -----------------------------------------
duration_habit_SPCIS$X <- NULL
DH_SPCIS_matched <- merge(duration_habit_SPCIS, SPCIS_names_ss[c("Name_submitted", "Name_matched")], 
                          by.x = c("AcceptedTaxonName", "Name_matched"),
                          by.y = c("Name_submitted", "Name_matched"), all.y = TRUE)

DH_SPCIS_matched <- DH_SPCIS_matched[c("Name_matched", "Duration", "Growth.Habit")]

# drop duplicates 
DH_SPCIS_matched <- DH_SPCIS_matched[!duplicated(DH_SPCIS_matched$Name_matched),]

DH_SPCIS_long <- reshape(DH_SPCIS_matched,
                         varying = c("Duration", "Growth.Habit"),
                         idvar = "Name_matched",
                         timevar = "TraitNameAbr", 
                         v.names = "mean",
                         times = c("Duration", "Growth.Habit"),
                         direction = "long")

colnames(DH_SPCIS_long)[1] <- "sps_try_match"

DH_SPCIS_long$mean <- ifelse(grepl(",", DH_SPCIS_long$mean), gsub(",.*", "", DH_SPCIS_long$mean), DH_SPCIS_long$mean)
DH_SPCIS_long$mean <- gsub("Forb/herb", "Forb", DH_SPCIS_long$mean)
DH_SPCIS_long$mean <- gsub("Subshrub", "Shrub", DH_SPCIS_long$mean)
                   
DH_SPCIS_long$sd <- NA
DH_SPCIS_long$n_studies <- NA
DH_SPCIS_long$n_species <- NA
DH_all <- DH_SPCIS_long
DH_all$source <- "SPCIS_data"

# remove attributes
DH_all <- lapply(DH_all, unname)
DH_all <- as.data.frame(DH_all)
# ---------------------------------------------- groot  ---------------------------------------------------

# drop rooting depth <- to be combined with rooting depth data 
groot_SPCIS_nodepth <- groot_SPCIS[!groot_SPCIS$traitName == "Rooting_depth",]

# species-level dataset 
groot_SPCIS_nodepth <- groot_SPCIS_nodepth[c("genus_species", "traitName", "meanSpecies", "entriesStudySite")]
colnames(groot_SPCIS_nodepth) <- c("sps_try_match", "TraitNameAbr", "mean", "n_studies")
groot_SPCIS_nodepth$sd <- NA

groot_SPCIS_w_lit <- rbind(groot_SPCIS_nodepth, root_lit_avg)
groot_SPCIS_w_lit$n_species <- 1

# need a weighted average accounting for number of studies in groot vs root_lit
groot_SPCIS_w_lit_splits <- split(groot_SPCIS_w_lit, list(groot_SPCIS_w_lit$sps_try_match, groot_SPCIS_w_lit$TraitNameAbr), drop = TRUE)

groot_lit_all <- do.call(rbind, lapply(groot_SPCIS_w_lit_splits, groot_weighted_avg))

# get genus-level estimates 
groot_genus_sp <- as.data.frame(t(sapply(strsplit(groot_lit_all$sps_try_match, " "), "[")))
colnames(groot_genus_sp) <- c("Genus", "Species")
SPCIS_groot_genus <- cbind(groot_lit_all, groot_genus_sp)
groot_genus_df <- SPCIS_groot_genus[c("Genus","TraitNameAbr", "mean", "n_studies")]

# split dataset for summary stats (i.e., genus-level estimates)
groot_genus_df_splits <- split(groot_genus_df, list(groot_genus_df$Genus, groot_genus_df$TraitNameAbr), drop = TRUE)
SPCIS_groot_genera_avgs <- do.call(rbind, lapply(groot_genus_df_splits, sum_stats_genus, "Genus","TraitNameAbr", "mean", "n_studies"))

# full GRoot dataset 
groot_all <- rbind(groot_lit_all, SPCIS_groot_genera_avgs)
groot_all$source <- "GRoot or literature"

groot_all[c("mean", "sd")] <- apply(groot_all[c("mean", "sd")], 2, function(x) round(x, 4))

# ------------------------------- rooting depth  ---------------------------------------------

# groot rooting depth  
groot_depth <- groot_SPCIS[groot_SPCIS$traitName == "Rooting_depth",]
groot_depth <- groot_depth[c("genus_species", "traitName", "meanSpecies", "entriesStudySite")]
colnames(groot_depth) <- c("sps_try_match", "TraitNameAbr", "mean", "n_studies")
groot_depth$TraitNameAbr <- "rootingdepth_m"
groot_depth$sd <- NA

# literature rooting depth 
depth_lit <- depth_lit[c("species", "TraitNameAbr", "mean", "n_obs")]
colnames(depth_lit) <- c("sps_try_match", "TraitNameAbr", "mean", "n_studies")
depth_lit$sd <- NA

# RSIP rooting depth  
rootdepth_SPCIS$TraitNameAbr <- "rootingdepth_m"
rootdepth_SPCIS$g_s <- gsub(" ", "_", rootdepth_SPCIS$Name_matched)
rootdepth_genus_sps <- data.frame(rootdepth_SPCIS[grep("_", rootdepth_SPCIS$g_s),])
rootdepth_genus_sps <- rootdepth_genus_sps[c("g_s", "Dr", "TraitNameAbr")]
rootdepth_genus_sps$g_s <- gsub("_", " ", rootdepth_genus_sps$g_s)
colnames(rootdepth_genus_sps) <- c("sps_try_match", "mean", "TraitNameAbr")
rootdepth_genus_sps$sd <- NA
rootdepth_genus_sps$n_studies <- NA

# bind together root depth datasets 
root_depth_all <- do.call(rbind, list(groot_depth, depth_lit, rootdepth_genus_sps))

# remove NAs 
root_depth_all <- root_depth_all[!is.na(root_depth_all$mean),]
  
# split RSIP dataset for species-level summary stats 
rootdepth_SPCIS_splits <- split(root_depth_all, 
                                list(root_depth_all$sps_try_match), drop = TRUE)

# ------- max rooting depth ------- #
SPCIS_rootdepth_max <- do.call(rbind, lapply(rootdepth_SPCIS_splits, max_root_function_2, "sps_try_match"))

# genus-level summary stats - max rooting depth  
genus_sp <- as.data.frame(t(sapply(strsplit(SPCIS_rootdepth_max$sps_try_match, " "), "[")))
colnames(genus_sp) <- c("Genus", "Species")
SPCIS_depth_genus <- cbind(SPCIS_rootdepth_max, genus_sp)
genus_df <- SPCIS_depth_genus[c("Genus", "TraitNameAbr", "mean", "n_studies")]
genus_df$mean <- as.numeric(genus_df$mean)

# split dataset for summary stats 
genus_df_splits <- split(genus_df, list(genus_df$Genus), drop = TRUE)
SPCIS_depth_genera_avgs <- do.call(rbind, lapply(genus_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_studies"))

# ------- mean rooting depth ------- #
SPCIS_rootdepth_mean <- do.call(rbind, lapply(rootdepth_SPCIS_splits, mean_root_depth_function))
# genus-level summary stats - mean rooting depth  
genus_sp_mean <- as.data.frame(t(sapply(strsplit(SPCIS_rootdepth_mean$sps_try_match, " "), "[")))
colnames(genus_sp_mean) <- c("Genus", "Species")
SPCIS_depth_mean_genus <- cbind(SPCIS_rootdepth_mean, genus_sp)
genus_mean_df <- SPCIS_depth_mean_genus[c("Genus", "TraitNameAbr", "mean", "n_studies")]
genus_mean_df$mean <- as.numeric(genus_mean_df$mean)

# split dataset for summary stats 
genus_mean_df_splits <- split(genus_mean_df, list(genus_mean_df$Genus), drop = TRUE)
SPCIS_depth_mean_genera_avgs <- do.call(rbind, lapply(genus_mean_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_studies"))

# full rooting depth dataset 
depth_all <- do.call(rbind, list(SPCIS_rootdepth_max, SPCIS_depth_genera_avgs, SPCIS_rootdepth_mean, SPCIS_depth_mean_genera_avgs))
depth_all$source <- "GRoot, RSIP, literature"

depth_all[c("mean", "sd")] <- apply(depth_all[c("mean", "sd")], 2, as.numeric)

depth_all[c("mean", "sd")] <- apply(depth_all[c("mean", "sd")], 2, function(x) round(x, 4))

# ---------------------- combine datasets (long format) -----------------------------

traits_dat <- rbind(TRY_all, DH_all,fungal_all, groot_all, depth_all)
traits_dat <- traits_dat[!duplicated(traits_dat),]

# ---------------------- create wide format dataset ---------------------------------

# get full SPCIS list for merging 
SPCIS_species <- data.frame(Species_name = SPCIS_names$Name_matched)

# get continuous traits into wide format 
traits_cont <- traits_dat[c("sps_try_match", "TraitNameAbr", "mean")][!traits_dat$TraitNameAbr %in% c("Mycorrhizal.type", "Duration", "Growth.Habit"),]
traits_cont$mean <- as.numeric(traits_cont$mean)
traits_cont$mean <- round(traits_cont$mean, 4)

traits_cont_wide <- reshape(traits_cont, idvar = "sps_try_match", timevar = "TraitNameAbr", direction = "wide")
colnames(traits_cont_wide) <- gsub("mean.", "", colnames(traits_cont_wide))

# get categorical traits into wide format
traits_cat <- traits_dat[c("sps_try_match","TraitNameAbr", "mean")][traits_dat$TraitNameAbr %in% c("Mycorrhizal.type", "Duration", "Growth.Habit"),]
traits_cat_wide <- reshape(traits_cat, idvar = "sps_try_match", timevar = "TraitNameAbr", direction = "wide")
colnames(traits_cat_wide)[c(2,3,4)] <- c("Duration", "Growth.Habit", "Mycorrhizal.type") 

# merge datasets 
SPCIS_traits_wide <- merge(SPCIS_species, traits_cont_wide, by.x = "Species_name", by.y = "sps_try_match", all.x = TRUE)
SPCIS_traits_wide <- merge(traits_cat_wide, SPCIS_traits_wide, by.x = "sps_try_match", by.y = "Species_name", all.x = TRUE)

# drop duplicates 
SPCIS_traits_wide <- SPCIS_traits_wide[!duplicated(SPCIS_traits_wide),]

write.csv(traits_dat, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_traits_sumstats.csv" )
write.csv(SPCIS_traits_wide, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_traits.csv" )
