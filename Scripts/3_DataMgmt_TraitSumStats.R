# Trait shifts 
# goal: get species-level and genus-level averages for traits 
# goal: combine trait databases 

TRY_SPCIS <- read.csv ("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_TRY_22398.csv")
fungal_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/fungalroot_SPCIS.csv")
groot_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/groot_SPCIS.csv")
rootdepth_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/rootdepth_SPCIS.csv")
duration_habit_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/duration_growthhabit_SPCIS.csv")

SPCIS_names <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TNRS.csv")

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

max_root_function <- function(df, species_col, value){
  out <- data.frame(sps_try_match = df[[species_col]][1],
                    TraitNameAbr = "max_rooting_depth_m", 
                    mean = max(df[[value]], na.rm = TRUE), 
                    sd = NA,
                    n_studies = nrow(df),
                    n_species = NA)
  return(out)
}


# ---------------------------------------------- TRY  -------------------------------------------------------
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
                                                                               
# TRY - species-level summary statistics 
# split dataset for summary stats 
TRY_sps_trait_splits <- split(TRY_SPCIS_continous_study, 
                              list(TRY_SPCIS_continous_study$spcis_try_match, 
                                   TRY_SPCIS_continous_study$TraitNameAbr), drop = TRUE)

TRY_SPCIS_trait_avgs <- do.call(rbind, lapply(TRY_sps_trait_splits, sum_stats_sps, "spcis_try_match", "TraitNameAbr", "StdValue"))
TRY_SPCIS_trait_avgs$n_species <- 1

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
TRY_all$source <- "TRY"
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
# species-level dataset 
groot_SPCIS <- groot_SPCIS[c("genus_species", "traitName", "meanSpecies", "entriesStudySite")]
colnames(groot_SPCIS) <- c("sps_try_match", "TraitNameAbr", "mean", "n_studies")
groot_SPCIS$sd <- NA
groot_SPCIS$n_species <- 1

# get genus-level estimates 
groot_genus_sp <- as.data.frame(t(sapply(strsplit(groot_SPCIS$sps_try_match, " "), "[")))
colnames(groot_genus_sp) <- c("Genus", "Species")
SPCIS_groot_genus <- cbind(groot_SPCIS, groot_genus_sp)
groot_genus_df <- SPCIS_groot_genus[c("Genus","TraitNameAbr", "mean", "n_studies")]

# split dataset for summary stats (i.e., genus-level estimates)
groot_genus_df_splits <- split(groot_genus_df, list(groot_genus_df$Genus, groot_genus_df$TraitNameAbr), drop = TRUE)
SPCIS_groot_genera_avgs <- do.call(rbind, lapply(groot_genus_df_splits, sum_stats_genus, "Genus","TraitNameAbr", "mean", "n_studies"))

# full SPCIS_TRY dataset 
groot_all <- rbind(groot_SPCIS, SPCIS_groot_genera_avgs)
groot_all$source <- "GRoot"

groot_all[c("mean", "sd")] <- apply(groot_all[c("mean", "sd")], 2, function(x) round(x, 4))

# ------------------------------- rooting depth  ---------------------------------------------
# both species and genus level values in dataset 

rootdepth_SPCIS$TraitNameAbr <- "rootingdepth_m"

# species-level means 
rootdepth_SPCIS$g_s <- gsub(" ", "_", rootdepth_SPCIS$Name_matched)
rootdepth_genus_sps <- data.frame(rootdepth_SPCIS[grep("_", rootdepth_SPCIS$g_s),])
rootdepth_genus_sps$g_s <- NULL

# split dataset for species-level summary stats 
rootdepth_SPCIS_splits <- split(rootdepth_genus_sps, 
                                list(rootdepth_genus_sps$Name_matched), drop = TRUE)

SPCIS_rootdepth_max <- do.call(rbind, lapply(rootdepth_SPCIS_splits, max_root_function, "Name_matched", "Dr"))
SPCIS_rootdepth_max$mean <- gsub(Inf, NA, SPCIS_rootdepth_max$mean)

# genus-level summary stats  
genus_sp <- as.data.frame(t(sapply(strsplit(SPCIS_rootdepth_max$sps_try_match, " "), "[")))
colnames(genus_sp) <- c("Genus", "Species")
SPCIS_depth_genus <- cbind(SPCIS_rootdepth_max, genus_sp)
genus_df <- SPCIS_depth_genus[c("Genus", "TraitNameAbr", "mean", "n_studies")]
genus_df$mean <- as.numeric(genus_df$mean)

# split dataset for summary stats 
genus_df_splits <- split(genus_df, list(genus_df$Genus), drop = TRUE)
SPCIS_depth_genera_avgs <- do.call(rbind, lapply(genus_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_studies"))

# full rooting depth dataset 
depth_all <- rbind(SPCIS_rootdepth_max, SPCIS_depth_genera_avgs)
depth_all$source <- "RootDepth"

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

write.csv(traits_dat,"/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_traits_sumstats.csv" )
write.csv(SPCIS_traits_wide,"/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_traits.csv" )
