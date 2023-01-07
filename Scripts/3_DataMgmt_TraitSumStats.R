# Trait shifts 
# goal: get species-level and genus-level averages for traits 
# goal: combine trait databases 

TRY_SPCIS <- read.csv ("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_TRY_22398.csv")
fungal_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/fungalroot_SPCIS.csv")
groot_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/groot_SPCIS.csv")
rootdepth_SPCIS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/rootdepth_SPCIS.csv")
SPCIS_names <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TNRS.csv")

# ------------------------------- functions ---------------------------------------------
sum_stats_sps <- function(df, species_col, trait_col, value){
  out <-  data.frame(sps_try_match = df[[species_col]][1],
                     TraitNameAbr = df[[trait_col]][1], 
                     mean = mean(df[[value]], na.rm = TRUE), 
                     sd = sd(df[[value]], na.rm = TRUE), 
                     n_obs = nrow(df))
  return(out)
}

sum_stats_genus <- function(df, species_col, trait_col, value, n_obs){
  out <-  data.frame(sps_try_match = df[[species_col]][1],
                     TraitNameAbr = df[[trait_col]][1], 
                     mean = mean(df[[value]], na.rm = TRUE), 
                     sd = sd(df[[value]], na.rm = TRUE), 
                     n_obs = sum(df[[n_obs]], na.rm = TRUE),
                     n_species = nrow(df))
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
# TRY - species-level summary statistics 
# split dataset for summary stats 
TRY_sps_trait_splits <- split(TRY_SPCIS_continous, 
                              list(TRY_SPCIS_continous$spcis_try_match, 
                                   TRY_SPCIS_continous$TraitNameAbr), drop = TRUE)

TRY_SPCIS_trait_avgs <- do.call(rbind, lapply(TRY_sps_trait_splits, sum_stats_sps, "spcis_try_match", "TraitNameAbr", "StdValue"))
TRY_SPCIS_trait_avgs$n_species <- 1

# TRY - genus-level summary statistics 
TRY_genus_sp <- as.data.frame(t(sapply(strsplit(TRY_SPCIS_trait_avgs$sps_try_match, " "), "[")))
colnames(TRY_genus_sp) <- c("Genus", "Species")
TRY_SPCIS_trait_avgs_genus <- cbind(TRY_SPCIS_trait_avgs, TRY_genus_sp)
TRY_genus_df <- TRY_SPCIS_trait_avgs_genus[c("Genus","TraitNameAbr", "mean", "n_obs")]

# split genus-level dataset for summary stats 
TRY_genus_df_splits <- split(TRY_genus_df, list(TRY_genus_df$Genus, TRY_genus_df$TraitNameAbr), drop = TRUE)
TRY_SPCIS_genera_avgs <- do.call(rbind, lapply(TRY_genus_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_obs"))

# full SPCIS_TRY dataset 
TRY_all <- rbind(TRY_SPCIS_trait_avgs, TRY_SPCIS_genera_avgs)
TRY_all$X <- NULL
TRY_all$source <- "TRY"

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
fungal_SPCIS_matched$n_obs <- NA
fungal_SPCIS_matched$n_species <- NA
fungal_all <- fungal_SPCIS_matched
fungal_all$source <- "FungalRoot"

# ---------------------------------------------- groot  ---------------------------------------------------
# species-level dataset 
groot_SPCIS <- groot_SPCIS[c("genus_species", "traitName", "meanSpecies")]
colnames(groot_SPCIS) <- c("sps_try_match", "TraitNameAbr", "mean")
groot_SPCIS$sd <- NA
groot_SPCIS$n_obs <- NA
groot_SPCIS$n_species <- 1

# get genus-level estimates 
groot_genus_sp <- as.data.frame(t(sapply(strsplit(groot_SPCIS$sps_try_match, " "), "[")))
colnames(groot_genus_sp) <- c("Genus", "Species")
SPCIS_groot_genus <- cbind(groot_SPCIS, groot_genus_sp)
groot_genus_df <- SPCIS_groot_genus[c("Genus","TraitNameAbr", "mean", "n_obs")]

# split dataset for summary stats (i.e., genus-level estimates)
groot_genus_df_splits <- split(groot_genus_df, list(groot_genus_df$Genus, groot_genus_df$TraitNameAbr), drop = TRUE)
SPCIS_groot_genera_avgs <- do.call(rbind, lapply(groot_genus_df_splits, sum_stats_genus, "Genus","TraitNameAbr", "mean", "n_obs"))

# full SPCIS_TRY dataset 
groot_all <- rbind(groot_SPCIS, SPCIS_groot_genera_avgs)
groot_all$n_obs <- NA
groot_all$source <- "GRoot"

# ------------------------------- rooting depth  ---------------------------------------------
# both species and genus level values in dataset 

rootdepth_SPCIS$TraitNameAbr <- "rootingdepth_m"

# get full species dataset for species-level means 
rootdepth_SPCIS$g_s <- gsub(" ", "_", rootdepth_SPCIS$Name_matched)
rootdepth_genus_sps <- data.frame(rootdepth_SPCIS[grep("_", rootdepth_SPCIS$g_s),])
rootdepth_genus_sps$g_s <- NULL

# split dataset for species-level summary stats 
rootdepth_SPCIS_splits <- split(rootdepth_genus_sps, 
                                list(rootdepth_genus_sps$Name_matched), drop = TRUE)

SPCIS_rootdepth_avgs <- do.call(rbind, lapply(rootdepth_SPCIS_splits, sum_stats_sps, "Name_matched","TraitNameAbr", "Dr"))
SPCIS_rootdepth_avgs$n_species <- 1

# get genus-level dataset for genus-level summary stats  
genus_sp <- as.data.frame(t(sapply(strsplit(SPCIS_rootdepth_avgs$sps_try_match, " "), "[")))
colnames(genus_sp) <- c("Genus", "Species")
SPCIS_depth_genus <- cbind(SPCIS_rootdepth_avgs, genus_sp)
genus_df <- SPCIS_depth_genus[c("Genus", "TraitNameAbr", "mean", "n_obs")]

# split dataset for summary stats 
genus_df_splits <- split(genus_df, list(genus_df$Genus), drop = TRUE)
SPCIS_depth_genera_avgs <- do.call(rbind, lapply(genus_df_splits, sum_stats_genus, "Genus", "TraitNameAbr", "mean", "n_obs"))

# full rooting depth dataset 
depth_all <- rbind(SPCIS_rootdepth_avgs, SPCIS_depth_genera_avgs)
depth_all$source <- "RootDepth"

# --------------------------- combine datasets ------------------------------------

traits_dat <- rbind(TRY_all, fungal_all, groot_all, depth_all)

write.csv(traits_dat,"/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/traits_sumstats.csv" )

