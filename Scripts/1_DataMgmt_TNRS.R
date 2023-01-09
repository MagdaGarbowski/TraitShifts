# Trait Shifts - data management (1)

# goal: clean and run datasets (TRY, fungal_root, groot, rooting_depth, SPCIS) through TNRS for matching
# note: this script takes over an hour to run
# note: after being cleaned, TRY data needs to be broken up into batches to run through TNRS
# resulting datasets: 
# (1) fungalroot, groot, root_depth datasets of SPCIS species 
# (2) simplified and cleaned TRY dataset (traits only with simplified AccSpeciesNames and SpeciesNames)
# (3) TRY species list from (2) with associated TNRS matched names
# (4) cleaned SPCIS dataset (i.e., varieties, subspecies, spaces dropped from "AcceptedTaxonName")
# (5) SPCIS "AcceptedTaxonName" species list from (4) with associated TNRS matched names 

# ------------------------------- packages ---------------------------------------------------

install.packages("remotes")
remotes::install_github("EnquistLab/RTNRS")
library(data.table)
library(TNRS)
library(stringr)

# --------------------------------- data ------------------------------------------------------

TRY_data <- fread("/Users/MagdaGarbowski 1/TraitShifts/Data/TRY_22398.txt", 
                  select = c("SpeciesName", "AccSpeciesName", "AccSpeciesID", "ObservationID",
                             "ObsDataID","TraitID", "TraitName", "OrigValueStr", "StdValue", "UnitName", "ErrorRisk"), quote = "")

fungalroot_data <- fread("/Users/MagdaGarbowski 1/TraitShifts/Data/fungal_root.csv") # tab three of full database (nph16569-sup-0002-tabless1-s4.xlsx)

groot_data <- fread("/Users/MagdaGarbowski 1/TraitShifts/Data/GRooTAggregateSpeciesVersion.csv")

rootdepth_data <- fread("/Users/MagdaGarbowski 1/TraitShifts/Data/rooting_depth.csv", 
                        select = c("ID", "Species", "Dr"))

spcis_data <- fread("/Users/MagdaGarbowski 1/TraitShifts/Data/FULLDatabase_10272022.csv")

# ------------------------------- functions for TNRS -------------------------------------------

# data prep for TNRS
TNRS_data_prep <- function(df, colname){
  df_2 = unique(df[colname])
  df_out <- data.frame(ID = 1:nrow(df_2), colname = df_2[1])
  return(df_out)
}

# TNRS function 
TNRS_function <- function(df) {
  dat = df
  df_out <- TNRS(dat,sources = c("usda", "tropicos", "wcvp", "wfo"), mode = "resolve", matches = "best")
}

# write csv function for TRY batches 
write_csv_function <- function(df, filename){
  id = df$ID[1]
  write.csv(df, file = paste0("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TNRS_TRY/",filename,id,".csv"))
}

# ------------------------------- SPCIS data cleanup ---------------------------------------------
# drop nonvascular and lichenous 
spcis_data <- as.data.frame(spcis_data[!spcis_data$Growth.Habit %in% c("Nonvascular", "Lichenous"),])

# drop varieties and subspecies
spcis_data$AcceptedTaxonName <- gsub(" var.*| ssp.*", "", spcis_data$AcceptedTaxonName)

# trim white space before and after entries
spcis_data$AcceptedTaxonName <- str_trim(spcis_data$AcceptedTaxonName, "both")

# prep for TNRS 
spcis_names_prepped <- TNRS_data_prep(spcis_data, "AcceptedTaxonName")

# SPCIS 
# takes about 5 minutes 
spcis_names_TNRS  <- TNRS(spcis_names_prepped,
                          sources = c("usda", "tropicos", "wcvp", "wfo"), mode = "resolve", matches = "best")

# ----------------------------- SPCIS Duration and Growth.Habit  --------------------------------------

# merge SPCIS names with TNRS names 
spcis_names_TNRS_ss <- spcis_names_TNRS[c("Name_submitted", "Name_matched")]
spcis_data_w_matched_names <- merge(spcis_data, spcis_names_TNRS_ss, by.x = "AcceptedTaxonName", "Name_submitted", all.x = TRUE)

spcis_duration_habit <- spcis_data_w_matched_names[c("AcceptedTaxonName", "Name_matched", "Duration", "Growth.Habit")]
spcis_duration_habit <- unique(spcis_duration_habit)

# ------------------------------- fungal root data cleanup ---------------------------------------------
fungalroot_data <- as.data.frame(fungalroot_data)
fungalroot_prepped <- TNRS_data_prep(fungalroot_data, "Genus")

# run through TNRS
# takes about 10 minutes 
fungalroot_names_TNRS <- TNRS(fungalroot_prepped, 
                              sources = c("usda", "tropicos", "wcvp", "wfo"), mode = "resolve", matches = "best")
# find and drop mismatches 
fungalroot_mismatches <- fungalroot_names_TNRS$Name_submitted[!fungalroot_names_TNRS$Name_submitted %in% c(fungalroot_names_TNRS$Name_matched)]

fungalroot_data <- fungalroot_data[!fungalroot_data$Genus %in% fungalroot_mismatches,]

# get SPCIS species 
fungalroot_SPCIS <- fungalroot_data[fungalroot_data$Genus %in% spcis_names_TNRS$Genus_matched,]

# ------------------------------- groot data cleanup ---------------------------------------------
groot_data <- as.data.frame(groot_data)
groot_data$genus_species <- paste(groot_data$genusTNRS, groot_data$speciesTNRS, sep = " ")

# prep for TNRS 
groot_data_prepped <- TNRS_data_prep(groot_data, "genus_species")

# run through TNRS
# takes about 5 minutes 
groot_names_TNRS <- TNRS(groot_data_prepped, sources = c("usda", "tropicos", "wcvp", "wfo"), mode = "resolve", matches = "best")

# find and drop mismatches 
groot_mismatches <- groot_names_TNRS$Name_submitted[!groot_names_TNRS$Name_submitted %in% c(groot_names_TNRS$Name_matched)]
groot_data <- groot_data[!groot_data$genus_species %in% groot_mismatches,]

# get SPCIS species 
groot_data_SPCIS <- groot_data[groot_data$genus_species %in% spcis_names_TNRS$Name_matched,]

# ------------------------------- root depth cleanup ---------------------------------------------
rootdepth_data <- as.data.frame(rootdepth_data)

# drop ssp., var., Var, sp., L 
rootdepth_data$Species <- gsub(" var.*| Var.*| ssp.*| sp.*| L.*", "", rootdepth_data$Species)

# prep for TNRS 
rootdepth_data_prepped <- TNRS_data_prep(rootdepth_data, "Species")

# run through TNRS 
# takes about 2 minutes 
rootdepth_names_TNRS <- TNRS(rootdepth_data_prepped, sources = c("usda", "tropicos", "wcvp", "wfo"), mode = "resolve", matches = "best")

# find mismatchs 
# do not drop these mismatches - most look like mismatches in spelling 
rootdepth_mismatches <- rootdepth_names_TNRS[!rootdepth_names_TNRS$Name_submitted %in% c(rootdepth_names_TNRS$Name_matched),]

# merge TNRS matches with original root_depth dataset 
rootdepth_data_matched <- merge(rootdepth_data, rootdepth_names_TNRS[c("Name_submitted", "Name_matched")],
                                by.x = "Species", by.y = "Name_submitted")

# get SPCIS species 
rootdepth_data_SPCIS <- rootdepth_data_matched[rootdepth_data_matched$Name_matched %in% spcis_names_TNRS$Name_matched,]

# ------------------------------- TRY clean up  ---------------------------------------------
# drop all entries that are not traits 
TRY_traits <- as.data.frame(TRY_data[!is.na(TRY_data$TraitID),])

# remove obvious mismatches 
TRY_ss <- TRY_traits[!TRY_traits$SpeciesName %in% c("", "??","?", "-"),]
TRY_ss <- TRY_ss[!grepl("[0-9]|\\(", TRY_ss$SpeciesName),] 

# drop values > 4 std away from mean for each trait based on error risk
TRY_ss <- TRY_ss[!TRY_ss$ErrorRisk > 3.99,]

# drop subspecies and varieties from  SpeciesName and AccSpeciesName columns
TRY_ss[c("SpeciesName", "AccSpeciesName")] <- apply(TRY_ss[c("SpeciesName", "AccSpeciesName")], 2, 
                                                    function(x) gsub(" var.*| Var.*| subvar.*| subsp.*| SUBSP*| SSP*| spp.*| ssp.*| ssp*| Sp.*| sp.", "", x))

# add latin encoding so str_trim works  
TRY_ss[c("SpeciesName", "AccSpeciesName")] <- apply(TRY_ss[c("SpeciesName", "AccSpeciesName")], 2, 
                                                    function(x) {Encoding(x) <- "latin1"; return(x)})

# trim white space before and after entries - does not work in a function?
TRY_ss$SpeciesName <- str_trim(TRY_ss$SpeciesName)
TRY_ss$AccSpeciesName <- str_trim(TRY_ss$AccSpeciesName)

# prep data for TNRS for both SpeciesName and AccSpecies name columns
# will need to match on both of these columns 
TRY_SpeciesNames_prepped <- TNRS_data_prep(TRY_ss, "SpeciesName")
TRY_AccSpeciesNames_prepped <- TNRS_data_prep(TRY_ss, "AccSpeciesName")

# split up into smaller datasets to run through TNRS 
TRY_SpeciesNames_prepped_ls <-  split(TRY_SpeciesNames_prepped, rep(1:76, each=1000))
TRY_AccSpeciesNames_prepped_ls <-  split(TRY_AccSpeciesNames_prepped, rep(1:47, each=1000))

# run smaller TRY datasets through TNRS 
# output files into folder in working directory 
# this takes over an hour! 
lapply(lapply(TRY_SpeciesNames_prepped_ls[1:25], TNRS_function), write_csv_function, "TRY_SpNames_")
lapply(lapply(TRY_SpeciesNames_prepped_ls[26:50], TNRS_function), write_csv_function,"TRY_SpNames_")
lapply(lapply(TRY_SpeciesNames_prepped_ls[51:76], TNRS_function), write_csv_function,"TRY_SpNames_")
lapply(lapply(TRY_AccSpeciesNames_prepped_ls[1:25], TNRS_function), write_csv_function,"TRY_AccSpNames_")
lapply(lapply(TRY_AccSpeciesNames_prepped_ls[26:47], TNRS_function), write_csv_function,"TRY_AccSpNames_")

# --------------------------------- bind TRY TNRS files  ------------------------------------------
# set file directory with TRY TNRS files 
filedirectory = "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TNRS_TRY/"

# pull AccSpeciesNames TRY files from working directory 
file_list_Acc <- list.files(path = filedirectory, pattern="AccSpNames", full.names = TRUE)
TryAccSps_ls <- lapply(file_list_Acc, read.csv, stringsAsFactor = FALSE)
TryAccSps_ls_2 <- do.call(rbind, TryAccSps_ls)
TRY_AccSpeciesNames_TNRS <- (TryAccSps_ls_2)

# pull SpeciesNames files from working directory 
file_list_Sps <- list.files(path = filedirectory, pattern="_SpNames_", full.names = TRUE)
TrySps_ls <- lapply(file_list_Sps, read.csv, stringsAsFactor = FALSE)
TrySps_ls_2 <- do.call(rbind, TrySps_ls)
TRY_SpeciesNames_TNRS <- (TrySps_ls_2)

# ------------------------------- write datasets  ---------------------------------------------
write.csv(fungalroot_SPCIS, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/fungalroot_SPCIS.csv")
write.csv(groot_data_SPCIS, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/groot_SPCIS.csv")
write.csv(rootdepth_data_SPCIS, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/rootdepth_SPCIS.csv")
write.csv(spcis_duration_habit, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/duration_growthhabit_SPCIS.csv")

write.csv(TRY_ss, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_TraitsOnly_22398.csv")
write.csv(TRY_SpeciesNames_TNRS,  "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_SpeciesName_TNRS_22398.csv")
write.csv(TRY_AccSpeciesNames_TNRS,  "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_AccSpeciesName_TNRS_22398.csv")
write.csv(spcis_data, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022.csv")
write.csv(spcis_names_TNRS, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TNRS.csv")

