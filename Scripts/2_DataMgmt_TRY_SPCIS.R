# Trait Shifts - data management (2)
# goal: merge cleaned and simplified TRY data with cleaned SPCIS data base on both AccSpeciesNames and SpeciesNames in TRY
# resulting dataset: SPCIS_10272022_TRY_22398.csv - TRY traits for SPCIS species 

# --------------------------------- data ------------------------------------------------------

TRY_ss  <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_TraitsOnly_22398.csv")
spcis_data <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022.csv")
TRY_SpeciesNames_list_TNRS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_SpeciesName_TNRS_22398.csv")
TRY_AccSpeciesName_list_TNRS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/TRY_AccSpeciesName_TNRS_22398.csv")
spcis_list_TNRS <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TNRS.csv")

# ----------------------------- subset TNRS dataframes  ----------------------------------------

out_TNRS_dfs <- lapply(list(TRY_SpeciesNames_list_TNRS = TRY_SpeciesNames_list_TNRS,
                            TRY_AccSpeciesName_list_TNRS = TRY_AccSpeciesName_list_TNRS,
                            spcis_list_TNRS = spcis_list_TNRS),
                       function(x) {xx = x[c("ID", "Name_submitted","Accepted_name", "Accepted_species")]; return(xx)})

list2env(out_TNRS_dfs, .GlobalEnv)

# ----------------------------- individual edits of TNRS  ----------------------------------------

spcis_list_TNRS$Accepted_species <- ifelse(spcis_list_TNRS$Name_submitted == "Eriogonum nudum", "Eriogonum nudum", 
                                           ifelse(spcis_list_TNRS$Name_submitted == "Viburnum nudum", "Viburnum nudum", spcis_list_TNRS$Accepted_species))

spcis_list_TNRS$Accepted_name <- ifelse(spcis_list_TNRS$Name_submitted == "Eriogonum nudum", "Eriogonum nudum", 
                                           ifelse(spcis_list_TNRS$Name_submitted == "Viburnum nudum", "Viburnum nudum", spcis_list_TNRS$Accepted_name))

# ------------------------------- merge datasets -------------------------------------------------
# merge TRY simplified with TRY TNRS SpeciesNames
TRY_ss_SpeciesName_merge <- merge(TRY_ss, TRY_SpeciesNames_list_TNRS, 
                                  by.x = "SpeciesName", by.y = "Name_submitted", all.x = TRUE)

# now merge with TRY TNRS AccSpeciesNames
TRY_ss_SpeciesName_AccName_merge <- merge(TRY_ss_SpeciesName_merge, TRY_AccSpeciesName_list_TNRS, 
                                          by.x = "AccSpeciesName", by.y = "Name_submitted", all.x = TRUE)
# subset relevant columns
TRY_ss_SpeciesName_AccName <- TRY_ss_SpeciesName_AccName_merge[c("SpeciesName","AccSpeciesName",  "TraitID", "TraitName", "OrigValueStr", "StdValue", "UnitName",
                                                                 "ObservationID", "ObsDataID", "Accepted_species.x", "Accepted_species.y")]
# rename "Accepted_species" columns 
colnames(TRY_ss_SpeciesName_AccName)[10:11] <- c("Accepted_species.species",  "Accepted_species.Acc")

# drop blanks
TRY_ss_SpeciesName_AccName <- TRY_ss_SpeciesName_AccName[!(TRY_ss_SpeciesName_AccName$Accepted_species.Acc == "" & 
                                                             TRY_ss_SpeciesName_AccName$Accepted_species.species == ""),]

# keep rows in which either Accepted_species.species or Accepted_species.Acc matches SPCIS Accepted_species
TRY_SPCIS <- TRY_ss_SpeciesName_AccName[(TRY_ss_SpeciesName_AccName$Accepted_species.species %in% c(spcis_list_TNRS$Accepted_species)) |
                                    (TRY_ss_SpeciesName_AccName$Accepted_species.Acc %in% c(spcis_list_TNRS$Accepted_species)),]  

# get unique names 
TRY_SPCIS_unique <- as.data.frame(unique(TRY_SPCIS[c("SpeciesName", "AccSpeciesName","Accepted_species.species",  "Accepted_species.Acc")]))

# add na to blanks 
TRY_SPCIS_unique$Accepted_species.species <- ifelse(TRY_SPCIS_unique$Accepted_species.species == "", NA, TRY_SPCIS_unique$Accepted_species.species)
TRY_SPCIS_unique$Accepted_species.Acc <- ifelse(TRY_SPCIS_unique$Accepted_species.Acc == "", NA, TRY_SPCIS_unique$Accepted_species.Acc)

# function to select match in SPCIS from TRY "Accepted_species.species" or "Accepted_species.Acc" 
select_function <- function(df){
  df$try_spcis_match <- ifelse(df$Accepted_species.species[1] %in% spcis_list_TNRS$Accepted_species, 
                               paste(df$Accepted_species.species),
                               paste(df$Accepted_species.Acc))
  colnames(df) <- c("SpeciesName", "AccSpeciesName", "Accepted_species.species", "Accepted_species.Acc", "spcis_try_match")
  return(df)
}

# split mismatches TRY_SPCIS into list by rows to run select_function 
TRY_SPCIS_list <- split(TRY_SPCIS_unique, rownames(TRY_SPCIS_unique))
TRY_SPCIS_list_out <- do.call(rbind, lapply(TRY_SPCIS_list, select_function))

# merge TRY data with SPCIS matched names dataset 
TRY_SPCIS_df <- merge(TRY_SPCIS[,c(1:9)], TRY_SPCIS_list_out[,c(1,2,5)], by = c("SpeciesName", "AccSpeciesName"), all.x = TRUE)

# -------------------------- write SPCIS_TRY dataset   -------------------------

write.csv(TRY_SPCIS_df, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_TRY_22398.csv")
