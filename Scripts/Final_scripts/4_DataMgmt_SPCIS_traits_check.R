# Trait shifts 
# merge trait dataset back with SPCIS database 

library(data.table)

traits <- read.csv("/Users/magdagarbowski/TraitShifts/Generated_Data/SPCIS_traits.csv")
spcis_data <- fread("/Users/magdagarbowski/TraitShifts/Data/FULLDatabase_10272022.csv")

traits$X <- NULL

# ------------------------------- SPCIS data cleanup ---------------------------------------------
# drop nonvascular and lichenous 
spcis_data <- as.data.frame(spcis_data[!spcis_data$Growth.Habit %in% c("Nonvascular", "Lichenous"),])

# drop varieties and subspecies
spcis_data$AcceptedTaxonName <- gsub(" var.*| ssp.*", "", spcis_data$AcceptedTaxonName)

# check mismatches 
unique(traits$sps_try_match)[!unique(traits$sps_try_match) %in% unique(spcis_data$AcceptedTaxonName)]
unique(spcis_data$AcceptedTaxonName)[!unique(spcis_data$AcceptedTaxonName) %in% unique(traits$sps_try_match)]
