# Trait-Shifts 
# goal: create native invasive plots

library(data.table)
library(ggplot2)
library(gridExtra)
library(usmap)
library(maptools)
library(rgdal)
library(tidyr)

traits <- read.csv("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_traits_sumstats.csv")
spcis_data <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Data/FULLDatabase_10272022.csv"))

# drop nonvascular and lichenous 
spcis_data <- as.data.frame(spcis_data[!spcis_data$Growth.Habit %in% c("Nonvascular", "Lichenous"),])

# drop varieties and subspecies
spcis_data$AcceptedTaxonName <- gsub(" var.*| ssp.*", "", spcis_data$AcceptedTaxonName)
spcis_data <- spcis_data[!is.na(spcis_data$Plot),]

traits$X <- NULL

# get invasion status for unique species in SPCIS 

spcis_unique <- unique(spcis_data[c("AcceptedTaxonName", "NativeStatus")])

spcis_unique_ecoregion <- unique(spcis_data[c("AcceptedTaxonName", "NativeStatus", "EcoRegionLevelI")])

# -------------------- >5000 with missing traits  ----------------------------------

sps_out <- aggregate(list(n_plot_occurance = spcis_data$PctCov), by = list(AcceptedTaxonName = spcis_data$AcceptedTaxonName), FUN = length)
sps_out_5000 <- sps_out[sps_out$n_plot_occurance > 5000,]

# get all possible species by trait combinations 
traits_all <- c("heightveg_m", "LDMC_g/g", "leafN_mg/g", "leafP_mg/g", "seedmass_mg", "SLA_mm2/mg", "SSD_g/cm3",
                "Mean_Root_diameter","Root_N_concentration", "Specific_root_length", "max_rooting_depth_m")

sps_5000_traits <- expand.grid(species = sps_out_5000$AcceptedTaxonName, TraitNameAbr = traits_all)

sps_5000_w_traits <- merge(sps_5000_traits, traits[c("sps_try_match", "TraitNameAbr", "mean", "n_studies", "n_species")], 
                              by.x = c("species", "TraitNameAbr"), by.y = c("sps_try_match", "TraitNameAbr"), all.x = TRUE)

sps_5000_w_traits_missing <- sps_5000_w_traits[is.na(sps_5000_w_traits$n_studies),]

write.csv(sps_5000_w_traits_missing, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/abundant_sps_totallymissing_data.csv")


# -------------------- traits by species occurrence  ----------------------------------
# species by occurrence counts 
sps_out <- aggregate(list(n_plot_occurance = spcis_data$PctCov), by = list(AcceptedTaxonName = spcis_data$AcceptedTaxonName), FUN = length)
sps_occurance_traits <- merge(sps_out, traits[c("sps_try_match", "TraitNameAbr", "mean", "n_studies", "n_species")], 
                              by.x = "AcceptedTaxonName", by.y = "sps_try_match", all.x = TRUE)

# drop categorical traits 
sps_occurance_traits <- sps_occurance_traits[!sps_occurance_traits$TraitNameAbr %in% c("Duration", "Growth.Habit","Mycorrhizal.type"),]
sps_occurance_traits$mean <- as.numeric(sps_occurance_traits$mean)

# subset traits of interest
sps_occurance_traits <- sps_occurance_traits[sps_occurance_traits$TraitNameAbr %in% c("heightveg_m","LDMC_g/g","leafN_mg/g","leafP_mg/g","Mean_Root_diameter","Root_tissue_density",
                                               "max_rooting_depth_m","seedmass_mg","SLA_mm2/mg","Specific_root_length","SSD_g/cm3"),]

sps_occurance_traits_splits <- split(sps_occurance_traits, sps_occurance_traits$TraitNameAbr, drop = TRUE)

plot_occurance_function <- function(df, n_obs_val, y_min, y_max, x_min, x_max){
  plot_out <- ggplot(df, aes(x = n_plot_occurance, y = n_studies)) + 
    geom_point(alpha = 0.5) +
    geom_text(data = df[(df$n_studies > n_obs_val | df$n_plot_occurance > 5000),],
              aes(x = n_plot_occurance, y = n_studies, label = AcceptedTaxonName),
              size = 3, angle = 45, hjust = -0.1) + 
    labs(title = paste(df$TraitNameAbr[[1]])) + 
    ylim(y_min, y_max) + 
    xlim(x_min, x_max)
  return(plot_out)
}

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/occurances.pdf", width = 16, height = 10)
plot_occurance_function(sps_occurance_traits_splits$heightveg_m, 200, -10, 400, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$`LDMC_g/g`, 200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$`leafN_mg/g`, 200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$`leafP_mg/g`,200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$max_rooting_depth_m, 200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$Mean_Root_diameter, 200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$Root_tissue_density, 200, -10, 200, 0, 19000)
plot_occurance_function(sps_occurance_traits_splits$seedmass_mg, 200, -10, 200, 0, 19000) 
plot_occurance_function(sps_occurance_traits_splits$`SLA_mm2/mg`, 200, -10, 200, 0, 19000) 
plot_occurance_function(sps_occurance_traits_splits$Specific_root_length, 200, -10, 200, 0, 19000) 
plot_occurance_function(sps_occurance_traits_splits$`SSD_g/cm3`, 200, -10, 200, 0, 19000) 
dev.off()


# create tables of species with >7500 plot occurrences and < 5% of quantile range for n_obs

occurance_tb_function <- function(df) {
  df_2 <- df[(df$n_plot_occurance > 5000 & df$n_studies < 3),]
  return(df_2)
}

occ_tab_out <- do.call(rbind, lapply(sps_occurance_traits_splits, occurance_tb_function))
write.csv(occ_tab_out, "/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/abundant_sps_missing_data.csv")

# ------------------------------- traits ---------------------------------------------
traits_ss <- traits[traits$TraitNameAbr %in% c("heightveg_m","LDMC_g/g","leafN_mg/g","leafP_mg/g","Mean_Root_diameter","Root_tissue_density",
                                               "max_rooting_depth_m","seedmass_mg","SLA_mm2/mg","Specific_root_length","SSD_g/cm3"),]

traits_ss <- traits_ss[,c(1:3)]

# ------------------------------- traits table  ---------------------------------------------
spcis_traits <- merge(spcis_unique, traits, by.x = "AcceptedTaxonName", by.y = "sps_try_match", all.y = TRUE)
spcis_traits <- spcis_traits[!is.na(spcis_traits$NativeStatus),]

out <- aggregate(mean ~ NativeStatus + TraitNameAbr, data = spcis_traits, FUN = length)
out_wide <- reshape(out, idvar = "TraitNameAbr", timevar = "NativeStatus", direction = "wide" )
out_wide$all <- rowSums(out_wide[,2:4])
colnames(out_wide) <- c("Trait", "I", "N", "NI", "all")
out_wide <- out_wide[order(-out_wide$all),]

out_wide_top20 <- out_wide[1:20,]

write.csv(out_wide,  "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/SPCIS_trait_numbers.csv")

# ------------------------------- overall spcis dataset ---------------------------------------------

spcis_traits <- merge(spcis_unique, traits_ss, by.x = "AcceptedTaxonName", by.y = "sps_try_match", all.y = TRUE)
spcis_traits$NativeStatus <- ifelse(spcis_traits$NativeStatus == "NI", "N", spcis_traits$NativeStatus)
spcis_traits$mean <- as.numeric(spcis_traits$mean)

spcis_traits_splits <- split(spcis_traits, spcis_traits$TraitNameAbr)

plot_function <- function(df, x_min, x_max){
  N_number <- nrow(df[df$NativeStatus == "N",])
  I_number <- nrow(df[df$NativeStatus == "I",])
  trait = df$TraitNameAbr[[1]]
  
  plot_out <- ggplot(df, aes(x = mean, group = NativeStatus, fill = NativeStatus)) + 
    geom_density(adjust = 1.5, alpha = 0.5) + 
    scale_fill_manual(values = c("#E69F00","#999999")) + 
    xlim(x_min, x_max) + 
    labs(x = trait, y = "density",
         title = paste(trait, "\n","natives = ",N_number, "\n","invasives = ", I_number)) + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  ggplot_data <- ggplot_build(plot_out + stat_density())$data[[2]]
  
  Btectorum_x_value <- mean(df[df$AcceptedTaxonName == "Bromus tectorum",]$mean)
  Atridentata_x_value <- mean(df[df$AcceptedTaxonName == "Artemisia tridentata",]$mean)
  
  plot_out_2 <- plot_out + 
    annotate("text", x = Btectorum_x_value, y = 0, label = "*BRTE", hjust = 0, size = 4, angle = 60) + 
    annotate("text", x = Atridentata_x_value, y = 0, label = "*ARTR", hjust = 0, size = 4, angle = 60)
  
  return(plot_out_2)
} 

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/native_vs_invasive.pdf", width = 16, height = 10)

grid.arrange(plot_function(spcis_traits_splits$heightveg_m, 0, 5),
             plot_function(spcis_traits_splits$`SLA_mm2/mg`, 0, 85),
             plot_function(spcis_traits_splits$`leafN_mg/g`, 0, 46),
             plot_function(spcis_traits_splits$max_rooting_depth_m, 0, 12), 
             plot_function(spcis_traits_splits$Specific_root_length, 0, 500),
             plot_function(spcis_traits_splits$Root_tissue_density, 0, 1),
              ncol = 3)

dev.off()

# ------------------------------- EcoRegion ---------------------------------------------

spcis_traits <- merge(spcis_unique, traits_ss, by.x = "AcceptedTaxonName", by.y = "sps_try_match", all.y = TRUE)

spcis_traits_eco <- merge(spcis_traits, spcis_unique_ecoregion, by = c("AcceptedTaxonName", "NativeStatus"), all.x = TRUE)

spcis_traits_eco <- spcis_traits_eco[!is.na(spcis_traits_eco$NativeStatus),]
spcis_traits_eco$NativeStatus <- ifelse(spcis_traits_eco$NativeStatus == "NI", "N", spcis_traits_eco$NativeStatus)
spcis_traits_eco$mean <- as.numeric(spcis_traits_eco$mean)
spcis_traits_eco <- spcis_traits_eco[!is.na(spcis_traits_eco$EcoRegionLevelI),]

spcis_traits_eco_ss <- spcis_traits_eco[spcis_traits_eco$EcoRegionLevelI %in% c("EASTERN TEMPERATE FORESTS", "GREAT PLAINS", "NORTH AMERICAN DESERTS", "NORTHERN FORESTS"),]

spcis_traits_splits <- split(spcis_traits_eco_ss, list(spcis_traits_eco_ss$TraitNameAbr), drop = TRUE)

plot_function_2 <- function(df, x_min, x_max){
  eco_splits <- split(df, df$EcoRegionLevelI)
  
  plot_function_inner <- function(dat){
  trait = dat$TraitNameAbr[[1]]
  ecoregion = dat$EcoRegionLevelI[[1]]
  
  plot_out <- ggplot(dat, aes(x = mean, group = NativeStatus, fill = NativeStatus)) + 
    geom_density(adjust = 1.5, alpha = 0.5) + 
    scale_fill_manual(values = c("#E69F00","#999999")) + 
    xlim(x_min, x_max) + 
    labs(x = trait, title = paste(trait, ecoregion)) + 
    theme_bw() + 
    theme(legend.position = "bottom")
  
  return(plot_out)
  }
  out_plot <- do.call("grid.arrange", lapply(eco_splits, plot_function_inner))
}

plot_function_2(spcis_traits_splits$heightveg_m, 0, 8)
plot_function_2(spcis_traits_splits$`SLA_mm2/mg`, 0, 80)
plot_function_2(spcis_traits_splits$`leafN_mg/g`, 0, 52)
plot_function_2(spcis_traits_splits$Root_tissue_density, 0, 1)
plot_function_2(spcis_traits_splits$max_rooting_depth_m, 0, 7)
plot_function_2(spcis_traits_splits$Mean_Root_diameter, 0, 3)

# ------------------------------- Deserts and Grasslands EcoRegion ---------------------------------------------

spcis_data_DG <- spcis_data[spcis_data$EcoRegionLevelI %in% c("NORTH AMERICAN DESERTS", "GREAT PLAINS", "MEDITERRANEAN CALIFORNIA", "NORTHWESTERN FORESTED MOUNTAINS"),]
spcis_data_DG <- spcis_data_DG[spcis_data_DG$Dataset == "AIM",]

spcis_data_DG_lat_long <- unique(spcis_data_DG[c("Plot", "Lat", "Long", "EcoRegionLevelI")])

spcis_data_DG <- spcis_data_DG[grep(" ", spcis_data_DG$AcceptedTaxonName),]

sps_out_DG <- aggregate(list(n_plot_occurance = spcis_data_DG$PctCov), by = list(AcceptedTaxonName = spcis_data_DG$AcceptedTaxonName, EcoRegionLevelI = spcis_data_DG$EcoRegionLevelI), FUN = length)
sps_out_DG <- sps_out_DG[sps_out_DG$EcoRegionLevelI %in% c("NORTH AMERICAN DESERTS", "GREAT PLAINS", "NORTHWESTERN FORESTED MOUNTAINS"),]


sps_out_DG_agg <- aggregate(sps_out_DG$n_plot_occurance, by = list(sps_out_DG$AcceptedTaxonName), FUN = "sum")
sps_out_DG_agg_50 <- sps_out_DG_agg[sps_out_DG_agg$x > 100,]
colnames(sps_out_DG_agg_50) <- c("AcceptedTaxonName", "n_plot_occurance")
sps_out_DG_agg_50$EcoRegionLevelI <- "GP_NAD_NFM"

sps_out_DG_50 <- sps_out_DG[sps_out_DG$n_plot_occurance > 100,]

sps_out_D_50 <- sps_out_DG_50[sps_out_DG_50$EcoRegionLevelI == "NORTH AMERICAN DESERTS",]
sps_out_G_50 <- sps_out_DG_50[sps_out_DG_50$EcoRegionLevelI == "GREAT PLAINS",]
sps_out_F_50 <- sps_out_DG_50[sps_out_DG_50$EcoRegionLevelI == "NORTHWESTERN FORESTED MOUNTAINS",]

# get all possible species by trait combinations 
traits_all <- c("heightveg_m", "LDMC_g/g", "leafN_mg/g", "leafP_mg/g", "seedmass_mg", "SLA_mm2/mg", "SSD_g/cm3",
                "Mean_Root_diameter","Root_N_concentration", "Specific_root_length", "max_rooting_depth_m")

### DESERTS ###
sps_D_traits <- expand.grid(species = sps_out_D_50$AcceptedTaxonName, TraitNameAbr = traits_all)

sps_D_traits_all <- merge(sps_D_traits, traits[c("sps_try_match", "TraitNameAbr", "mean")], 
                           by.x = c("species", "TraitNameAbr"), by.y = c("sps_try_match", "TraitNameAbr"), all.x = TRUE)

sps_D_traits_all <- sps_D_traits_all[!is.na(sps_D_traits_all$species),]

# split by trait 
sps_D_traits_splits <- split(sps_D_traits_all, sps_D_traits_all$TraitNameAbr)

### GREAT PLAINS ###

sps_G_traits <- expand.grid(species = sps_out_G_50$AcceptedTaxonName, TraitNameAbr = traits_all)

sps_G_traits_all <- merge(sps_G_traits, traits[c("sps_try_match", "TraitNameAbr", "mean")], 
                          by.x = c("species", "TraitNameAbr"), by.y = c("sps_try_match", "TraitNameAbr"), all.x = TRUE)

sps_G_traits_all <- sps_G_traits_all[!is.na(sps_G_traits_all$species),]

# split by trait 
sps_G_traits_splits <- split(sps_G_traits_all, sps_G_traits_all$TraitNameAbr)


### FORESTED MOUNTAINS ###

sps_F_traits <- expand.grid(species = sps_out_F_50$AcceptedTaxonName, TraitNameAbr = traits_all)

sps_F_traits_all <- merge(sps_F_traits, traits[c("sps_try_match", "TraitNameAbr", "mean")], 
                          by.x = c("species", "TraitNameAbr"), by.y = c("sps_try_match", "TraitNameAbr"), all.x = TRUE)

sps_F_traits_all <- sps_F_traits_all[!is.na(sps_F_traits_all$species),]

# split by trait 
sps_F_traits_splits <- split(sps_F_traits_all, sps_F_traits_all$TraitNameAbr)

### ALL THREE ECOREGIONS TOGETHER ### 

sps_DFG_traits <- expand.grid(species = sps_out_DG_agg_50$AcceptedTaxonName, TraitNameAbr = traits_all)

sps_DFG_traits <- unique(sps_DFG_traits)

sps_DFG_traits_all <- merge(sps_DFG_traits, traits[c("sps_try_match", "TraitNameAbr", "mean")], 
                          by.x = c("species", "TraitNameAbr"), by.y = c("sps_try_match", "TraitNameAbr"), all.x = TRUE)

sps_DFG_traits_all <- sps_DFG_traits_all[!is.na(sps_DFG_traits_all$species),]

# split by trait 
sps_DFG_traits_splits <- split(sps_DFG_traits_all, sps_DFG_traits_all$TraitNameAbr)

summary_function <- function(df, df_2){
  
  total_sps = nrow(df_2)
  
  n_sps_missing <- nrow(df[is.na(df$mean),])
  n_sps_present <- nrow(df[!is.na(df$mean),])
  trait <- df$TraitNameAbr[1]
  
  out <- data.frame(trait = trait, 
                    Missing_Species = n_sps_missing, 
                    Existing_Coverage = round(n_sps_present/total_sps, 2) * 100,
                    EcoRegionLevel1 = df_2$EcoRegionLevelI[1])
  return(out)
  
}

deserts_out <- do.call(rbind, lapply(sps_D_traits_splits, summary_function, sps_out_D_50))
great_plains_out <- do.call(rbind, lapply(sps_G_traits_splits, summary_function, sps_out_G_50))
forests_out <- do.call(rbind, lapply(sps_F_traits_splits, summary_function, sps_out_F_50))
all_out <- do.call(rbind, lapply(sps_DFG_traits_splits, summary_function, sps_out_DG_agg_50))
  
out_DFG <- do.call(rbind, list(deserts_out, great_plains_out,forests_out))

write.csv(all_out,  "/Users/MagdaGarbowski 1/Desktop/DFG_100_traits.csv")

# ------------------------------ world map plot --------------------------------------------#
spcis_data_AIM <- spcis_data[spcis_data$Dataset == "AIM",]
state <- map_data("state")
states <- state[state$region %in% c("colorado", "wyoming", "utah", "arizona", "nevada", "oregon", "washington",
                                    "california", "idaho","new mexico", "montana", "south dakota", "north dakota"),]
                                    
AIM_plot <- ggplot(data=states, mapping=aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill="white") + 
  geom_polygon(color="black", fill=NA) + 
  ggtitle("AIM Dataset")  + 
  theme_bw() + 
  theme(legend.text = element_text(size = 18), 
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        legend.title = element_text(size = 18))+
  geom_point(aes(x = Long, y = Lat, group = EcoRegionLevelI, color = EcoRegionLevelI),
             size = 1,
             shape = 16,
             alpha = 0.5,
             data = spcis_data_AIM) + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
  labs(x = "Latitude", y = "Longitude")

pdf(file = "/Users/MagdaGarbowski 1/Desktop/AIM_plot.pdf", width = 16, height = 12)
AIM_plot
dev.off()


spcis_data_AIM