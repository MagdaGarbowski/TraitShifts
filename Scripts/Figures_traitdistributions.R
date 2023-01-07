
# Trait shifts 
# plotting! 

library(data.table)
library(psych)
library(ggplot2)
library(gridExtra)

# ----------------------------------- plotting CWM -------------------------------------

CWM <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM.csv", 
                           select = c("Plot", "Year", "TraitNameAbr", "CWM")))
colnames(CWM) <- c("Plot", "Year", "Trait", "CWM")

CWM_plot_trait <- paste(CWM$Plot, CWM$Year, CWM$Trait, sep = "_")

plot_info <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

# drop duplicates - these ought to be dropped sooner 
plot_info  <- plot_info[!duplicated(plot_info),]

# merge the datasets 
CWM_plot_info <- merge(CWM, plot_info, by = c("Plot", "Year", "Trait"), all.x = TRUE)

CWM_plot_info$Inv_level <- ifelse(CWM_plot_info$I_relcov >= 66, "high", 
                                  ifelse((CWM_plot_info$I_relcov >= 33 & CWM_plot_info$I_relcov < 66), "med", "low"))

# drop ecoregions without a lot of plots 
CWM_plot_info <- CWM_plot_info[!CWM_plot_info$EcoRegionLevelI %in% c("TAIGA", "TUNDRA", "WATER", "TROPICAL WET FORESTS"),]

# select out plots with > 80% coverage of traits 
CWM_80 <- CWM_plot_info[CWM_plot_info$trait_cov >= 80,]

# split by trait 
CWM_80_splits <- split(CWM_80, CWM_80$Trait)

density_plot_function <- function(df,trait, x_min, x_max){
  n_plots <- nrow(df)
  plot_out <- ggplot(df, aes(x = CWM, group = Inv_level, fill = Inv_level)) + 
    geom_density(adjust = 1.5, alpha = 0.4) + 
    scale_fill_manual(breaks = c("low", "med", "high"),
                      labels = c("low (<33%)", "med (33-66%)", "high (>66%)"),
                      values = c( "#CCCCCC", "#4f5157","#3c5e27")) + 
    xlim(x_min, x_max) + 
    labs(x = trait, title = paste(trait, "\n", "n_plots =", n_plots)) + 
    theme_bw()
  return(plot_out)
}

height_gen <- density_plot_function(CWM_80_splits$height_gen, "Height-gen (m)", 0, 3)
height_veg <- density_plot_function(CWM_80_splits$height_veg, "Height-veg (m)", 0, 5)
LDMC <- density_plot_function(CWM_80_splits$LDMC, "LDMC (g g-1)", 0, 0.6)
leaf_area <- density_plot_function(CWM_80_splits$leaf_area, "Leaf Area (mm2)", 0, 3000)
leaf_N <- density_plot_function(CWM_80_splits$leafN, "leaf_N (mg/g)", 0, 40)
leaf_P <- density_plot_function(CWM_80_splits$leafP, "leaf_P (mg/g)", 0, 4)
seed_mass <- density_plot_function(CWM_80_splits$seed_mass, "seed_mass (mg)", 0, 10)
SLA <- density_plot_function(CWM_80_splits$SLA, "SLA (mm2 mg-1)", 0, 50)
SSD <- density_plot_function(CWM_80_splits$SSD, "SSD", 0, 1)



pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_1.pdf", width = 10, height = 6)
grid.arrange(height_veg, leaf_area, SSD,seed_mass, ncol = 2)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_2.pdf", width = 10, height = 6)
grid.arrange(SLA, LDMC, leaf_N, leaf_P, ncol = 2)
dev.off()

# ----------------------------------- ecoregion plots  -------------------------------------


density_eco_plot_function <- function(dat,trait, x_min, x_max){
  
  eco_splits <- split(dat, dat$EcoRegionLevelI, drop = TRUE)
  inner_fun <- function(df){
    n_plots <- nrow(df)
    eco_reg <- df$EcoRegionLevelI
    plot_out <- ggplot(df, aes(x = CWM, group = Inv_level, fill = Inv_level)) + 
      geom_density(adjust = 1.5, alpha = 0.4) + 
      scale_fill_manual(breaks = c("low", "med", "high"),
                        labels = c("low (<33%)", "med (33-66%)", "high (>66%)"),
                        values = c( "#CCCCCC", "#4f5157","#3c5e27")) + 
      xlim(x_min, x_max) + 
      labs(x = trait, title = paste(eco_reg, "\n",trait, "\n",  "n_plots =", n_plots)) + 
      theme_bw() + 
      theme(legend.position = "none")
    return(plot_out)
  }
  plots_out <- lapply(eco_splits, inner_fun)
  plots_print <- do.call("grid.arrange", c(plots_out, ncol = 3))
  return(plots_print)
}

height_veg_eco <- density_eco_plot_function(CWM_80_splits$height_veg, "Height_veg (m)", 0, 15)
LDMC_eco <- density_eco_plot_function(CWM_80_splits$LDMC, "LDMC (g g-1)", 0, 0.6)
leaf_area_eco <- density_eco_plot_function(CWM_80_splits$leaf_area, "leaf_area (mm2)", 0, 3000)
leafN_eco <- density_eco_plot_function(CWM_80_splits$leafN, "leaf_N (mg/g)", 0, 40)
leafP_eco <- density_eco_plot_function(CWM_80_splits$leafP, "leaf_P (mg/g)", 0, 4)
seed_mass_eco <- density_eco_plot_function(CWM_80_splits$seed_mass, "seed_mass (mg)", 0, 10)
SLA_eco <- density_eco_plot_function(CWM_80_splits$SLA, "SLA (mm2 mg-1)", 0, 50)
SSD_eco <- density_eco_plot_function(CWM_80_splits$SSD, "SSD", 0, 1)

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/height_veg_eco.pdf", width = 10, height = 10)
density_eco_plot_function(CWM_80_splits$height_veg, "Height_veg", 0, 8)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/seed_mass_eco.pdf", width = 10, height = 10)
density_eco_plot_function(CWM_80_splits$seed_mass, "seed_mass", 0, 10)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/SLA_eco.pdf", width = 10, height = 10)
density_eco_plot_function(CWM_80_splits$SLA, "SLA", 0, 50)
dev.off()




# ----------------------------------- which plots have >80 of several traits  -------------------------------------

CWM_80_plot_splits <- split(CWM_80, CWM_80$Plot, drop = TRUE)

trait_count_function <- function(df){
  yrs = length(levels(as.factor(df$Year)))
  
  inner_function <- function(df){
  out <- data.frame(Plot = df$Plot[1],
                    Year = df$Year[1],
                    Inv_level = df$Inv_level[1],
                    EcoRegionLevelI = df$EcoRegionLevelI[1],
                    trait_levels = length(levels(as.factor(df$Trait))),
                    traits = toString(levels(as.factor(df$Trait))))
  return(out)
  }
  
  if(yrs == 1){
    out <- inner_function(df)
    return(out)
  }
  if(yrs > 1){
    df_yr_splits <- split(df, list(df$Year))
    out_all <- do.call(rbind, lapply(df_yr_splits, inner_function))
    return(out_all)
  }
}

cov_80_traits_df <- do.call(rbind, lapply(CWM_80_plot_splits, trait_count_function))

cov_80_3plus_traits <- cov_80_traits_df[cov_80_traits_df$trait_levels > 2,]

cov_80_traits_df_ecosplits <- split(cov_80_3plus_traits, 
                                    list(cov_80_3plus_traits$EcoRegionLevelI,
                                         cov_80_3plus_traits$Inv_level, drop = TRUE))

bytrait_count_function <- function(df){
  out <- data.frame(EcoRegionLevelI = df$EcoRegionLevelI[1],
                    Inv_level = df$Inv_level[1],
                    LHS = nrow(df[((grepl("height_veg", df$traits)) & (grepl("SLA", df$traits)) & (grepl("seed_mass", df$traits))),]),
                    traits_3 = nrow(df),
                    traits_4 = nrow(df[df$trait_levels > 3,]),
                    traits_5 = nrow(df[df$trait_levels > 4,]),
                    traits_6 = nrow(df[df$trait_levels > 5,]),
                    traits_7 = nrow(df[df$trait_levels > 6,]),
                    traits_8 = nrow(df[df$trait_levels > 7,]))
  return(out)
}

# export this table 
trait_counts_df <- do.call(rbind, lapply(cov_80_traits_df_ecosplits, bytrait_count_function))
           
# work with deserts to try to get a multivariate plot
#
#
#
#
dat <- fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wTRY.csv", 
             select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean"))

# ecoregions PCA 

# deserts 
deserts <- cov_80_3plus_traits[(cov_80_3plus_traits$EcoRegionLevelI == "NORTH AMERICAN DESERTS" & 
                                              cov_80_3plus_traits$trait_levels > 6),]
           
deserts_CWM_80 <- CWM_80[CWM_80$Plot %in% deserts$Plot,]

# eastern temperate forests 
east_temp_forests <- cov_80_3plus_traits[(cov_80_3plus_traits$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS" & 
                                  cov_80_3plus_traits$trait_levels > 6),]

east_temp_forests_CWM_80 <- CWM_80[CWM_80$Plot %in% east_temp_forests$Plot,]

# great plains 
great_plains <- cov_80_3plus_traits[(cov_80_3plus_traits$EcoRegionLevelI == "GREAT PLAINS" & 
                                            cov_80_3plus_traits$trait_levels > 6),]

great_plains_CWM_80 <- CWM_80[CWM_80$Plot %in% great_plains$Plot,]

# data prep 

out_CWM_80 <- lapply(list(deserts_CWM_80 = deserts_CWM_80,
                          east_temp_forests_CWM_80 = east_temp_forests_CWM_80,
                          great_plains_CWM_80 = great_plains_CWM_80),
                       function(x) {xx = x[c("Plot", "Year", "Trait", "CWM", "Inv_level")];
                       return(xx)})

list2env(out_CWM_80, .GlobalEnv)

deserts_CWM_80_splits <- split(deserts_CWM_80, list(deserts_CWM_80$Plot, deserts_CWM_80$Year), drop = TRUE)
east_temp_forests_CWM_80_splits <- split(east_temp_forests_CWM_80, list(east_temp_forests_CWM_80$Plot, east_temp_forests_CWM_80$Year), drop = TRUE)
great_plains_CWM_80_splits <- split(great_plains_CWM_80, list(great_plains_CWM_80$Plot, great_plains_CWM_80$Year), drop = TRUE)


# go from long to wide for each dataset 

all_traits <- levels(as.factor(deserts_CWM_80$Trait))

wide_function <- function(df){
  df_traits <- levels(as.factor(df$Trait))
  missing_traits <- all_traits[!all_traits %in% df_traits]
  row_n = nrow(df)
  
  if(row_n == 9){
    df_wide <- reshape(df, 
                       idvar = c("Plot", "Year", "Inv_level"),
                       timevar = c("Trait"), 
                       direction = "wide")
    return(df_wide)
  }
  
  if(row_n < 9){
    missing_df <- data.frame(Plot = df$Plot[1], 
                             Year = df$Year[1], 
                             Trait = missing_traits, 
                             CWM = NA, 
                             Inv_level = df$Inv_level[1])
    
    df_all <- rbind(df, missing_df)
    df_wide <- reshape(df_all, 
                       idvar = c("Plot", "Year", "Inv_level"),
                       timevar = c("Trait"), 
                       direction = "wide")
    
    return(df_wide)
  }
}

desert_wide <- do.call(rbind, lapply(deserts_CWM_80_splits, wide_function))
east_temp_forests_wide <- do.call(rbind, lapply(east_temp_forests_CWM_80_splits, wide_function))
great_plains_wide <- do.call(rbind, lapply(great_plains_CWM_80_splits, wide_function))


# log transform height, leaf area, seed mass, SLA

desert_wide[c("CWM.log.height_veg", "CWM.log.leaf_area", "CWM.log.seed_mass", "CWM.log.SLA")] <- 
  apply(desert_wide[c("CWM.height_veg", "CWM.leaf_area", "CWM.seed_mass", "CWM.SLA")], 2, function(x) {log(x); return(x)})

trait_counts <- sapply(desert_wide, function(x) sum(is.na(x)))

desert_wide_ss <- desert_wide[c("Plot", "Year", "Inv_level", "CWM.log.height_veg","CWM.log.seed_mass", "CWM.log.SLA",
                                "CWM.LDMC", "CWM.leafN", "CWM.leafP")]

# get complete cases for now... 4512 plots PRETTY GOOD! 
desert_wide_ss_cc <- desert_wide_ss[complete.cases(desert_wide_ss),]

########

# log transform height, leaf area, seed mass, SLA

east_temp_forests_wide[c("CWM.log.height_veg", "CWM.log.leaf_area", "CWM.log.seed_mass", "CWM.log.SLA")] <- 
  apply(east_temp_forests_wide[c("CWM.height_veg", "CWM.leaf_area", "CWM.seed_mass", "CWM.SLA")], 2, function(x) {log(x); return(x)})

trait_counts <- sapply(east_temp_forests_wide, function(x) sum(is.na(x)))

east_temp_forests_wide_ss <- east_temp_forests_wide[c("Plot", "Year", "Inv_level", "CWM.log.height_veg","CWM.log.seed_mass", "CWM.log.SLA",
                                "CWM.LDMC", "CWM.leafN", "CWM.leafP")]

# get complete cases for now... 3032 plots PRETTY GOOD! 
east_temp_forests_wide_ss_cc <- east_temp_forests_wide_ss[complete.cases(east_temp_forests_wide_ss),]

# log transform height, leaf area, seed mass, SLA

great_plains_wide[c("CWM.log.height_veg", "CWM.log.leaf_area", "CWM.log.seed_mass", "CWM.log.SLA")] <- 
  apply(great_plains_wide[c("CWM.height_veg", "CWM.leaf_area", "CWM.seed_mass", "CWM.SLA")], 2, function(x) {log(x); return(x)})

trait_counts <- sapply(great_plains_wide, function(x) sum(is.na(x)))

great_plains_wide_ss <- great_plains_wide[c("Plot", "Year", "Inv_level", "CWM.log.height_veg","CWM.log.seed_mass", "CWM.log.SLA",
                                "CWM.LDMC", "CWM.leafN", "CWM.leafP")]

# get complete cases for now... 899 plots PRETTY GOOD! 
great_plains_wide_ss_cc <- great_plains_wide_ss[complete.cases(great_plains_wide_ss),]


# ---------------------------------------- PCA from psych to prcomp obj for plotting ----------------------------

pca.obj <- function(df, principal_obj, pc_names){
  princ.load <- as.matrix.data.frame(principal_obj$loadings, rownames = TRUE)
  colnames(princ.load) <- pc_names
  princ.scores<-as.data.frame(as.matrix(principal_obj$scores, rownames = TRUE))
  colnames(princ.scores) <- pc_names
  princ.scores<-as.matrix(princ.scores)
  pca_var.obj <- list(sdev = NULL, rotation = princ.load, 
                      center=NULL, scale = NULL, x = princ.scores) 
  class(pca_var.obj) <- "prcomp"
  return(pca_var.obj)
}


plot_function <- function(df, ecoregion, x_lab, y_lab){
  plot_out <- ggplot(df, aes(x = PC1, y = PC2, color = Inv_level)) +
    geom_point(alpha = 0.4) +
    stat_ellipse() +      
    scale_color_manual(breaks = c("low", "med", "high"),
                       labels = c("low (<33%)", "med (33-66%)", "high (>66%)"),
                       values = c( "#CCCCCC", "#4f5157","#3c5e27")) + 
    labs (x = x_lab,
          y = y_lab,
          title = ecoregion) + 
    theme_bw() + 
    theme(legend.position = "none")+
    ylim(-2.5, 7)
  return(plot_out)
}

# Need to edit labels and PC axes for each set of plots 


deserts_pca <- principal(desert_wide_ss_cc[,-c(1:3)], 3, scores = TRUE)
deserts_pca_obj <- pca.obj(desert_wide_ss_cc[,-c(1:3)], deserts_pca, c("PC1","PC2","PC3"))
deserts_pca_df <- as.data.frame(deserts_pca_obj$x)
deserts_pca_df$Inv_level <- desert_wide_ss_cc$Inv_level
deserts <- plot_function(deserts_pca_df, "deserts \n n_plots = 4512",
                         "PC1, 33% var. explained \n leaf_N, leaf_P, -LDMC",
                         "PC2, 31% var. explained \n log_height, log_seedmass")

great_plains_pca <- principal(great_plains_wide_ss_cc[,-c(1:3)], 3, scores = TRUE)
great_plains_pca_obj <- pca.obj(great_plains_wide_ss_cc[,-c(1:3)], great_plains_pca, c("PC1","PC2","PC3"))
great_plains_pca_df <- as.data.frame(great_plains_pca_obj$x)
great_plains_pca_df$Inv_level <- great_plains_wide_ss_cc$Inv_level
gp <- plot_function(great_plains_pca_df, "great_plains \n n_plots = 899",
                    "PC1, 34% var. explained \n log_SLA, -LDMC, leaf_P",
                    "PC2, 26% var. explained \n log_height, log_seedmass")

east_temp_forests_pca <- principal(east_temp_forests_wide_ss_cc[,-c(1:3)], 3, scores = TRUE)
east_temp_forests_pca_obj <- pca.obj(east_temp_forests_wide_ss_cc[,-c(1:3)], east_temp_forests_pca, c("PC1","PC2","PC3"))
east_temp_forests_pca_df <- as.data.frame(east_temp_forests_pca_obj$x)
east_temp_forests_pca_df$Inv_level <- east_temp_forests_wide_ss_cc$Inv_level
e.t.forests <- plot_function(east_temp_forests_pca_df, "east_temp_forests \n n_plots = 3032",
                             "PC1, 31% var. explained \n log_height, log_seedmass",
                             "PC2, 27% var. explained \n leaf_N, leaf_P")

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/des_gp_etf_PCA.pdf", width = 10, height = 4)
grid.arrange(deserts, gp, e.t.forests, ncol = 3)
dev.off()


