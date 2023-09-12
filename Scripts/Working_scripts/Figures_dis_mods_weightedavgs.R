# Trait shifts 
# plotting! 

library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(cowplot)
library(MuMIn)
library(gridExtra)
library(merTools)

# ----------------------------------- load data -------------------------------------

CWM <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM.csv", 
                           select = c("Plot", "Year", "TraitNameAbr", "CWM")))

dat <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wtraits_2.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean")))

plot_info <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage_2.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

SPCIS_ecoregion_IV <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Data/SPCIS_ecoreg_03182023.csv", 
                                          select = c("Plot", "US_L4CODE", "US_L4NAME")))

colnames(CWM) <- c("Plot", "Year", "Trait", "CWM")

cb_palette <- c("#006633",  "#E69F00","#F0E442", "#D55E00", "#009E73", "#56B4E9", "black")

# ----------------------------------- data prep -------------------------------------

# drop duplicates - these ought to be dropped sooner 
plot_info  <- plot_info[!duplicated(plot_info),]

# keep only last year of plots that have multiple samplings
plot_info_splits <- split(plot_info, plot_info$Plot)

max_year_function <- function(df){
  max_year = max(df$Year)
  df_out <- df[df$Year == max_year,]
  return(df_out)
}

plot_info_max_yrs <- do.call(rbind, lapply(plot_info_splits, max_year_function))

# merge the datasets and subset plots with >80% coverage by trait 
CWM_plot_info <- merge(CWM, plot_info_max_yrs, by = c("Plot", "Year", "Trait"))
CWM_plot_info_80 <- CWM_plot_info[CWM_plot_info$trait_cov > 80,]
CWM_plot_info_80 <-  CWM_plot_info_80[CWM_plot_info_80$EcoRegionLevelI %in% c("EASTERN TEMPERATE FORESTS", "GREAT PLAINS", "NORTHERN FORESTS",
                                                              "MEDITERRANEAN CALIFORNIA", "NORTH AMERICAN DESERTS", "NORTHWESTERN FORESTED MOUNTAINS"),]

# ------ three datasets: full, >2% invasion plots, plots that have traits for invaders and native ------
# all analyses will be run using these datasets 

# full
CWM_plot_info_80 <- CWM_plot_info_80[!CWM_plot_info_80$CWM == -Inf,]

# 2% > invasion 
CWM_plot_info_80_2percent <- CWM_plot_info_80[CWM_plot_info_80$I_relcov > 1.999,]

# create "bins" based on plot numbers, not absolute values of invasion
# this didn't go well because so many plots have zero invasion 
# decided just to "eyeball" it 

CWM_plot_info_80$Inv_level <- ifelse(CWM_plot_info_80$I_relcov < 33, "low", 
                                  ifelse(CWM_plot_info_80$I_relcov >= 33 & CWM_plot_info_80$I_relcov < 66, "med", "high"))

CWM_plot_info_80_2percent$Inv_level <- ifelse(CWM_plot_info_80_2percent$I_relcov < 33, "low", 
                                           ifelse(CWM_plot_info_80_2percent$I_relcov >= 33 & CWM_plot_info_80_2percent$I_relcov < 66, "med", "high"))

# merge with eoregion level 4 

CWM_plot_info_80_2 <- merge(CWM_plot_info_80, SPCIS_ecoregion_IV, by = "Plot", all.x = TRUE)

# ----------------------------------- density plots  -------------------------------------

# split by trait 
CWM_80_splits <- split(CWM_plot_info_80_2, CWM_plot_info_80_2$Trait)

density_plot_function <- function(df,trait, x_min, x_max, level_labs){
  n_plots <- nrow(df)
  plot_out <- ggplot(df, aes(x = CWM, group = Inv_level, fill = Inv_level)) + 
    geom_density(adjust = 1.5, alpha = 0.5) + 
    scale_fill_manual(breaks = c("low", "med", "high"),
                      labels = level_labs,
                      values = c("#999999", "#56B4E9", "#E69F00")) + 
    xlim(x_min, x_max) + 
    labs(x = trait, title = paste(trait, "\n", "n_plots =", n_plots)) + 
    theme_bw() + 
    theme(legend.position = "none") 
  return(plot_out)
}

# size gradient
height_veg_densityplot <- density_plot_function(CWM_80_splits$heightveg_m, "Height-veg (m)", 0, 6, c("low (<33%)", "med (33-66%)", "high (>66%)"))
maxheight_veg_densityplot <- density_plot_function(CWM_80_splits$max_heightveg_m, "Max-Height-veg (m)", 0, 6, c("low (<33%)", "med (33-66%)", "high (>66%)"))

root_depth_densityplot <- density_plot_function(CWM_80_splits$max_rooting_depth_m, "max_root_depth (m)", 0, 12, c("low (<33%)", "med (33-66%)", "high (>66%)"))
root_depth_avg_densityplot <- density_plot_function(CWM_80_splits$mean_rooting_depth_m, "mean_root_depth (m)", 0, 6, c("low (<33%)", "med (33-66%)", "high (>66%)"))

# leaf conservation 
SLA_densityplot <- density_plot_function(CWM_80_splits$`SLA_mm2/mg`, "SLA (mm2 mg-1)", 0, 50, c("low (<33%)", "med (33-66%)", "high (>66%)"))
LDMC_densityplot <- density_plot_function(CWM_80_splits$`LDMC_g/g`, "LDMC (g g-1)", 0, 0.6, c("low (<33%)", "med (33-66%)", "high (>66%)"))
leaf_N_densityplot <- density_plot_function(CWM_80_splits$`leafN_mg/g`, "leaf_N (mg/g)", 0, 40, c("low (<33%)", "med (33-66%)", "high (>66%)"))
leaf_P_densityplot <- density_plot_function(CWM_80_splits$`leafP_mg/g`, "leaf_P (mg/g)", 0, 4, c("low (<33%)", "med (33-66%)", "high (>66%)"))

# roots - conservation and collaboration 
root_diam_densityplot <- density_plot_function(CWM_80_splits$Mean_Root_diameter, "root_diam (mm)", 0, 1, c("low (<33%)", "med (33-66%)", "high (>66%)"))
root_N_densityplot <- density_plot_function(CWM_80_splits$Root_N_concentration, "root_N", 0, 20, c("low (<33%)", "med (33-66%)", "high (>66%)"))
RTD_densityplot <- density_plot_function(CWM_80_splits$Root_tissue_density, "RTD", 0, 0.75, c("low (<33%)", "med (33-66%)", "high (>66%)"))
SRL_densityplot <- density_plot_function(CWM_80_splits$Specific_root_length, "SRL", 0, 350, c("low (<33%)", "med (33-66%)", "high (>66%)"))

# other 
seed_mass_densityplot <- density_plot_function(CWM_80_splits$seedmass_mg, "seed_mass (mg)", 0, 10, c("low (<33%)", "med (33-66%)", "high (>66%)"))
SSD_densityplot <- density_plot_function(CWM_80_splits$`SSD_g/cm3`, "SSD", 0, 1, c("low (<33%)", "med (33-66%)", "high (>66%)"))
RDMC_densityplot <- density_plot_function(CWM_80_splits$Root_dry_matter_content, "RDMC", 0, 0.45, c("low (<33%)", "med (33-66%)", "high (>66%)"))
RMF_densityplot <- density_plot_function(CWM_80_splits$Root_mass_fraction, "RMF", 0, 0.75, c("low (<33%)", "med (33-66%)", "high (>66%)"))
woodiness_densityplot <- density_plot_function(CWM_80_splits$Woodiness, "Woodiness", 0, 1, c("low (<33%)", "med (33-66%)", "high (>66%)"))
annual_densityplot <- density_plot_function(CWM_80_splits$Annual, "Annual_prop", 0, 1, c("low (<33%)", "med (33-66%)", "high (>66%)"))
leaf_area_densityplot <- density_plot_function(CWM_80_splits$leafarea_mm2, "Leaf Area (mm2)", 0, 6600, c("low (<33%)", "med (33-66%)", "high (>66%)"))
root_P_densityplot <- density_plot_function(CWM_80_splits$Root_P_concentration, "root_P", 0, 4, c("low (<33%)", "med (33-66%)", "high (>66%)"))

density_legend <- get_legend(height_veg_densityplot)

# --------------------------------- models - continuous responses -------------------------------------

mod_function <-function(dat){
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  trait <- dat$Trait[1]
  
  if(trait == "SLA_mm2/mg"| trait == "heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_rooting_depth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    mod <- lmer(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI), data = dat)
  }
  
  if(trait != "SLA_mm2/mg"& trait != "heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_rooting_depth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    mod <- lmer(CWM ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI), data = dat)
  }

  summary <- summary(mod)
  pvalue <- summary$coefficients[[10]]
  R2_marginal <- round(r.squaredGLMM(mod)[1],3)
  R2_conditional <- round(r.squaredGLMM(mod)[2],3)
  
  ls <- list(dat, R2_marginal, R2_conditional, mod, pvalue)
  return(ls)
}


height_mod <- mod_function(CWM_80_splits$heightveg_m)
height_maxmod <- mod_function(CWM_80_splits$max_heightveg_m)
root_depth_mod <- mod_function(CWM_80_splits$max_rooting_depth_m)
root_avg_depth_mod <- mod_function(CWM_80_splits$mean_rooting_depth_m)

SLA_mod <- mod_function(CWM_80_splits$`SLA_mm2/mg`)
LDMC_mod <- mod_function(CWM_80_splits$`LDMC_g/g`)
leaf_N_mod <- mod_function(CWM_80_splits$`leafN_mg/g`) # boundary fit singular 
leaf_P_mod <- mod_function(CWM_80_splits$`leafP_mg/g`)

root_diam_mod <- mod_function(CWM_80_splits$Mean_Root_diameter)
RTD_mod <- mod_function(CWM_80_splits$Root_tissue_density)
root_N_mod <- mod_function(CWM_80_splits$Root_N_concentration)
root_P_mod <- mod_function(CWM_80_splits$Root_P_concentration) # boundary fit singular 
SRL_mod <- mod_function(CWM_80_splits$Specific_root_length)

# --------------------------------- updated models (level 4 ecoregions) - continuous responses -------------------------------------

mod_function_updated <-function(dat){
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  trait <- dat$Trait[1]
  
  if(trait == "SLA_mm2/mg"| trait == "heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_rooting_depth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    mod_2 <- lmer(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
  }
  
  if(trait != "SLA_mm2/mg"& trait != "heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_rooting_depth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    mod_2 <- lmer(CWM ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
  }
  
  summary <- summary(mod_2)
  pvalue <- summary$coefficients[[10]]
  R2_marginal <- round(r.squaredGLMM(mod_2)[1],3)
  R2_conditional <- round(r.squaredGLMM(mod_2)[2],3)
  
  ls <- list(dat, R2_marginal, R2_conditional, mod_2, pvalue)
  return(ls)
}

height_mod <- mod_function_updated(CWM_80_splits$heightveg_m)
maxheight_mod <- mod_function_updated(CWM_80_splits$max_heightveg_m) # failed to converge 

root_depth_mod <- mod_function_updated(CWM_80_splits$max_rooting_depth_m)
root_avg_depth_mod <- mod_function_updated(CWM_80_splits$mean_rooting_depth_m)

SLA_mod <- mod_function_updated(CWM_80_splits$`SLA_mm2/mg`)
LDMC_mod <- mod_function_updated(CWM_80_splits$`LDMC_g/g`)
leaf_N_mod <- mod_function_updated(CWM_80_splits$`leafN_mg/g`) 
leaf_P_mod <- mod_function_updated(CWM_80_splits$`leafP_mg/g`)

root_diam_mod <- mod_function_updated(CWM_80_splits$Mean_Root_diameter)
RTD_mod <- mod_function_updated(CWM_80_splits$Root_tissue_density)
root_N_mod <- mod_function_updated(CWM_80_splits$Root_N_concentration)
root_P_mod <- mod_function_updated(CWM_80_splits$Root_P_concentration) # boundary fit singular 
SRL_mod <- mod_function_updated(CWM_80_splits$Specific_root_length)


# ----------------------------------- model plots  -------------------------------------

mod_plot_function <- function(ls){
  # data for plot 
  df = ls[[1]]
  R2_marginal <- ls[[2]]
  R2_conditional <- ls[[3]]
  mod = ls[[4]]
  trait = df$Trait[[1]]
  
  xmin = min(df$I_relcov_scaled)
  xmax = max(df$I_relcov_scaled)

  # new data for predictions
  new_dat_full <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                         EcoRegionLevelI = c(levels(as.factor(df$EcoRegionLevelI))))

  predict_int_full <- predictInterval(mod, new_dat_full, n.sims = 2000,
                                   stat = "median",
                                   type = "linear.prediction",
                                   include.resid.var = TRUE)
  
  new_dat_main <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.0005), 
                              EcoRegionLevelI = c(levels(as.factor(df$EcoRegionLevelI))))
  
  predict_int_main <- predictInterval(mod, new_dat_main, n.sims = 2000,
                                      stat = "median",
                                      which = "fixed",
                                      type = "linear.prediction",
                                      include.resid.var = TRUE)
  
  # data from model 
  if(trait == "SLA_mm2/mg"| trait == "heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_rooting_depth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    predicted_dat_full_exp <- exp(predict_int_full)
    predicted_dat_main_exp <- exp(predict_int_main)
    out_full <- cbind(new_dat_full, predicted_dat_full_exp)
    out_main <- cbind(new_dat_main, predicted_dat_main_exp)
  }
  
  if(trait != "SLA_mm2/mg"& trait != "heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_rooting_depth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    out_full <- cbind(new_dat_full, predict_int_full)
    out_main <- cbind(new_dat_main, predict_int_main)
  }
  
  plot <- ggplot() +
    geom_point(data = df, aes(x = I_relcov_scaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), alpha = 0.2) +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS", ], aes(x = I_relcov_scaled, y = fit), color = "#006633") +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI =="GREAT PLAINS",], aes(x = I_relcov_scaled, y = fit), color = "#E69F00", linewidth = 1)  +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",],aes(x = I_relcov_scaled, y = fit), color = "#F0E442", linewidth = 1)  +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI =="NORTH AMERICAN DESERTS",],aes(x = I_relcov_scaled, y = fit), color = "#D55E00", linewidth = 1)  +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI =="NORTHERN FORESTS",], aes(x = I_relcov_scaled, y = fit), color = "#009E73", linewidth = 1)  +
    geom_smooth(data = out_full[out_full$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",], aes(x = I_relcov_scaled, y = fit), color = "#56B4E9", linewidth = 1) +  
    geom_smooth(data = out_main, aes(x = I_relcov_scaled, y = fit), color = "black", linewidth = 1) + 
    geom_ribbon(data = out_main, aes(x = I_relcov_scaled, ymin = lwr, ymax = upr), alpha = 0.2) + 
    labs(title = paste(df$Trait[1], "\n", "R2_m =", R2_marginal,"; R2_c =", R2_conditional),
         x = "Scaled invasion level", y = paste(df$Trait[1], "(CWM)")) +
    theme_bw() + 
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0))+
    scale_fill_manual(values = cb_palette) + 
    scale_color_manual(values = cb_palette) + 
    theme(legend.position = "none", 
          legend.text = element_text(size = 8), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank()) + 
    guides(fill = guide_legend(nrow = 3, byrow = TRUE), 
           color = guide_legend(override.aes = list(alpha = 1)))
return(plot)
}

# ----------------------------------- updated model plots (w ER level 4) -------------------------------------

mod_plot_function_updated <- function(ls){
  # data for plot 
  df = ls[[1]]
  R2_marginal <- ls[[2]]
  R2_conditional <- ls[[3]]
  mod = ls[[4]]
  trait = df$Trait[[1]]
  
  # prediction intervals (mod 1)
  xmin = min(df$I_relcov_scaled)
  xmax = max(df$I_relcov_scaled)
  
  # updated - predictions from model using model coefficients 
  new_dat_full_2 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                                EcoRegionLevelI = c(levels(as.factor(df$EcoRegionLevelI))))
  coefs <- coef(mod)
  
  preds_new_dat_full_2 = new_dat_full_2
  preds_new_dat_full_2$y <- ifelse(preds_new_dat_full_2$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS", 
                                   (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[1,2]]) + coefs$EcoRegionLevelI[[1,1]],
                                   ifelse(preds_new_dat_full_2$EcoRegionLevelI == "GREAT PLAINS", 
                                          (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[2,2]]) + coefs$EcoRegionLevelI[[2,1]],
                                          ifelse(preds_new_dat_full_2$EcoRegionLevelI == "MEDITERRANEAN CALIFORNIA", 
                                                 (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[3,2]]) + coefs$EcoRegionLevelI[[3,1]],
                                                 ifelse(preds_new_dat_full_2$EcoRegionLevelI == "NORTH AMERICAN DESERTS", 
                                                        (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[4,2]]) + coefs$EcoRegionLevelI[[4,1]],
                                                        ifelse(preds_new_dat_full_2$EcoRegionLevelI == "NORTHERN FORESTS", 
                                                               (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[5,2]]) + coefs$EcoRegionLevelI[[5,1]],
                                                               (preds_new_dat_full_2$I_relcov_scaled*coefs$EcoRegionLevelI[[6,2]]) + coefs$EcoRegionLevelI[[6,1]])))))
                                                               
                                                 
  preds_new_dat_full_2$y_exp <- exp(preds_new_dat_full_2$y)
  
  # overall trend line 
  mean_intercept <- mean(coefs$EcoRegionLevelI$`(Intercept)`)
  mean_slope <- mean(coefs$EcoRegionLevelI$I_relcov_scaled)
  
  overall_I_relcov_scaled <- seq(xmin,xmax,0.005)
  overall_ys <- (overall_I_relcov_scaled*mean_slope) + mean_intercept
  y_exp <- exp(overall_ys)
                                  
  overall_line_dat <- data.frame(I_relcov_scaled = overall_I_relcov_scaled, 
                                 EcoRegionLevelI = "Overall", 
                                 y = overall_ys, 
                                 y_exp = y_exp)
  
  
  preds_new_dat_full_2 <- rbind(preds_new_dat_full_2, overall_line_dat)
  
  # data from model 
  if(trait == "SLA_mm2/mg"| trait == "heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_rooting_depth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    preds_new_dat_full_2$y_plot <- preds_new_dat_full_2$y_exp
  }
  
  if(trait != "SLA_mm2/mg"& trait != "heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_rooting_depth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    preds_new_dat_full_2$y_plot <- preds_new_dat_full_2$y
  }
  
  plot <- ggplot() +
    geom_point(data = df, aes(x = I_relcov_scaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.8, alpha = 0.4) +
    geom_point(data = preds_new_dat_full_2, aes(x = I_relcov_scaled, y = y_plot, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.4) +
    labs(title = paste(df$Trait[1], "\n", "R2_m =", R2_marginal,"; R2_c =", R2_conditional),
         x = "Scaled invasion level", y = paste(df$Trait[1], "(CWM)")) +
    theme_bw() + 
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0))+
    scale_color_manual(values = cb_palette) + 
    scale_fill_manual(values = cb_palette) +
    theme(legend.position = "none", 
          legend.text = element_text(size = 8), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank()) + 
    guides(fill = guide_legend(nrow = 3, byrow = TRUE), 
           color = guide_legend(override.aes = list(alpha = 1)))
  return(plot)
}

height_mod_plot <- mod_plot_function_updated(height_mod) + scale_y_continuous(trans = "log2", limits = c(0.04, 46))
maxheight_mod_plot <- mod_plot_function_updated(maxheight_mod) + scale_y_continuous(trans = "log2", limits = c(0.04, 46))
root_depth_mod_plot <- mod_plot_function_updated(root_depth_mod) + scale_y_continuous(trans = "log2", limits = c(0.12, 64))
root_avg_depth_mod_plot <- mod_plot_function_updated(root_avg_depth_mod) + scale_y_continuous(trans = "log2", limits = c(0.12, 32))

LDMC_mod_plot <- mod_plot_function_updated(LDMC_mod) + ylim(0.1, 0.6)
leaf_N_mod_plot <- mod_plot_function_updated(leaf_N_mod) + scale_y_continuous(trans = "log2")
leaf_P_mod_plot <- mod_plot_function_updated(leaf_P_mod) + scale_y_continuous(trans = "log2")
SLA_mod_plot <- mod_plot_function_updated(SLA_mod) + scale_y_continuous(trans = "log2")

root_diam_mod_plot <- mod_plot_function_updated(root_diam_mod) + scale_y_continuous(trans = "log2", limits = c(0.10, 1.25))
root_N_mod_plot <- mod_plot_function_updated(root_N_mod) + scale_y_continuous(trans = "log2", limits = c(2.5, 40))
RTD_mod_plot <- mod_plot_function_updated(RTD_mod) + ylim (0, 0.6)
SRL_mod_plot <- mod_plot_function_updated(SRL_mod) + scale_y_continuous(trans = "log2", limits = c(6,520))

# -------------------- data prep - weighted avgs (invasives vs. natives) ------------------------------

# merge datasets to select plots from dat that have >80 for specific traits 
dat_80_cov <- merge(dat, CWM_plot_info_80, by.x = c("Plot", "TraitNameAbr", "Year"), 
                    by.y = c("Plot", "Trait", "Year"))

dat_80_cov$NativeStatus <- gsub("NI", "N", dat_80_cov$NativeStatus)
dat_80_cov <- dat_80_cov[!dat_80_cov$TraitNameAbr %in% c("Duration", "Growth.Habit", "Mycorrhizal.type"),]

dat_80_splits <- split(dat_80_cov, list(dat_80_cov$Plot), drop = TRUE) 

# return just plots with both N and I in them 

inv_nat_present <- function(df){
  inv_levels_length <- length(levels(as.factor(df$NativeStatus)))
  
  if(inv_levels_length == 2){
    return(df)
  }
}

N_I_plots <- lapply(dat_80_splits, inv_nat_present)
N_I_plots[sapply(N_I_plots, is.null)] <- NULL

CWM_function_N_I <- function(df){
  df_splits <- split(df, df$TraitNameAbr)
  
  splits_function <- function(df_out){
    df_n <- df_out[df_out$NativeStatus == "N",]
    df_i <- df_out[df_out$NativeStatus == "I",]
    
    cwm_function_2 <- function(df){
      df$mean <- as.numeric(df$mean)
      tot_cov = sum(df$PctCov_100)
      df$rel_100_ss = df$PctCov_100/tot_cov * 100
      df$weighted_mean = (df$mean * df$rel_100_ss)/100
      cwm = sum(df$weighted_mean)
    }
    
    cwm_i_n <- as.data.frame(do.call(cbind, lapply(list(df_i, df_n), cwm_function_2)))
    colnames(cwm_i_n) <- c("CWM_I", "CWM_N")
    
    cwm_i_n$Plot <- df_out$Plot[1]
    cwm_i_n$Year <- df_out$Year[1]
    cwm_i_n$TraitNameAbr <- df_out$TraitNameAbr[1]
    cwm_i_n$inv_relcov <- sum(df_out[df_out$NativeStatus == "I",]$rel_100)
    cwm_i_n$nat_relcov <- sum(df_out[df_out$NativeStatus == "N",]$rel_100)
    return(cwm_i_n)
  }
  out <- do.call(rbind, lapply(df_splits, splits_function))
  return(out)
}

inv_nat_cwm <- do.call(rbind, lapply(N_I_plots, CWM_function_N_I))

# drop plots that do not have weighted avgs for both natives and invasive 
inv_nat_cwm <- inv_nat_cwm[!(inv_nat_cwm$inv_relcov == 0 | inv_nat_cwm$nat_relcov == 0),]

# merge with ecoregion 
eco_r <- unique(CWM_plot_info_80_2[c("Plot", "Year", "EcoRegionLevelI", "US_L4NAME")])
inv_nat_cwm_er <- merge(inv_nat_cwm, eco_r, by = c("Plot", "Year"))

# split by trait before models and plots 
inv_nat_cwm_er_splits <- split(inv_nat_cwm_er, inv_nat_cwm_er$TraitNameAbr)

# need averages for the different ecoregions? 
# before these were model estimates 
plot_dat <- function(df){
  Trait = df$TraitNameAbr[1]
  dat_long <- reshape(df, 
                      varying = c("CWM_I", "CWM_N"),
                      v.names = "CWM",
                      idvar = c("Plot"),
                      timevar = "NativeStatus",
                      times = c("I", "N"),
                      direction = "long")
  
  dat_long <- dat_long[!dat_long$CWM == -Inf,]
  
  if(Trait == "SLA_mm2/mg"| Trait == "heightveg_m" | Trait == "leafP_mg/g" | 
     Trait == "leafarea_mm2" | Trait == "leafN_mg/g" | Trait == "max_rooting_depth_m"| 
     Trait == "Mean_Root_diameter" | Trait == "Root_N_concentation"| Trait =="Specific_root_length"){
  mod <- lmer(log(CWM) ~ NativeStatus*EcoRegionLevelI + (1|Plot),
              data = dat_long)
  
  est_means_out <- as.data.frame(emmeans(mod, ~EcoRegionLevelI+NativeStatus))
  
  est_means_out_wide <- reshape(est_means_out,  idvar = "EcoRegionLevelI",
                                timevar = "NativeStatus",
                                direction = "wide")
  
  est_means_out_wide <- as.data.frame(est_means_out_wide)
  est_means_out_wide$CWM_I <- exp(est_means_out_wide$emmean.I)
  est_means_out_wide$CWM_N <- exp(est_means_out_wide$emmean.N)
  }
  
  if(Trait != "SLA_mm2/mg"&Trait != "heightveg_m" & Trait != "leafP_mg/g" & 
     Trait != "leafarea_mm2" & Trait != "leafN_mg/g" & Trait != "max_rooting_depth_m"& 
     Trait != "Mean_Root_diameter" & Trait != "Root_N_concentation"& Trait !="Specific_root_length"){
    mod <- lmer(CWM ~ NativeStatus*EcoRegionLevelI + (1|Plot),
                data = dat_long)
    est_means_out <- as.data.frame(emmeans(mod, ~EcoRegionLevelI+NativeStatus))
    est_means_out_wide <- reshape(est_means_out,  idvar = "EcoRegionLevelI",
                                  timevar = "NativeStatus",
                                  direction = "wide")
    est_means_out_wide$CWM_I <- est_means_out_wide$emmean.I
    est_means_out_wide$CWM_N <- est_means_out_wide$emmean.N
    
  }
  return(est_means_out_wide)
} 

# need averages for the different ecoregions? 
# update: just calculate means and sd 
plot_dat_updated <- function(df){
  Trait = df$TraitNameAbr[1]
  dat_long <- reshape(df, 
                      varying = c("CWM_I", "CWM_N"),
                      v.names = "CWM",
                      idvar = c("Plot"),
                      timevar = "NativeStatus",
                      times = c("I", "N"),
                      direction = "long")
  
  dat_long <- dat_long[!dat_long$CWM == -Inf,]
  
  # split by ecoregion to get means of weighted averages for native and invasion in each ER
  
  splits <- split(dat_long, dat_long$EcoRegionLevelI)
  
  ER_means_sd <- function(df_ER){
    df_ER_split <- split(df_ER, df_ER$NativeStatus)
    mean_sd_function <- function(df_ss){
      out <- data.frame(
        TraitNameAbr = df_ss$TraitNameAbr[1],
        EcoRegionLevelI = df_ss$EcoRegionLevelI[1],
        NativeStatus = df_ss$NativeStatus[1],
        weighted_mean = mean(df_ss$CWM, na.rm = TRUE),
        sd = sd(df_ss$CWM, na.rm = TRUE))
      return(out)
    }
    out_both <- do.call(rbind, lapply(df_ER_split, mean_sd_function))
    return(out_both)
  }
  out_all <- do.call(rbind, lapply(splits, ER_means_sd))
  out_all_wide <- reshape(out_all, 
                          idvar = c("EcoRegionLevelI", "TraitNameAbr"),
                          timevar = "NativeStatus", direction = "wide")
  return(out_all_wide)
} 

# get data for plots 
height <- plot_dat_updated(inv_nat_cwm_er_splits$heightveg_m)
root_depth <- plot_dat_updated(inv_nat_cwm_er_splits$max_rooting_depth_m) 
root_avg_depth <- plot_dat_updated(inv_nat_cwm_er_splits$mean_rooting_depth_m) 

SLA <- plot_dat_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
LDMC <- plot_dat_updated(inv_nat_cwm_er_splits$`LDMC_g/g`)
leafN <- plot_dat_updated(inv_nat_cwm_er_splits$`leafN_mg/g`)
leafP <- plot_dat_updated(inv_nat_cwm_er_splits$`leafP_mg/g`)

RTD <- plot_dat_updated(inv_nat_cwm_er_splits$Root_tissue_density)
Root_diam <- plot_dat_updated(inv_nat_cwm_er_splits$Mean_Root_diameter)
RootN <- plot_dat_updated(inv_nat_cwm_er_splits$Root_N_concentration)
SRL <- plot_dat_updated(inv_nat_cwm_er_splits$Specific_root_length)

# plotting 
plot_function <- function(df, means_df, trait){
  n_plots = nrow(df)
  plot <-   ggplot(df, aes (x = CWM_N, y = CWM_I)) + 
    geom_point(aes(color = EcoRegionLevelI), size = 1, alpha = 0.3) + 
    geom_point(data = means_df,
               aes(x = CWM_N, y = CWM_I, fill = EcoRegionLevelI),
               shape = 21, size = 5, alpha = 1) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9"))+
    scale_fill_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9")) +
    ylab(paste("Introduced Weighted Average", trait)) + 
    xlab(paste("Native Weighted Average", trait)) + 
    labs(title = paste(trait, "\n", "n_plots =", n_plots)) + 
    theme_bw() + 
    theme(legend.position = "none")
  return(plot)
}


# plotting 
plot_function_updated <- function(df, means_df, trait){
  n_plots = nrow(df)
  plot <-   ggplot(df) + 
    geom_point(aes(x = CWM_N, y = CWM_I, color = EcoRegionLevelI), size = 1, alpha = 0.3) + 
    geom_point(data = means_df,
               aes(x = weighted_mean.N, y = weighted_mean.I, fill = EcoRegionLevelI),
               shape = 21, size = 5, alpha = 1) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9"))+
    scale_fill_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9")) +
    geom_errorbar(data = means_df, aes(x = weighted_mean.N,
                                       ymin = weighted_mean.I - sd.I, 
                                       ymax = weighted_mean.I + sd.I), width = 0) + 
    geom_errorbarh(data = means_df, aes( y = weighted_mean.I, 
                                      xmin = weighted_mean.N - sd.N, 
                                      xmax = weighted_mean.N + sd.N), height = 0) + 
    ylab(paste("Introduced Weighted Average - ", trait)) + 
    xlab(paste("Native Weighted Average - ", trait)) + 
    labs(title = paste(trait, "\n", "n_plots =", n_plots)) + 
    theme_bw() + 
    theme(legend.position = "none")
  return(plot)
}

height_plot <- plot_function_updated(inv_nat_cwm_er_splits$heightveg_m, height, "Height") + 
  scale_y_continuous( limits = c(-5, 40), expand = c(0,0)) + 
  scale_x_continuous( limits = c(-5, 40), expand = c(0,0))

SLA_plot <- plot_function_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA, "SLA_mm2/mg")+ 
  scale_y_continuous(limits = c(4, 60), expand = c(0,0)) + 
  scale_x_continuous( limits = c(4, 60), expand = c(0,0))

LDMC_plot <- plot_function_updated(inv_nat_cwm_er_splits$`LDMC_g/g`, LDMC,  "LDMC_g/g") + 
  scale_y_continuous(limits = c(0.08, 0.5), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.08, 0.5), expand = c(0,0))

leafN_plot <- plot_function_updated(inv_nat_cwm_er_splits$`leafN_mg/g`, leafN, "leafN_mg/g") + 
  scale_y_continuous(limits = c(8, 50), expand = c(0,0)) + 
  scale_x_continuous(limits = c(8, 50), expand = c(0,0))

leafP_plot <- plot_function_updated(inv_nat_cwm_er_splits$`leafP_mg/g`, leafP, "leafP_mg/g") + 
  scale_y_continuous(limits = c(1, 8), expand = c(0,0)) + 
  scale_x_continuous(limits = c(1, 8), expand = c(0,0))

root_depth_plot <- plot_function_updated(inv_nat_cwm_er_splits$max_rooting_depth_m, root_depth, "max_rooting_depth")+ 
  scale_y_continuous(limits = c(-2, 30), expand = c(0,0)) + 
  scale_x_continuous(limits = c(-2, 30), expand = c(0,0))

root_avg_depth_plot <- plot_function_updated(inv_nat_cwm_er_splits$mean_rooting_depth_m, root_avg_depth, "rooting_avg_depth")+ 
  scale_y_continuous(limits = c(-2, 15), expand = c(0,0)) + 
  scale_x_continuous(limits = c(-2, 15), expand = c(0,0))

RTD_plot <- plot_function_updated(inv_nat_cwm_er_splits$Root_tissue_density, RTD, "RTD")  + 
  scale_y_continuous(limits = c(0.01, 0.75), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.01, 0.75), expand = c(0,0))

Root_diam_plot <- plot_function_updated(inv_nat_cwm_er_splits$Mean_Root_diameter, Root_diam, "Root diameter")  + 
  scale_y_continuous(limits = c(-0.25, 4), expand = c(0,0)) + 
  scale_x_continuous( limits = c(-0.25, 4), expand = c(0,0))

RootN_plot <- plot_function_updated(inv_nat_cwm_er_splits$Root_N_concentration, RootN, "RootN")  + 
  scale_y_continuous(limits = c(5, 25), expand = c(0,0)) + 
  scale_x_continuous( limits = c(5, 25), expand = c(0,0))

SRL_plot <- plot_function_updated(inv_nat_cwm_er_splits$Specific_root_length, SRL, "SRL")+ 
  scale_y_continuous(limits = c(10, 660), expand = c(0,0)) + 
  scale_x_continuous(limits = c(10, 660), expand = c(0,0))

legend <- get_legend(height_plot)

# updated models for weighted averages between native and non-native species WITHIN plots 
# Need to get CWM into one column 

weighted_avg_mod <- function(df){
  TraitNameAbr <- df$TraitNameAbr[1]
  dat_long <- reshape(df, 
                      varying = c("CWM_I", "CWM_N"),
                      v.names = "CWM",
                      idvar = c("Plot"),
                      timevar = "NativeStatus",
                      times = c("I", "N"),
                      direction = "long")
  
  
  # models
  if(TraitNameAbr == "SLA_mm2/mg"| TraitNameAbr == "heightveg_m" | TraitNameAbr == "leafP_mg/g" | 
     TraitNameAbr == "leafarea_mm2" | TraitNameAbr == "leafN_mg/g" | TraitNameAbr == "max_rooting_depth_m"| 
     TraitNameAbr == "Mean_Root_diameter" | TraitNameAbr == "Root_N_concentation"| TraitNameAbr =="Specific_root_length"){
    mod_out <- lmer(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME), data = dat_long)
  }
  
  if(TraitNameAbr != "SLA_mm2/mg"& TraitNameAbr != "heightveg_m" & TraitNameAbr != "leafP_mg/g" & 
     TraitNameAbr != "leafarea_mm2" & TraitNameAbr != "leafN_mg/g" & TraitNameAbr != "max_rooting_depth_m"& 
     TraitNameAbr != "Mean_Root_diameter" & TraitNameAbr != "Root_N_concentation"& TraitNameAbr !="Specific_root_length"){
    mod_out <- lmer(CWM ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME), data = dat_long)
  }
  
  sum_mod <- summary(mod_out)
  return(sum_mod)
} 

weighted_avg_mod(inv_nat_cwm_er_splits$heightveg_m)
weighted_avg_mod(inv_nat_cwm_er_splits$max_heightveg_m)
weighted_avg_mod(inv_nat_cwm_er_splits$max_rooting_depth_m)
weighted_avg_mod(inv_nat_cwm_er_splits$mean_rooting_depth_m)

weighted_avg_mod(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
weighted_avg_mod(inv_nat_cwm_er_splits$`LDMC_g/g`)
weighted_avg_mod(inv_nat_cwm_er_splits$`leafN_mg/g`)
weighted_avg_mod(inv_nat_cwm_er_splits$`leafP_mg/g`)

weighted_avg_mod(inv_nat_cwm_er_splits$Specific_root_length)
weighted_avg_mod(inv_nat_cwm_er_splits$Mean_Root_diameter)
weighted_avg_mod(inv_nat_cwm_er_splits$Root_N_concentration)
weighted_avg_mod(inv_nat_cwm_er_splits$Root_tissue_density)
weighted_avg_mod(inv_nat_cwm_er_splits$Root_P_concentration)


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_size_plots_June2023.pdf", height = 10, width = 12)
plot_grid(plot_grid(height_veg_densityplot, height_mod_plot, height_plot,
                    root_depth_densityplot, root_depth_mod_plot, root_depth_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_leaf_plots_June2023.pdf", height = 10, width = 12)
plot_grid(plot_grid(SLA_densityplot, SLA_mod_plot, SLA_plot,
                    LDMC_densityplot, LDMC_mod_plot, LDMC_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_leafnutrients_plots_June2023.pdf", height = 10, width = 12)
plot_grid(plot_grid(leaf_N_densityplot, leaf_N_mod_plot, leafN_plot,
                    leaf_P_densityplot, leaf_P_mod_plot, leafP_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_rootcollab_plots_June2023.pdf", height = 10, width = 12)
plot_grid(plot_grid(SRL_densityplot, SRL_mod_plot, SRL_plot,
                    root_diam_densityplot, root_diam_mod_plot, Root_diam_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_rootcons_plots_June2023.pdf", height = 10, width = 12)
plot_grid(plot_grid(root_N_densityplot, root_N_mod_plot, RootN_plot,
                    RTD_densityplot, RTD_mod_plot, RTD_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

# Not using or need to review 
pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_size_plots_2.pdf", height = 10, width = 12)
plot_grid(plot_grid(height_veg_densityplot, height_mod_plot, height_plot,
                    root_depth_avg_densityplot, root_depth_mod_plot, root_avg_depth_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()
