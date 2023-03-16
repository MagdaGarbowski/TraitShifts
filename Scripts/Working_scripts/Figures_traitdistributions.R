# Trait shifts 
# plotting! 

library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(cowplot)
library(MuMIn)
library(gridExtra)
library(merTools)

# ----------------------------------- load data -------------------------------------

CWM <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM.csv", 
                           select = c("Plot", "Year", "TraitNameAbr", "CWM")))

colnames(CWM) <- c("Plot", "Year", "Trait", "CWM")

plot_info <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

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

# get time-series dataset (for later use?)

time_series_function <- function(df){
  years = length(levels(as.factor(df$Year)))
  if(years > 1){
    return(df)
  }
}

plot_info_time_series <- do.call(rbind, lapply(plot_info_splits, time_series_function))

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

CWM_plot_info_80$Inv_level <- ifelse(CWM_plot_info_80$I_relcov < 2, "low", 
                                  ifelse(CWM_plot_info_80$I_relcov >= 2 & CWM_plot_info_80$I_relcov < 20, "med", "high"))

CWM_plot_info_80_2percent$Inv_level <- ifelse(CWM_plot_info_80_2percent$I_relcov < 15, "low", 
                                           ifelse(CWM_plot_info_80_2percent$I_relcov >= 15 & CWM_plot_info_80_2percent$I_relcov < 45, "med", "high"))

# ----------------------------------- density plots  -------------------------------------

cb_palette <- c("#006633",  "#E69F00", "#56B4E9","#D55E00", "#F0E442","#009E73",  "#0072B2",  "#CC79A7","#999999" )


# split by trait 
CWM_80_splits <- split(CWM_plot_info_80, CWM_plot_info_80$Trait)

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
    theme(legend.position = "bottom") 
  return(plot_out)
}

# size gradient
height_veg <- density_plot_function(CWM_80_splits$heightveg_m, "Height-veg (m)", 0, 3, c("low (<2%)", "med (2-20%)", "high (>20%)"))
root_depth <- density_plot_function(CWM_80_splits$max_rooting_depth_m, "max_root_depth (m)", 0, 12, c("low (<2%)", "med (2-20%)", "high (>20%)"))

# leaf conservation 
SLA <- density_plot_function(CWM_80_splits$`SLA_mm2/mg`, "SLA (mm2 mg-1)", 0, 50, c("low (<2%)", "med (2-20%)", "high (>20%)"))
LDMC <- density_plot_function(CWM_80_splits$`LDMC_g/g`, "LDMC (g g-1)", 0, 0.6, c("low (<2%)", "med (2-20%)", "high (>20%)"))
leaf_N <- density_plot_function(CWM_80_splits$`leafN_mg/g`, "leaf_N (mg/g)", 0, 40, c("low (<2%)", "med (2-20%)", "high (>20%)"))
leaf_P <- density_plot_function(CWM_80_splits$`leafP_mg/g`, "leaf_P (mg/g)", 0, 4, c("low (<2%)", "med (2-20%)", "high (>20%)"))

# roots - conservation and collaboration 
root_diam <- density_plot_function(CWM_80_splits$Mean_Root_diameter, "root_diam (mm)", 0, 1, c("low (<2%)", "med (2-20%)", "high (>20%)"))
root_N <- density_plot_function(CWM_80_splits$Root_N_concentration, "root_N", 0, 20, c("low (<2%)", "med (2-20%)", "high (>20%)"))
RTD <- density_plot_function(CWM_80_splits$Root_tissue_density, "RTD", 0, 0.75, c("low (<2%)", "med (2-20%)", "high (>20%)"))
SRL <- density_plot_function(CWM_80_splits$Specific_root_length, "SRL", 0, 350, c("low (<2%)", "med (2-20%)", "high (>20%)"))

# other 
seed_mass <- density_plot_function(CWM_80_splits$seedmass_mg, "seed_mass (mg)", 0, 10, c("low (<2%)", "med (2-20%)", "high (>20%)"))
SSD <- density_plot_function(CWM_80_splits$`SSD_g/cm3`, "SSD", 0, 1, c("low (<2%)", "med (2-20%)", "high (>20%)"))
RDMC <- density_plot_function(CWM_80_splits$Root_dry_matter_content, "RDMC", 0, 0.45, c("low (<2%)", "med (2-20%)", "high (>20%)"))
RMF <- density_plot_function(CWM_80_splits$Root_mass_fraction, "RMF", 0, 0.75, c("low (<2%)", "med (2-20%)", "high (>20%)"))
woodiness <- density_plot_function(CWM_80_splits$Woodiness, "Woodiness", 0, 1, c("low (<2%)", "med (2-20%)", "high (>20%)"))
annual <- density_plot_function(CWM_80_splits$Annual, "Annual_prop", 0, 1, c("low (<2%)", "med (2-20%)", "high (>20%)"))
leaf_area <- density_plot_function(CWM_80_splits$leafarea_mm2, "Leaf Area (mm2)", 0, 3000, c("low (<2%)", "med (2-20%)", "high (>20%)"))
root_P <- density_plot_function(CWM_80_splits$Root_P_concentration, "root_P", 0, 4, c("low (<2%)", "med (2-20%)", "high (>20%)"))


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_1.pdf", width = 10, height = 6)
grid.arrange(height_veg, leaf_area, SSD, seed_mass, ncol = 2)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_2.pdf", width = 10, height = 6)
grid.arrange(SLA, LDMC, leaf_N, leaf_P, ncol = 2)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_3.pdf", width = 10, height = 6)
grid.arrange(SRL, RDMC,root_N, root_P,  ncol = 2)
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_overall_4.pdf", width = 10, height = 6)
grid.arrange(root_depth, RMF, RTD, ncol = 2)
dev.off()

# --------------------------------- models -------------------------------------

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

  R2_marginal <- round(r.squaredGLMM(mod)[1],3)
  R2_conditional <- round(r.squaredGLMM(mod)[2],3)
  
  ls <- list(dat, R2_marginal, R2_conditional, mod)
  return(ls)
}


height_mod <- mod_function(CWM_80_splits$heightveg_m)
root_depth <- mod_function(CWM_80_splits$max_rooting_depth_m)

SLA_mod <- mod_function(CWM_80_splits$`SLA_mm2/mg`)
LDMC_mod <- mod_function(CWM_80_splits$`LDMC_g/g`)
leaf_N_mod <- mod_function(CWM_80_splits$`leafN_mg/g`) # boundary fit singular 
leaf_P_mod <- mod_function(CWM_80_splits$`leafP_mg/g`)

root_diam_mod <- mod_function(CWM_80_splits$Mean_Root_diameter)
RTD_mod <- mod_function(CWM_80_splits$Root_tissue_density)
root_N_mod <- mod_function(CWM_80_splits$Root_N_concentration)
SRL_mod <- mod_function(CWM_80_splits$Specific_root_length)

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
  
  # effects coeffs
  coeffs_rand <- coef(mod)$EcoRegionLevelI
  coeffs_rand$EcoRegionLevelI <- rownames(coeffs_rand)
  coeffs_mod <- summary(mod)$coefficients
  
  # new data for predictions
  new_dat <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.00025), 
                         EcoRegionLevelI = c(levels(as.factor(df$EcoRegionLevelI))))

  predict_int <- predictInterval(mod, new_dat, n.sims = 2000,
                                   stat = "median",
                                   type = "linear.prediction",
                                   which = "fixed",
                                   include.resid.var = TRUE)

  # data from model 
  if(trait == "SLA_mm2/mg"| trait == "heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_rooting_depth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    predicted_dat_exp <- exp(predicted_dat)
    out <- cbind(new_dat, predict_int)
  }
  
  if(trait != "SLA_mm2/mg"& trait != "heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_rooting_depth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    out <- do.call(cbind, list(new_dat, predict_int))
  }
  
  plot <- ggplot() +
    geom_point(data = df, aes(x = I_relcov_scaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), alpha = 0.2) +
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="EASTERN TEMPERATE FORESTS",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="EASTERN TEMPERATE FORESTS",]$`(Intercept)`), color = "#006633", linewidth = 1) + 
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="GREAT PLAINS",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="GREAT PLAINS",]$`(Intercept)`), color = "#E69F00", linewidth = 1)  +
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",]$`(Intercept)`), color = "#F0E442", linewidth = 1)  +
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTH AMERICAN DESERTS",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTH AMERICAN DESERTS",]$`(Intercept)`), color = "#D55E00", linewidth = 1)  +
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTHERN FORESTS",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTHERN FORESTS",]$`(Intercept)`), color = "#009E73", linewidth = 1)  +
    geom_abline(slope = coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",]$I_relcov_scaled, 
                intercept = (coeffs_rand[coeffs_rand$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",]$`(Intercept)`), color = "#56B4E9", linewidth = 1)  +
    geom_abline(slope = coeffs_mod[2], intercept = coeffs_mod[1], color = "black", linewidth = 1.5) +   
    geom_ribbon(data = out, aes(x = I_relcov_scaled, ymin = lwr, ymax = upr), alpha = 0.2) + 
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
           color = guide_legend(override.aes = list(alpha = 1))) + scale_y_continuous(trans = "log2")
  plot
return(plot)
}

height_mod <- mod_plot_function(height_mod) + scale_y_continuous(trans = "log2") 
LDMC_mod <- mod_plot_function(LDMC_mod)
leaf_area_mod <- mod_plot_function(mods_out$leafarea_mm2, -0.65, 3.6) + scale_y_continuous(trans = "log2")
leaf_N_mod <- mod_plot_function(mods_out$`leafN_mg/g`, -0.65, 3.6) + scale_y_continuous(trans = "log2")
leaf_P_mod <- mod_plot_function(mods_out$`leafP_mg/g`, -0.65, 3.6) + scale_y_continuous(trans = "log2")
SLA_mod <- mod_plot_function(mods_out$`SLA_mm2/mg`, -0.65, 3.6) + scale_y_continuous(trans = "log2")
seed_mass_mod <- mod_plot_function(mods_out$seedmass_mg, -0.65, 3.6) + scale_y_continuous(trans = "log2")
SSD_mod <- mod_plot_function(mods_out$`SSD_g/cm3`, -0.65, 4.2)
root_diam_mod <- mod_plot_function(mods_out$Mean_Root_diameter, -0.65, 4.2) + scale_y_continuous(trans = "log2")
RDMC_mod <- mod_plot_function(mods_out$Root_dry_matter_content, -0.65, 2)
RMF_mod <- mod_plot_function(mods_out$Root_mass_fraction, -0.65, 2.5)
root_N_mod <- mod_plot_function(mods_out$Root_N_concentration,  -0.65, 2.8) + scale_y_continuous(trans = "log2")
root_P_mod <- mod_plot_function(mods_out$Root_P_concentration,  -0.65, 3.5) + scale_y_continuous(trans = "log2")
RTD_mod <- mod_plot_function(mods_out$Root_tissue_density, -0.65, 3.5)
root_depth_mod <- mod_plot_function(mods_out$max_rooting_depth_m, -0.65, 2.8) + scale_y_continuous(trans = "log2")
SRL_mod <- mod_plot_function(mods_out$Specific_root_length, -0.65, 3.2) + scale_y_continuous(trans = "log2", limits = c(4, 512))
annual_mod <- mod_plot_function(mods_out$Annual, -0.65, 4.2) + ylim(0, 1)
woody_mod <- mod_plot_function(mods_out$Woodiness, -0.65, 4.2) + ylim(0, 1)

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/TraitShifts_plots_01252023.pdf", width = 15, height = 8)
ggdraw(add_sub(plot_grid(SLA, SLA_mod, nrow = 1, align = "h", axis = "tb"),
               "SLA increases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(LDMC, LDMC_mod, nrow = 1, align = "h", axis = "tb"),
               "LDMC decreases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(leaf_N, leaf_N_mod, nrow = 1, align = "h", axis = "tb"),
               "Leaf N increases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(leaf_P, leaf_P_mod, nrow = 1, align = "h", axis = "tb"),
               "Leaf P increases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(height_veg, height_mod, nrow = 1, align = "h", axis = "tb"),
               "Height decreases with invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(SSD, SSD_mod, nrow = 1, align = "h", axis = "tb"),
               "SSD decreases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(RMF, RMF_mod, nrow = 1, align = "h", axis = "tb"),
               "RMF increases with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(root_diam, root_diam_mod, nrow = 1, align = "h", axis = "tb"),
               "Root diameter decreased with increasing invasion", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(root_N, root_N_mod, nrow = 1, align = "h", axis = "tb"),
               "No change", y = 0, vjust = -0.5))
ggdraw(add_sub(plot_grid(RTD, RTD_mod, nrow = 1, align = "h", axis = "tb"),
               "RTD is confusing (decreasing in some ecoregions?)", y = 0, vjust = -0.5))
density_eco_plot_function(CWM_80_splits$Root_tissue_density, "RTD", 0, 0.75)
ggdraw(add_sub(plot_grid(root_depth, root_depth_mod, nrow = 1, align = "h", axis = "tb"),
               "Max rooting depth is confusing", y = 0, vjust = -0.5))
density_eco_plot_function(CWM_80_splits$Rooting_depth, "root_depth (m)", 0, 4)
ggdraw(add_sub(plot_grid(SRL, SRL_mod, nrow = 1, align = "h", axis = "tb"),
               "SRL is very confusing", y = 0, vjust = -0.5))
density_eco_plot_function(CWM_80_splits$Specific_root_length, "SRL", 0, 350)
ggdraw(add_sub(plot_grid(seed_mass, seed_mass_mod, nrow = 1, align = "h", axis = "tb"),
               "Seed mass is very confusing", y = 0, vjust = -0.5))
density_eco_plot_function(CWM_80_splits$seedmass_mg, "seed_mass (mg)", 0, 20)
ggdraw(add_sub(plot_grid(leaf_area, leaf_area_mod, nrow = 1, align = "h", axis = "tb"),
               "Leaf area is confusing", y = 0, vjust = -0.5))
density_eco_plot_function(CWM_80_splits$leafarea_mm2, "leaf_area (mm2)", 0, 3000)

dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/TraitShifts_plots_2_01252023.pdf", width = 15, height = 12)
plot_grid(SLA_mod,leaf_N_mod,leaf_P_mod,RMF_mod,
          LDMC_mod,height_mod,SSD_mod,root_diam_mod,
          root_N_mod,RTD_mod,root_depth_mod,SRL_mod, ncol = 4, align = "h", axis = "l")
dev.off() 


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/TraitShifts_plots_3_01252023.pdf", width = 15, height = 12)
plot_grid(SLA_mod_ss,leaf_N_mod_ss,leaf_P_mod_ss,RMF_mod_ss,
          LDMC_mod_ss,height_mod_ss,SSD_mod_ss,root_diam_mod_ss,
          root_N_mod_ss,RTD_mod_ss,root_depth_mod_ss,SRL_mod_ss, ncol = 4, align = "h", axis = "l")
dev.off()  


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/traits_proportions.pdf", width = 12, height = 10)
plot_grid(woodiness, annual, woody_mod, annual_mod, ncol = 2, align = "v", axis = "b")
dev.off()

#
#
#
#
#
#


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/height.pdf", width = 15, height = 7)
plot_grid(height_veg, height_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/leaf_N.pdf", width = 15, height = 7)
plot_grid(leaf_N, leaf_N_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/root_diam.pdf", width = 15, height = 7)
plot_grid(root_diam, root_diam_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/SRL.pdf", width = 15, height = 7)
plot_grid(SRL, SRL_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/rooting_depth.pdf", width = 15, height = 7)
plot_grid(root_depth, root_depth_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/SLA.pdf", width = 15, height = 7)
plot_grid(SLA, SLA_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/RTD.pdf", width = 15, height = 7)
plot_grid(RTD, RTD_mod, ncol = 2, align = "h")
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Feb_2023_meeting/six_traits.pdf", width = 15, height = 10)
plot_grid(leaf_N_mod, SLA_mod, height_mod,RTD_mod, SRL_mod, root_depth_mod, ncol = 3, align = "h")
dev.off()




# level of invasion by ecosystem 
head(plot_info)
plot_ss <- plot_info[c("Plot", "Year", "I_relcov", "EcoRegionLevelI")]
plot_ss <- plot_ss[!duplicated(plot_ss),]

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/Inv_level.pdf", width = 15, height = 12)
ggplot(plot_ss, aes(x = I_relcov)) + 
  geom_density(adjust = 1.5, alpha = 0.4) + 
  facet_wrap(~EcoRegionLevelI, scales = "free")




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

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
           

# ----------------------------------- traits by eco-region counts  -------------------------------------
CWM_80_trait_splits <- split(CWM_80, list(CWM_80$Trait, CWM_80$EcoRegionLevelI), drop = TRUE)

trait_by_eco_function <- function(df){
  df_unique <- unique(df[c("Plot", "Trait", "EcoRegionLevelI", "Inv_level")])
  agg_dat <- as.data.frame(aggregate(Plot ~ Inv_level, df_unique, length))
  all_df <- data.frame(Inv_level = "all", 
                       Plot = sum(agg_dat$Plot))
  agg_dat_out <- rbind(agg_dat, all_df)
  agg_dat_out$Plot <- as.numeric(agg_dat_out$Plot)
  agg_dat_out$Trait <- df$Trait[1]
  agg_dat_out$EcoRegionLevelI <- df$EcoRegionLevelI[1]
  return(agg_dat_out)
}


trait_eco <- do.call(rbind, lapply(CWM_80_trait_splits, trait_by_eco_function))
trait_eco_wide <- reshape(trait_eco, idvar = c("EcoRegionLevelI","Trait"), timevar = "Inv_level", direction = "wide")
colnames(trait_eco_wide) <- gsub("Plot.", "", colnames(trait_eco_wide))
trait_eco_wide <- trait_eco_wide[c("EcoRegionLevelI","Trait", "all", "low", "med", "high")]

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


df = CWM_80_splits$`SLA_mm2/mg`[CWM_80_splits$`SLA_mm2/mg`$EcoRegionLevelI == "TEMPERATE SIERRAS",]
df = df[df$Inv_level %in% c("high", "low"),]
ggplot(df, aes(x = CWM, group = Inv_level, fill = Inv_level)) + 
  geom_density(adjust = 1.5, alpha = 0.4) + 
  scale_fill_manual(name = "Invasion level", 
                    breaks = c("low",  "high"),
                    labels = c("low (<33%)",  "high (>66%)"),
                    values = c("#999999",  "#E69F00")) + 
  xlim(0, 40) + 
  labs(x = "CWM - SLA") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16))

        