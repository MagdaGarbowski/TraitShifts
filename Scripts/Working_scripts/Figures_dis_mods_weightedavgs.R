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

dat <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wtraits.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean")))

plot_info <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

colnames(CWM) <- c("Plot", "Year", "Trait", "CWM")


cb_palette <- c("#006633",  "#E69F00", "#56B4E9","#D55E00", "#F0E442","#009E73",  "#0072B2",  "#CC79A7","#999999" )

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

CWM_plot_info_80$Inv_level <- ifelse(CWM_plot_info_80$I_relcov < 10, "low", 
                                  ifelse(CWM_plot_info_80$I_relcov >= 10 & CWM_plot_info_80$I_relcov < 30, "med", "high"))

CWM_plot_info_80_2percent$Inv_level <- ifelse(CWM_plot_info_80_2percent$I_relcov < 15, "low", 
                                           ifelse(CWM_plot_info_80_2percent$I_relcov >= 15 & CWM_plot_info_80_2percent$I_relcov < 45, "med", "high"))

# ----------------------------------- density plots  -------------------------------------

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
    theme(legend.position = "none") 
  return(plot_out)
}

# size gradient
height_veg_densityplot <- density_plot_function(CWM_80_splits$heightveg_m, "Height-veg (m)", 0, 6, c("low (<10%)", "med (10-30%)", "high (>30%)"))
root_depth_densityplot <- density_plot_function(CWM_80_splits$max_rooting_depth_m, "max_root_depth (m)", 0, 12, c("low (<10%)", "med (10-30%)", "high (>30%)"))

# leaf conservation 
SLA_densityplot <- density_plot_function(CWM_80_splits$`SLA_mm2/mg`, "SLA (mm2 mg-1)", 0, 50, c("low (<10%)", "med (10-30%)", "high (>30%)"))
LDMC_densityplot <- density_plot_function(CWM_80_splits$`LDMC_g/g`, "LDMC (g g-1)", 0, 0.6, c("low (<10%)", "med (10-30%)", "high (>30%)"))
leaf_N_densityplot <- density_plot_function(CWM_80_splits$`leafN_mg/g`, "leaf_N (mg/g)", 0, 40, c("low (<10%)", "med (10-30%)", "high (>30%)"))
leaf_P_densityplot <- density_plot_function(CWM_80_splits$`leafP_mg/g`, "leaf_P (mg/g)", 0, 4, c("low (<10%)", "med (10-30%)", "high (>30%)"))

# roots - conservation and collaboration 
root_diam_densityplot <- density_plot_function(CWM_80_splits$Mean_Root_diameter, "root_diam (mm)", 0, 1, c("low (<10%)", "med (10-30%)", "high (>30%)"))
root_N_densityplot <- density_plot_function(CWM_80_splits$Root_N_concentration, "root_N", 0, 20, c("low (<10%)", "med (10-30%)", "high (>30%)"))
RTD_densityplot <- density_plot_function(CWM_80_splits$Root_tissue_density, "RTD", 0, 0.75, c("low (<10%)", "med (10-30%)", "high (>30%)"))
SRL_densityplot <- density_plot_function(CWM_80_splits$Specific_root_length, "SRL", 0, 350, c("low (<10%)", "med (10-30%)", "high (>30%)"))

# other 
seed_mass_densityplot <- density_plot_function(CWM_80_splits$seedmass_mg, "seed_mass (mg)", 0, 10, c("low (<10%)", "med (10-30%)", "high (>30%)"))
SSD_densityplot <- density_plot_function(CWM_80_splits$`SSD_g/cm3`, "SSD", 0, 1, c("low (<10%)", "med (10-30%)", "high (>30%)"))
RDMC_densityplot <- density_plot_function(CWM_80_splits$Root_dry_matter_content, "RDMC", 0, 0.45, c("low (<10%)", "med (10-30%)", "high (>30%)"))
RMF_densityplot <- density_plot_function(CWM_80_splits$Root_mass_fraction, "RMF", 0, 0.75, c("low (<10%)", "med (10-30%)", "high (>30%)"))
woodiness_densityplot <- density_plot_function(CWM_80_splits$Woodiness, "Woodiness", 0, 1, c("low (<10%)", "med (10-30%)", "high (>30%)"))
annual_densityplot <- density_plot_function(CWM_80_splits$Annual, "Annual_prop", 0, 1, c("low (<10%)", "med (10-30%)", "high (>30%)"))
leaf_area_densityplot <- density_plot_function(CWM_80_splits$leafarea_mm2, "Leaf Area (mm2)", 0, 3000, c("low (<10%)", "med (10-30%)", "high (>30%)"))
root_P_densityplot <- density_plot_function(CWM_80_splits$Root_P_concentration, "root_P", 0, 4, c("low (<10%)", "med (10-30%)", "high (>30%)"))

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

  R2_marginal <- round(r.squaredGLMM(mod)[1],3)
  R2_conditional <- round(r.squaredGLMM(mod)[2],3)
  
  ls <- list(dat, R2_marginal, R2_conditional, mod)
  return(ls)
}


height_mod <- mod_function(CWM_80_splits$heightveg_m)
root_depth_mod <- mod_function(CWM_80_splits$max_rooting_depth_m)

SLA_mod <- mod_function(CWM_80_splits$`SLA_mm2/mg`)
LDMC_mod <- mod_function(CWM_80_splits$`LDMC_g/g`)
leaf_N_mod <- mod_function(CWM_80_splits$`leafN_mg/g`) # boundary fit singular 
leaf_P_mod <- mod_function(CWM_80_splits$`leafP_mg/g`)

root_diam_mod <- mod_function(CWM_80_splits$Mean_Root_diameter)
RTD_mod <- mod_function(CWM_80_splits$Root_tissue_density)
root_N_mod <- mod_function(CWM_80_splits$Root_N_concentration)
root_P_mod <- mod_function(CWM_80_splits$Root_P_concentration) # boundary fit singular 

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

  # new data for predictions
  new_dat_full <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.01), 
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

height_mod_plot <- mod_plot_function(height_mod) + scale_y_continuous(trans = "log2", limits = c(0.04, 46))
root_depth_mod_plot <- mod_plot_function(root_depth_mod) + scale_y_continuous(trans = "log2", limits = c(0.12, 64))

LDMC_mod_plot <- mod_plot_function(LDMC_mod) + ylim(0.1, 0.6)
leaf_N_mod_plot <- mod_plot_function(leaf_N_mod) + scale_y_continuous(trans = "log2")
leaf_P_mod_plot <- mod_plot_function(leaf_P_mod) + scale_y_continuous(trans = "log2")
SLA_mod_plot <- mod_plot_function(SLA_mod) + scale_y_continuous(trans = "log2")

root_diam_mod_plot <- mod_plot_function(root_diam_mod) + scale_y_continuous(trans = "log2", limits = c(0.10, 1.25))
root_N_mod_plot <- mod_plot_function(root_N_mod) + scale_y_continuous(trans = "log2", limits = c(2.5, 40))
root_P_mod_plot <- mod_plot_function(root_P_mod) + scale_y_continuous(trans = "log2", limits = c(1, 16))
RTD_mod_plot <- mod_plot_function(RTD_mod) + ylim (0, 0.6)
SRL_mod_plot <- mod_plot_function(SRL_mod) + scale_y_continuous(trans = "log2", limits = c(6,520))

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
inv_nat_cwm <- inv_nat_cwm[!(inv_nat_cwm$inv_relcov == 0 | inv_nat_cwm$nat_relcov ==0),]

# merge with ecoregion 
eco_r <- unique(CWM_plot_info_80[c("Plot", "Year", "EcoRegionLevelI")])
inv_nat_cwm_er <- merge(inv_nat_cwm, eco_r, by = c("Plot", "Year"))

# split by trait before models 
inv_nat_cwm_er_splits <- split(inv_nat_cwm_er, inv_nat_cwm_er$TraitNameAbr)

# need averages for the different ecoregions? 
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

# get data for plots 
height <- plot_dat(inv_nat_cwm_er_splits$heightveg_m)
root_depth <- plot_dat(inv_nat_cwm_er_splits$max_rooting_depth_m) # boundary fit singular 

SLA <- plot_dat(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
LDMC <- plot_dat(inv_nat_cwm_er_splits$`LDMC_g/g`)
leafN <- plot_dat(inv_nat_cwm_er_splits$`leafN_mg/g`)
leafP <- plot_dat(inv_nat_cwm_er_splits$`leafP_mg/g`)

RTD <- plot_dat(inv_nat_cwm_er_splits$Root_tissue_density)
Root_diam <- plot_dat(inv_nat_cwm_er_splits$Mean_Root_diameter)
RootN <- plot_dat(inv_nat_cwm_er_splits$Root_N_concentration)
SRL <- plot_dat(inv_nat_cwm_er_splits$Specific_root_length)

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

height_plot <- plot_function(inv_nat_cwm_er_splits$heightveg_m, height, "Height") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 60), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(0.1, 60), expand = c(0,0))

SLA_plot <- plot_function(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA, "SLA_mm2/mg")+ 
  scale_y_continuous(trans = "log2", limits = c(4, 60), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(4, 60), expand = c(0,0))

LDMC_plot <- plot_function(inv_nat_cwm_er_splits$`LDMC_g/g`, LDMC,  "LDMC_g/g") + 
  scale_y_continuous(limits = c(0.08, 0.5), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.08, 0.5), expand = c(0,0))

leafN_plot <- plot_function(inv_nat_cwm_er_splits$`leafN_mg/g`, leafN, "leafN_mg/g") + 
  scale_y_continuous(trans = "log2", limits = c(8, 57), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(8, 57), expand = c(0,0))

leafP_plot <- plot_function(inv_nat_cwm_er_splits$`leafP_mg/g`, leafP, "leafP_mg/g") + 
  scale_y_continuous(trans = "log2", limits = c(1, 8), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(1, 8), expand = c(0,0))

root_depth_plot <- plot_function(inv_nat_cwm_er_splits$max_rooting_depth_m, root_depth, "rooting_depth")+ 
  scale_y_continuous(trans = "log2", limits = c(0.1, 57), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(0.1, 57), expand = c(0,0))

RTD_plot <- plot_function(inv_nat_cwm_er_splits$Root_tissue_density, RTD, "RTD")  + 
  scale_y_continuous(limits = c(0.01, 0.75), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.01, 0.75), expand = c(0,0))

Root_diam_plot <- plot_function(inv_nat_cwm_er_splits$Mean_Root_diameter, Root_diam, "Root diameter")  + 
  scale_y_continuous(limits = c(0.1, 0.55), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.1, 0.55), expand = c(0,0))

RootN_plot <- plot_function(inv_nat_cwm_er_splits$Root_N_concentration, RootN, "RootN")  + 
  scale_y_continuous(limits = c(5, 25), expand = c(0,0)) + 
  scale_x_continuous( limits = c(5, 25), expand = c(0,0))

RootP_plot <- plot_function(inv_nat_cwm_er_splits$Root_P_concentration, RootP, "RootP")  + 
  scale_y_continuous(limits = c(1, 15), expand = c(0,0)) + 
  scale_x_continuous( limits = c(1, 15), expand = c(0,0))

SRL_plot <- plot_function(inv_nat_cwm_er_splits$Specific_root_length, SRL, "SRL")+ 
  scale_y_continuous(trans = "log2", limits = c(10, 660), expand = c(0,0)) + 
  scale_x_continuous(trans = "log2", limits = c(10, 660), expand = c(0,0))

legend <- get_legend(height_plot)

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_size_plots.pdf", height = 10, width = 12)
plot_grid(plot_grid(height_veg_densityplot, height_mod_plot, height_plot,
                    root_depth_densityplot, root_depth_mod_plot, root_depth_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_leaf_plots.pdf", height = 10, width = 12)
plot_grid(plot_grid(SLA_densityplot, SLA_mod_plot, SLA_plot,
                    LDMC_densityplot, LDMC_mod_plot, LDMC_plot,
                    leaf_N_densityplot, leaf_N_mod_plot, leafN_plot,
                    RTD_densityplot, RTD_mod_plot, RTD_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_leafnuts_plots.pdf", height = 10, width = 12)
plot_grid(plot_grid(leaf_N_densityplot, leaf_N_mod_plot, leafN_plot,
                    leaf_P_densityplot, leaf_P_mod_plot, leafP_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_rootcollab_plots.pdf", height = 10, width = 12)
plot_grid(plot_grid(SRL_densityplot, SRL_mod_plot, SRL_plot,
                    root_diam_densityplot, root_diam_mod_plot, Root_diam_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/sets_rootcons_plots.pdf", height = 10, width = 12)
plot_grid(plot_grid(root_N_densityplot, root_N_mod_plot, RootN_plot,
                    RTD_densityplot, RTD_mod_plot, RTD_plot,
                    ncol = 3, rel_widths = c(1,1,1)),
          plot_grid(density_legend, legend, ncol = 2, rel_widths = c(1,2)), ncol = 1, rel_heights = c(6,1))
dev.off()
