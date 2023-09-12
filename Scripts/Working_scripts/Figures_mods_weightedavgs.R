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
library(brms)

# ----------------------------------- load data -------------------------------------

CWM <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_CWM_1.csv", 
                           select = c("Plot", "Year", "TraitNameAbr", "CWM")))

dat <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wtraits_3.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean")))

plot_info <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage_3.csv", 
                                 select = c("Plot", "Zone", "Year", "Trait", "trait_cov", "N_relcov", "I_relcov", "NI_relcov", "UNKstat_relcov", "EcoRegionLevelI")))

SPCIS_ecoregion_IV <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Data/SPCIS_ecoreg_03182023.csv", 
                                          select = c("Plot", "US_L4CODE", "US_L4NAME")))

basal_area_plots <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Data/SPCIS_BA_plots_08102023.csv"))

colnames(CWM) <- c("Plot", "Year", "Trait", "CWM")

cb_palette <- c("#006633",  "#E69F00","#F0E442", "#D55E00", "#009E73", "#56B4E9", "black")

# ----------------------------------- data prep -------------------------------------

# drop duplicates - these ought to be dropped sooner 
plot_info  <- plot_info[!duplicated(plot_info),]

# drop basal area plots 
plot_info <- plot_info[!plot_info$Plot %in% c(basal_area_plots$Plot),]

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

# ----------------------------------- split by trait for plotting  -------------------------------------

CWM_80_splits <- split(CWM_plot_info_80_2, CWM_plot_info_80_2$Trait)

# --------------------------------- models - continuous responses -------------------------------------

mod_function <-function(dat){
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  trait <- dat$Trait[1]
  
  if(trait == "SLA_mm2/mg"| trait == "max_975_heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_975_rootdepth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    mod <- lmer(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
  }
  
  if(trait != "SLA_mm2/mg"& trait != "max_975_heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_975_rootdepth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    mod <- lmer(CWM ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI)+ (1|US_L4CODE), data = dat)
  }

  summary <- summary(mod)
  pvalue <- summary$coefficients[[10]]
  R2_marginal <- as.numeric(round(r.squaredGLMM(mod)[1],3))
  R2_conditional <- as.numeric(round(r.squaredGLMM(mod)[2],3))
  
  ls <- list(dat, R2_marginal, R2_conditional, mod, pvalue)
  return(ls)
}


max_height_mod <- mod_function(CWM_80_splits$max_975_heightveg_m)
root_maxdepth_mod <- mod_function(CWM_80_splits$max_975_rootdepth_m)

SLA_mod <- mod_function(CWM_80_splits$`SLA_mm2/mg`)
LDMC_mod <- mod_function(CWM_80_splits$`LDMC_g/g`)
leaf_N_mod <- mod_function(CWM_80_splits$`leafN_mg/g`)
leaf_P_mod <- mod_function(CWM_80_splits$`leafP_mg/g`)

root_diam_mod <- mod_function(CWM_80_splits$Mean_Root_diameter)
RTD_mod <- mod_function(CWM_80_splits$Root_tissue_density)
root_N_mod <- mod_function(CWM_80_splits$Root_N_concentration)
SRL_mod <- mod_function(CWM_80_splits$Specific_root_length)

# ----------------------------------- updated model plots (w ER level 4) -------------------------------------

mod_plot_function_updated <- function(ls, y1, y2, y3, title_lab){
  
  # data for plot 
  dat = ls[[1]]
  R2_marginal <- ls[[2]]
  R2_conditional <- ls[[3]]
  mod = ls[[4]]
  trait = dat$Trait[[1]]
  R2_M_label = paste("Marginal_R^2 ==",R2_marginal)
  R2_C_label = paste("Conditional_R^2 ==",R2_conditional)
  n_plot = nrow(dat)
  n_label = paste("n ==", n_plot)
  
  
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  
  # get EcoRegion1 lines by simulating data, setting random-effect of ER4 to average, and predicting lines 
  # prediction intervals 
  # get min and max for x-axis 
  
  xmin = min(dat$I_relcov_scaled)
  xmax = max(dat$I_relcov_scaled)
  
  # find average intercept for ecoregion 4 
  coefs <- coef(mod)
  avg_ER4 <- mean(coefs$US_L4CODE$`(Intercept)`)
  L4_coeffs <- as.data.frame(coefs$US_L4CODE)
  
  target.index <- which(abs(L4_coeffs$`(Intercept)` - avg_ER4) == min(abs(L4_coeffs$`(Intercept)`  - avg_ER4)))
  ER4_avg_value <- L4_coeffs$`(Intercept)`[target.index]
  ER4_avg_L4_code <- rownames(L4_coeffs)[L4_coeffs$`(Intercept)` == ER4_avg_value]
  
  new_dat <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                                EcoRegionLevelI = c(levels(as.factor(dat$EcoRegionLevelI))),
                                US_L4CODE = ER4_avg_L4_code)
  
  predict_int <- predictInterval(mod, new_dat, n.sims = 2000,
                                        stat = "median",
                                        type = "linear.prediction",
                                        include.resid.var = TRUE)
  
  out_all <- cbind(new_dat, predict_int)
  
  # exponentiate correct traits 
  if(trait == "SLA_mm2/mg"| trait == "max_975_heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_975_rootdepth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    out_all <- cbind(out_all[1:3],exp(out_all[4:6]))
  }
  
  if(trait != "SLA_mm2/mg"& trait != "max_975_heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_975_rootdepth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    out_all <- cbind(new_dat, predict_int)
    
  }

  # get overall trend line same as above  
  overall_intercept <- mean(coefs$EcoRegionLevelI$`(Intercept)`)
  overall_slope <- mean(coefs$EcoRegionLevelI$I_relcov_scaled)
  
  new_main_dat <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                                EcoRegionLevelI = "overall",
                                US_L4CODE = ER4_avg_L4_code)
  
  predict_int_main <- predictInterval(mod, new_main_dat, n.sims = 5000,
                                        stat = "median",
                                        type = "linear.prediction",
                                        include.resid.var = TRUE)
  
  main_all <- cbind(new_main_dat, predict_int_main)
  
  # exponentiate overall trend line for correct traits 
  if(trait == "SLA_mm2/mg"| trait == "max_975_heightveg_m" | trait == "leafP_mg/g" | 
     trait == "leafarea_mm2" | trait == "leafN_mg/g" | trait == "max_975_rootdepth_m"| 
     trait == "Mean_Root_diameter" | trait == "Root_N_concentation"| trait =="Specific_root_length"){
    main_all <- cbind(main_all[1:3],exp(main_all[4:6]))
  }
  
  if(trait != "SLA_mm2/mg"& trait != "max_975_heightveg_m" & trait != "leafP_mg/g" & 
     trait != "leafarea_mm2" & trait != "leafN_mg/g" & trait != "max_975_rootdepth_m"& 
     trait != "Mean_Root_diameter" & trait != "Root_N_concentation"& trait !="Specific_root_length"){
    main_all <- cbind(new_main_dat, predict_int_main)
  }
  
plot <-   ggplot() +
    geom_point(data = dat, aes(x = I_relcov_scaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.8, alpha = 0.4) +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS", ], aes(x = I_relcov_scaled, y = fit), color = "#006633") +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI =="GREAT PLAINS",], aes(x = I_relcov_scaled, y = fit), color = "#E69F00", linewidth = 1)  +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",],aes(x = I_relcov_scaled, y = fit), color = "#F0E442", linewidth = 1)  +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI =="NORTH AMERICAN DESERTS",],aes(x = I_relcov_scaled, y = fit), color = "#D55E00", linewidth = 1)  +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI =="NORTHERN FORESTS",], aes(x = I_relcov_scaled, y = fit), color = "#009E73", linewidth = 1)  +
    geom_smooth(data = out_all[out_all$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",], aes(x = I_relcov_scaled, y = fit), color = "#56B4E9", linewidth = 1) + 
    geom_smooth(data = main_all, aes(x = I_relcov_scaled, y = fit), color = "black", linewidth = 1) + 
    geom_ribbon(data = main_all, aes(x = I_relcov_scaled, ymin = lwr, ymax = upr), alpha = 0.2) + 
    labs(title = title_lab,
         x = "Scaled abundance of introduced species",
         y = title_lab) + 
    theme_bw() + 
    scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0))+
    scale_color_manual(values = cb_palette) + 
    scale_fill_manual(values = cb_palette) +
    theme(legend.position = "none", 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    guides(fill = guide_legend(nrow = 3, byrow = TRUE), 
           color = guide_legend(override.aes = list(alpha = 1, size = 4))) + 
  annotate(geom = "text", x = 0, y = y1, label = R2_M_label, parse = TRUE, hjust = -0.08, size = 4) + 
  annotate(geom = "text", x = 0, y = y2, label = R2_C_label, parse = TRUE, hjust = -0.08, size = 4) + 
  annotate(geom = "text", x = 0, y = y3, label = n_label, parse = TRUE, hjust = -0.2, size = 4) 
return(plot)
}

legend <- get_legend(maxheight_mod_plot)

maxheight_mod_plot <- mod_plot_function_updated(max_height_mod, 0.22, 0.15, 0.1, "CWM max height (m)") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 46)) 

maxroot_mod_plot <- mod_plot_function_updated(root_maxdepth_mod, 0.22, 0.15, 0.1, "CWM max rooting depth (m)") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 46)) 

SLA_mod_plot <- mod_plot_function_updated(SLA_mod, 5, 4, 3.1,  expression(paste("CWM SLA (m" ^{2}," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(3, 60)) 

LDMC_mod_plot <- mod_plot_function_updated(LDMC_mod, 0.17, 0.14, 0.105, expression(paste("CWM LDMC (g" ," g" ^{-1},")"))) + 
  ylim(0.1, 0.52) 

LeafN_mod_plot <- mod_plot_function_updated(leaf_N_mod, 11.1,9.6,8.3, expression(paste("CWM Leaf N (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(8, 50)) 

LeafP_mod_plot <- mod_plot_function_updated(leaf_P_mod, 0.72,0.6,0.5, expression(paste("CWM Leaf P (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(0.5, 5)) 

RootN_mod_plot <- mod_plot_function_updated(root_N_mod, 6,5.25,4.6, expression(paste("CWM Root N (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(4.5, 32))

RTD_mod_plot <- mod_plot_function_updated(RTD_mod, 0.1, 0.066, 0.028, expression(paste("CWM RTD (g" ," cm" ^{-3},")"))) + 
  ylim(0.025, 0.52) 

RootDiam_mod_plot <- mod_plot_function_updated(root_diam_mod, 0.133,0.115,0.101, "CWM root diameter (mm)") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 0.8))

SRL_mod_plot <- mod_plot_function_updated(SLA_mod, 4.25, 3.25, 2.51,   expression(paste("CWM SRL (m" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(2.5, 130)) 

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/mods_LES_plots.pdf", height = 8.35, width = 8.35)
plot_grid(plot_grid(SLA_mod_plot, LDMC_mod_plot, LeafN_mod_plot, LeafP_mod_plot,
                    ncol = 2, rel_widths = c(1,1,1,1), labels = c("a", "b", "c", "d")), legend, ncol = 1, rel_heights = c(6,0.8))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/mods_RootCons_plots.pdf", height = 4.175, width = 8.35)
plot_grid(plot_grid(RootN_mod_plot,RTD_mod_plot,
                    ncol = 2, rel_widths = c(1,1), labels = c("a", "b")), legend, ncol = 1, rel_heights = c(3,0.4))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/mods_RootCollab_plots.pdf", height = 4.175, width = 8.35)
plot_grid(plot_grid(SRL_mod_plot,RootDiam_mod_plot,
                    ncol = 2, rel_widths = c(1,1), labels = c("a", "b"), align = "v"), legend, ncol = 1, rel_heights = c(3,0.4))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/mods_size_plots.pdf", height = 4.175, width = 8.35)
plot_grid(plot_grid(maxheight_mod_plot,maxroot_mod_plot,
                    ncol = 2, rel_widths = c(1,1), labels = c("a", "b"), align = "v"), legend, ncol = 1, rel_heights = c(3,0.4))
dev.off()

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

# need averages for the different plots  
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
    mod_out <- lmer(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), data = dat_long)
  }
  
  if(TraitNameAbr != "SLA_mm2/mg"& TraitNameAbr != "heightveg_m" & TraitNameAbr != "leafP_mg/g" & 
     TraitNameAbr != "leafarea_mm2" & TraitNameAbr != "leafN_mg/g" & TraitNameAbr != "max_rooting_depth_m"& 
     TraitNameAbr != "Mean_Root_diameter" & TraitNameAbr != "Root_N_concentation"& TraitNameAbr !="Specific_root_length"){
    mod_out <- lmer(CWM ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), data = dat_long)
  }
  
  sum_mod <- summary(mod_out)
  return(sum_mod)
} 

weighted_avg_mod(inv_nat_cwm_er_splits$max_975_heightveg_m)
weighted_avg_mod(inv_nat_cwm_er_splits$max_975_rootdepth_m)

weighted_avg_mod(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
weighted_avg_mod(inv_nat_cwm_er_splits$`LDMC_g/g`)
weighted_avg_mod(inv_nat_cwm_er_splits$`leafN_mg/g`)
weighted_avg_mod(inv_nat_cwm_er_splits$`leafP_mg/g`)

weighted_avg_mod(inv_nat_cwm_er_splits$Specific_root_length)
weighted_avg_mod(inv_nat_cwm_er_splits$Mean_Root_diameter)
weighted_avg_mod(inv_nat_cwm_er_splits$Root_N_concentration)
weighted_avg_mod(inv_nat_cwm_er_splits$Root_tissue_density)

 
# get data for plots 
height <- plot_dat_updated(inv_nat_cwm_er_splits$max_975_heightveg_m)
root_depth <- plot_dat_updated(inv_nat_cwm_er_splits$max_975_rootdepth_m) 

SLA <- plot_dat_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
LDMC <- plot_dat_updated(inv_nat_cwm_er_splits$`LDMC_g/g`)
leafN <- plot_dat_updated(inv_nat_cwm_er_splits$`leafN_mg/g`)
leafP <- plot_dat_updated(inv_nat_cwm_er_splits$`leafP_mg/g`)

RTD <- plot_dat_updated(inv_nat_cwm_er_splits$Root_tissue_density)
Root_diam <- plot_dat_updated(inv_nat_cwm_er_splits$Mean_Root_diameter)
RootN <- plot_dat_updated(inv_nat_cwm_er_splits$Root_N_concentration)
SRL <- plot_dat_updated(inv_nat_cwm_er_splits$Specific_root_length)

# plotting 
plot_function_updated <- function(df, means_df, trait, title_lab, x_position, y_max){
  
  # labels for plots 
  n_plot = nrow(df)
  n_label = paste("n ==", n_plot)
  
  
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
    labs(title = title_lab) + 
    ylab(paste("Introduced species \n Community weighted average", trait)) + 
    xlab(paste("Native species \n Community weighted average", trait)) + 
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    annotate(geom = "text", x = x_position, y = y_max - (y_max * 0.1) , label = n_label, parse = TRUE, hjust = 0, size = 4) 
  return(plot)
}

height_plot <- plot_function_updated(inv_nat_cwm_er_splits$max_975_heightveg_m, height,
                                     "max height (m)", "Max height (m)", 0, 75) + 
  scale_y_continuous( limits = c(-5, 75), expand = c(0,0)) + 
  scale_x_continuous( limits = c(-5, 75), expand = c(0,0))

rootdepth_plot <- plot_function_updated(inv_nat_cwm_er_splits$max_975_rootdepth_m, root_depth,
                                        "max rooting depth (m)", "Max rooting depth (m)", 0, 30) + 
  scale_y_continuous( limits = c(-2, 30), expand = c(0,0)) + 
  scale_x_continuous( limits = c(-2, 30), expand = c(0,0))

SLA_plot <- plot_function_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA,
                                 "SLA", expression(paste("SLA (m" ^{2}," g" ^{-1},")")), 6, 65) + 
  scale_y_continuous( limits = c(4, 65), expand = c(0,0)) + 
  scale_x_continuous( limits = c(4, 65), expand = c(0,0))

LDMC_plot <- plot_function_updated(inv_nat_cwm_er_splits$`LDMC_g/g`, LDMC,
                                  "LDMC", expression(paste("LDMC (g" ," g" ^{-1},")")), 0.1, 0.55) + 
  scale_y_continuous( limits = c(0.08, 0.55), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.08, 0.55), expand = c(0,0))


leafN_plot <- plot_function_updated(inv_nat_cwm_er_splits$`leafN_mg/g`, leafN, 
                                    "leaf N", expression(paste("Leaf N (mg" ," g" ^{-1},")")), 10, 50) + 
  scale_y_continuous(limits = c(8, 50), expand = c(0,0)) + 
  scale_x_continuous(limits = c(8, 50), expand = c(0,0))

leafP_plot <- plot_function_updated(inv_nat_cwm_er_splits$`leafP_mg/g`, leafP, 
                                    "leaf P", expression(paste("Leaf P (mg" ," g" ^{-1},")")), 1.25, 8) + 
  scale_y_continuous(limits = c(1, 8), expand = c(0,0)) + 
  scale_x_continuous(limits = c(1,8), expand = c(0,0))

RTD_plot <- plot_function_updated(inv_nat_cwm_er_splits$Root_tissue_density, RTD, 
                                  "RTD",expression(paste("RTD (g" ," cm" ^{-3},")")), 0.05, 0.75)  + 
  scale_y_continuous(limits = c(0.01, 0.75), expand = c(0,0)) + 
  scale_x_continuous( limits = c(0.01, 0.75), expand = c(0,0))

Root_diam_plot <- plot_function_updated(inv_nat_cwm_er_splits$Mean_Root_diameter, Root_diam, 
                                        "root diameter", "Root diameter (mm)", -0.2, 3)  + 
  scale_y_continuous(limits = c(-0.25, 3), expand = c(0,0)) + 
  scale_x_continuous( limits = c(-0.25, 3), expand = c(0,0))

RootN_plot <- plot_function_updated(inv_nat_cwm_er_splits$Root_N_concentration, RootN, 
                                    "Root N", expression(paste("Root N (mg" ," g" ^{-1},")")), 5, 35) + 
  scale_y_continuous(limits = c(3, 35), expand = c(0,0)) + 
  scale_x_continuous( limits = c(3, 35), expand = c(0,0))

SRL_plot <- plot_function_updated(inv_nat_cwm_er_splits$Specific_root_length, SRL, 
                                  "SRL", expression(paste("SRL (m" ," g" ^{-1},")")), 10, 660)+ 
  scale_y_continuous(limits = c(5, 660), expand = c(0,0)) + 
  scale_x_continuous(limits = c(5, 660), expand = c(0,0))

legend_2 <- legend + theme(plot.background = element_rect(color = "black", linewidth = 1))

LES_title <- ggdraw() +
  draw_label("Leaf conservation gradient")

RootCons_title <- ggdraw() +
  draw_label("Root conservation gradient")

RootCollab_title <- ggdraw() +
  draw_label("Root collaboration gradient")

Size_title <- ggdraw() +
  draw_label("Size gradients")


LES_grid <- plot_grid(LES_title, plot_grid(SLA_plot, LDMC_plot, leafN_plot, leafP_plot, ncol = 2, rel_widths = c(1,1,1,1), labels = c("a", "b", "c", "d"), scale = 0.95),
                      ncol = 1, rel_heights = c(1,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCons_grid <- plot_grid(RootCons_title, plot_grid(RootN_plot, RTD_plot, ncol = 1, rel_heights = c(1,1), labels = c("e", "f"), scale = 0.95),
                           ncol = 1, rel_heights = c(1,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCollab_grid <- plot_grid(RootCollab_title, plot_grid(SRL_plot, Root_diam_plot, ncol = 1, rel_heights = c(1,1), labels = c("a", "b"), scale = 0.95),
                             ncol = 1, rel_heights = c(0.5,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

Size_grid <- plot_grid(Size_title, plot_grid(height_plot, rootdepth_plot, ncol = 1, rel_heights = c(1,1), labels = c("c", "d"), scale = 0.95),
                             ncol = 1, rel_heights = c(0.5,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/within_community_plots_conservationgradients_2.pdf", height = 10, width = 13)
plot_grid(plot_grid(LES_grid, RootCons_grid, nrow = 1, rel_widths = c(2,1)), legend, ncol = 1, rel_heights = c(10,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/within_community_plots_rootcollab_size.pdf", height = 10, width = 10)
plot_grid(plot_grid(RootCollab_grid, Size_grid, nrow = 1, rel_widths = c(1,1)), legend, ncol = 1, rel_heights = c(10,1))
dev.off()

#
#
#
#
#
#
#
#
# plotting 
plot_function_log_scale_updated <- function(df, means_df, trait, title_lab, x_position, y_max){
  
  # labels for plots 
  n_plot = nrow(df)
  n_label = paste("n ==", n_plot)
  
  
  plot <-   ggplot(df) + 
    geom_point(aes(x = CWM_N, y = CWM_I, color = EcoRegionLevelI), size = 1, alpha = 0.3) + 
    geom_point(data = means_df,
               aes(x = weighted_mean.N, y = weighted_mean.I, fill = EcoRegionLevelI),
               shape = 21, size = 5, alpha = 1) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9"))+
    scale_fill_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9"))+
    labs(title = title_lab) + 
    ylab(paste("Introduced species \n Community weighted average", trait)) + 
    xlab(paste("Native species \n Community weighted average", trait)) + 
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    annotate(geom = "text", x = x_position, y = y_max - (y_max * 0.1) , label = n_label, parse = TRUE, hjust = 0, size = 4) 
  return(plot)
}

SLA_log_plot <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA,
                                                "SLA", expression(paste("SLA (m" ^{2}," g" ^{-1},")")), 6, 65) + 
  scale_y_continuous(trans = "log2", limits = c(4, 65))+ 
  scale_x_continuous(trans = "log2", limits = c(4, 65))

SLA_log_plot_werror <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA,
                                                "SLA", expression(paste("SLA (m" ^{2}," g" ^{-1},")")), 6, 65) + 
  scale_y_continuous(trans = "log2", limits = c(4, 65))+ 
  scale_x_continuous(trans = "log2", limits = c(4, 65)) +
  geom_errorbar(data = SLA, aes(x = weighted_mean.N,
                                     ymin = weighted_mean.I - sd.I, 
                                     ymax = weighted_mean.I + sd.I), width = 0) + 
  geom_errorbarh(data = SLA, aes( y = weighted_mean.I, 
                                       xmin = weighted_mean.N - sd.N, 
                                       xmax = weighted_mean.N + sd.N), height = 0) 

SRL_log_plot <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$`Specific_root_length`, SRL,
                                                "SRL", expression(paste("SRL (m" ^{2}," g" ^{-1},")")), 6, 65) + 
  scale_y_continuous(trans = "log2", limits = c(4, 665))+ 
  scale_x_continuous(trans = "log2", limits = c(4, 665))

SRL_log_plot_werror <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$`Specific_root_length`, SRL,
                                                       "SRL", expression(paste("SRL (m" ^{2}," g" ^{-1},")")), 6, 65) + 
  scale_y_continuous(trans = "log2", limits = c(4, 665))+ 
  scale_x_continuous(trans = "log2", limits = c(4, 665)) +
  geom_errorbar(data = SRL, aes(x = weighted_mean.N,
                                ymin = weighted_mean.I - sd.I, 
                                ymax = weighted_mean.I + sd.I), width = 0) + 
  geom_errorbarh(data = SRL, aes( y = weighted_mean.I, 
                                  xmin = weighted_mean.N - sd.N, 
                                  xmax = weighted_mean.N + sd.N), height = 0) 

Height_log_plot <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$max_975_heightveg_m, height,
                                                "Height", expression(paste("Height (m))")), 6, 65) + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 85))+ 
  scale_x_continuous(trans = "log2", limits = c(0.1, 85))

Height_log_plot_werror <- plot_function_log_scale_updated(inv_nat_cwm_er_splits$max_975_heightveg_m, height,
                                                          "Height", expression(paste("Height (m)")), 6, 85) + 
  scale_y_continuous(trans = "log2", limits = c(0.01, 85))+ 
  scale_x_continuous(trans = "log2", limits = c(0.01, 85)) +
  geom_errorbar(data = height, aes(x = weighted_mean.N,
                                ymin = weighted_mean.I - sd.I, 
                                ymax = weighted_mean.I + sd.I), width = 0) + 
  geom_errorbarh(data = height, aes( y = weighted_mean.I, 
                                  xmin = weighted_mean.N - sd.N, 
                                  xmax = weighted_mean.N + sd.N), height = 0) 


# -------- difference in weighted mean of natives vs. introduced species 

weighted_avg_mod <- function(df){
  df$difference = df$CWM_I - df$CWM_N
  TraitNameAbr <- df$TraitNameAbr[1]
  mod_out <- lmer(difference ~ (1|EcoRegionLevelI/US_L4NAME), data = df)
  sum_mod <- summary(mod_out)
  return(sum_mod)
} 

out <- weighted_avg_mod (inv_nat_cwm_er_splits$max_975_heightveg_m)


# differences plot 

df = inv_nat_cwm_er_splits$max_975_heightveg_m
df$difference = df$CWM_I - df$CWM_N

plot <- ggplot(df) + 
  geom_point(aes(y = difference, x = EcoRegionLevelI, color = EcoRegionLevelI), size = 1, alpha = 0.3) +
  scale_color_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9")) +
  scale_fill_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9")) +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.text = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10), 
        legend.spacing.y = unit(0.15, "cm"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + coord_flip()
  
# now plot the model estimate onto this figure somehow...? 

box_plot_function <- function(df){
  df$difference = df$CWM_I - df$CWM_N
  box_plot <- ggplot(df, aes(y = difference, x = 0)) + geom_violin() +   
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + coord_flip()
  return(box_plot)
}

box_plot_function(inv_nat_cwm_er_splits$`SLA_mm2/mg`) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)





# ---------------- what % of the dataset for the traits we use is at the genus level -----------
dat_ss <- dat[dat$TraitNameAbr %in% c("LDMC_g/g", "leafN_mg/g", "leafP_mg/g", "max_975_heightveg_m",
                                      "max_975_rootdepth_m", "Mean_Root_diameter", "Root_N_concentration",
                                      "Root_tissue_density", "SLA_mm2/mg", "Specific_root_length"),]


# if space is in the "AcceptedTaxonName" then new column gets YES otherwise NO 

dat_ss$genus_only <- ifelse(grepl(" ", dat_ss$AcceptedTaxonName), "NO", "YES")
percent_rel_genus <- sum(dat_ss[dat_ss$genus_only == "YES",]$rel_100)/sum(dat_ss$rel_100)


# testing brms for overall models
dat_test <- CWM_80_splits$max_975_heightveg_m 
dat_test <- dat_test[dat_test$Year > 2010,]

dat_test <- dat_test[!is.na(dat_test$EcoRegionLevelI),]
dat_test$I_relcov_scaled <- scale(dat_test$I_relcov, center = TRUE, scale = TRUE)
brms_mod <- brm(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), 
                data = dat_test,
                chains = 2, cores = 4,
                warmup = 250, iter = 1000)

View(brms_mod$fit)
summary(brms_mod)
hypothesis(brms_mod, "Intercept = 0", alpha = 0.1)
posterior_samples(brms_mod)

# testing brms for within community analyses 
dat_com <- inv_nat_cwm_er_splits$`SLA_mm2/mg`
dat_com_height <- inv_nat_cwm_er_splits$max_975_heightveg_m
dat_com_LDMC <- inv_nat_cwm_er_splits$`LDMC_g/g`

dat_com_long <- reshape(dat_com , 
                    varying = c("CWM_I", "CWM_N"),
                    v.names = "CWM",
                    idvar = c("Plot"),
                    timevar = "NativeStatus",
                    times = c("I", "N"),
                    direction = "long")

dat_com_long_height <- reshape(dat_com_height , 
                        varying = c("CWM_I", "CWM_N"),
                        v.names = "CWM",
                        idvar = c("Plot"),
                        timevar = "NativeStatus",
                        times = c("I", "N"),
                        direction = "long")

dat_com_long_LDMC <- reshape(dat_com_LDMC , 
                        varying = c("CWM_I", "CWM_N"),
                        v.names = "CWM",
                        idvar = c("Plot"),
                        timevar = "NativeStatus",
                        times = c("I", "N"),
                        direction = "long")

# models
mod_com_out <- lmer(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), data = dat_com_long)

# brms models 
brms_mod_com_out <- brm(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), 
                        data = dat_com_long,
                        chains = 2, cores = 4,
                        warmup = 2000, iter = 3000, control = list(adapt_delta = 0.95))

brms_mod_com_out_nested <- brm(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME/Plot), 
                        data = dat_com_long,
                        chains = 2, cores = 4,
                        warmup = 2000, iter = 3000, control = list(adapt_delta = 0.95))

brms_mod_com_out_height <- brm(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), 
                        data = dat_com_long_height,
                        chains = 4, cores = 4,
                        warmup = 2000, iter = 3000, control = list(adapt_delta = 0.95))

brms_mod_com_out_LDMC <- brm(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), 
                               data = dat_com_long_LDMC,
                               chains = 4, cores = 4,
                             warmup = 2000, iter = 3000, control = list(adapt_delta = 0.95))

#### is the issue arising because there is only one plot per ecoregion level 4?

dat_com_SLA <- inv_nat_cwm_er_splits$`SLA_mm2/mg`

# split by ecoregion level 4 

ER4_splits <- split(dat_com_SLA, dat_com_SLA$US_L4NAME)

# keep ER4 with > 5 rows 
# this does not slove the problem 
ER4_2 <- do.call(rbind, ER4_splits[sapply(ER4_splits, function(x) nrow(x) > 5)])

dat_com_long <- reshape(ER4_2 , 
                        varying = c("CWM_I", "CWM_N"),
                        v.names = "CWM",
                        idvar = c("Plot"),
                        timevar = "NativeStatus",
                        times = c("I", "N"),
                        direction = "long")

# models
mod_com_out <- lmer(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME/Plot), data = dat_com_long)
