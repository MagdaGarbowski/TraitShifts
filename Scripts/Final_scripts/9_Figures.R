
# project: Trait-Shifts 
# objective: analyses & figures 
# author: Magda Garbowski 
# date: November 22, 2023

library(lme4)
library(lmerTest)
library(MuMIn)
library(merTools)
library(cowplot)
library(gridGraphics)
library(ggplot2)

options(scipen = 999)

# ------------------------------------------------------- DATA -----------------------------------------------------

CWM_plot_info_80 <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80.csv")
CWM_plot_info_80_2per <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80_2percent.csv")
inv_nat_cwm_er <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_cooccurring_weighted_means.csv")

cb_palette <- c("#009E73","#CC79A7","#E69F00", "#D55E00", "#56B4E9", "#006633", "black")

# ----------------------------------------------- FUNCTIONS -------------------------------------------------------
# ------------------------------------------------ models ---------------------------------------------------------

# model function - abundance gradients 
mod_function <-function(dat){
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  dat$unscaled <- dat$I_relcov_scaled * attr(dat$I_relcov_scaled, "scaled:scale") + attr(dat$I_relcov_scaled, "scaled:center")
  trait <- dat$Trait[1]
  
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    mod <- lmer(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    mod <- lmer(CWM ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI)+ (1|US_L4CODE), data = dat)
  }
  
  summary_mod <- summary(mod)
  coeffs_out <- coef(mod)
  pvalue <- summary_mod$coefficients[[10]]
  R2_marginal <- as.numeric(round(r.squaredGLMM(mod)[1],3))
  R2_conditional <- as.numeric(round(r.squaredGLMM(mod)[2],3))
  
  ls <- list(dat, R2_marginal, R2_conditional, mod, pvalue, summary_mod, coeffs_out)
  return(ls)
}

# model function - co-occurring native and introduced species trait differences
diff_mod_function <- function(df){
  trait <- df$Trait[1]
  
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    df$difference = log(df$CWM_I) - log(df$CWM_N)
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    df$difference = df$CWM_I - df$CWM_N
  }
  
  mod_out <- lmer(difference ~ (1|EcoRegionLevelI/US_L4NAME), data = df)
  sum_mod <- summary(mod_out)
  
  # get estimates on random effects 
  out_RE <- as.data.frame(REsim(mod_out, n.sims = 200))
  out_RE_ss <- out_RE[out_RE$groupFctr == "EcoRegionLevelI",]
  
  return(list(sum_mod, out_RE_ss))
}

# ----------------------------------------------- FUNCTIONS -------------------------------------------------------
# ------------------------------------------------ figures --------------------------------------------------------

# plotting function - abundance gradients 
mod_plot_function <- function(ls, y1, y2, y3, y4, title_lab){
  # data for plot 
  dat = ls[[1]]
  R2_marginal <- ls[[2]]
  R2_conditional <- ls[[3]]
  mod = ls[[4]]
  p = ls[[5]]
  p_val = paste(round(ls[[5]],4))
  trait = dat$Trait[[1]]
  
  scale_att_1 <- attr(dat$I_relcov_scaled, "scaled:scale")
  scale_center <- attr(dat$I_relcov_scaled, "scaled:center")
  dat$unscaled <- dat$I_relcov_scaled * attr(dat$I_relcov_scaled, "scaled:scale") + attr(dat$I_relcov_scaled, "scaled:center")
  xmin <- min(dat$I_relcov_scaled)
  xmax <- max(dat$I_relcov_scaled)
  
  R2_M_label = paste("Marginal_R^2 ==",R2_marginal)
  R2_C_label = paste("Conditional_R^2 ==",R2_conditional)
  n_plot = nrow(dat)
  n_label = paste("n ==", n_plot)
  p_label = ifelse(p < 0.0001, paste("p < 0.0001"), paste("p ==", p_val))
  
  # simulate data, set random-effect of ER4 to mean, get predictions 
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
                                 stat = "mean",
                                 type = "linear.prediction",
                                 include.resid.var = FALSE)
  
  out_all <- cbind(new_dat, predict_int)
  out_all$unscaled <- out_all$I_relcov_scaled * scale_att_1 + scale_center
  out_all$unscaled <- round(out_all$unscaled,6)
  
  # exponentiate traits 
  
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    out_all <- cbind(out_all[,c(1,2,3,7)],exp(out_all[4:6]))
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    out_all <- cbind(out_all[,c(1,2,3,7)],out_all[4:6])
  }
  
  # get overall trend line same as above  
  overall_intercept <- mean(coefs$EcoRegionLevelI$`(Intercept)`)
  overall_slope <- mean(coefs$EcoRegionLevelI$I_relcov_scaled)
  
  new_main_dat <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                              EcoRegionLevelI = "overall",
                              US_L4CODE = ER4_avg_L4_code)
  
  predict_int_main <- predictInterval(mod, new_main_dat, n.sims = 2000,
                                      stat = "mean",
                                      type = "linear.prediction",
                                      include.resid.var = FALSE)
  
  main_all <- cbind(new_main_dat, predict_int_main)
  main_all$unscaled <- main_all$I_relcov_scaled * scale_att_1 + scale_center
  main_all$unscaled <- round(main_all$unscaled, 6)
  
  # exponentiate overall trend line for correct traits 
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    main_all <- cbind(main_all[,c(1,2,3,7)], exp(main_all[4:6]))
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    main_all <- cbind(main_all[,c(1,2,3,7)], main_all[4:6])
  }
  
  # get min and max for x-axis 
  xmin_unscaled = min(dat$unscaled)
  xmax_unscaled = max(dat$unscaled)
  
  plot_2 <- ggplot() + 
    geom_point(data = dat, aes(x = unscaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.8, alpha = 0.2) + 
    geom_point(data = out_all[out_all$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS", ], aes(x = unscaled, y = fit), color = "#009E73", size = 0.5) + 
    geom_point(data = out_all[out_all$EcoRegionLevelI =="GREAT PLAINS",], aes(x = unscaled, y = fit), color = "#CC79A7",  size = 0.5)  +
    geom_point(data = out_all[out_all$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",],aes(x = unscaled, y = fit), color = "#E69F00",  size = 0.5)  +
    geom_point(data = out_all[out_all$EcoRegionLevelI =="NORTH AMERICAN DESERTS",],aes(x = unscaled, y = fit), color = "#D55E00", size = 0.5)  +
    geom_point(data = out_all[out_all$EcoRegionLevelI =="NORTHERN FORESTS",], aes(x = unscaled, y = fit), color = "#56B4E9",  size = 0.5)  +
    geom_point(data = out_all[out_all$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",], aes(x = unscaled, y = fit), color = "#006633",  size = 0.5) + 
    geom_point(data = main_all, aes(x = unscaled, y = fit), color = "black", size = 1) + 
    geom_ribbon(data = main_all, aes(x = unscaled, ymin = lwr, ymax = upr), alpha = 0.4) + 
    theme_bw() + 
    scale_x_continuous(limits = c(0, 100), expand = c(0,0))+
    scale_color_manual(values = cb_palette) + 
    scale_fill_manual(values = cb_palette) +
    labs(title = title_lab,
         x = "Relative abundance of introduced species",
         y = title_lab) +
    theme(legend.position = "none", 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    guides(fill = guide_legend(nrow = 3, byrow = TRUE), 
           color = guide_legend(override.aes = list(alpha = 1, size = 4))) + 
    annotate(geom = "text", x = 25, y = y1, label = R2_M_label, parse = TRUE, hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = 25, y = y2, label = R2_C_label, parse = TRUE, hjust = 0, size = 3.5) +
    annotate(geom = "text", x = 25, y = y3, label = p_label, parse = TRUE, hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = 25, y = y4, label = n_label, parse = TRUE, hjust = 0, size = 3.5) 
  return(plot_2)
}

# plotting function - co-occurring native and introduced trait differences
diff_mod_function_plot <- function(df, ann_y_1, ann_y_2){
  df_2 <- df
  df_2$EcoRegionLevelI <- "OVERALL"
  df_3 <- rbind(df, df_2)
  
  trait <- df_3$Trait[1]
  
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    df_3$difference = log(df_3$CWM_I) - log(df_3$CWM_N)
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    df_3$difference = df_3$CWM_I - df_3$CWM_N
  }
  
  mod_out <- lmer(difference ~ (1|EcoRegionLevelI/US_L4NAME), data = df_3[!df_3$EcoRegionLevelI == "OVERALL",])
  sum_mod <- summary(mod_out)
  p_val = round(sum_mod$coefficients[5],4)
  n = nrow(df)
  out_RE <- as.data.frame(REsim(mod_out, n.sims = 200))
  out_RE_ss <- out_RE[out_RE$groupFctr == "EcoRegionLevelI",]
  out_FE <- as.data.frame(FEsim(mod_out, n.sims = 200))
  out_FE$groupFctr <- "OVERALL"
  out_FE$groupID <- "OVERALL"
  out_effs_all <- rbind(out_FE, out_RE_ss)
  out_effs_all$mean_eff <- ifelse(out_effs_all$groupID == "OVERALL", out_effs_all$mean,
                                  out_effs_all$mean + (out_effs_all[out_effs_all$groupID == "OVERALL",]$mean)) 
  
  plot_out <-  ggplot(out_effs_all, aes(x = groupID, y = mean_eff)) + 
    geom_point(data = df_3, aes(x = EcoRegionLevelI, y = difference, color = EcoRegionLevelI), alpha = 0.1, position = position_jitter(width = 0.25)) + 
    geom_linerange(aes(ymin = mean_eff - (1.96 * sd), 
                       ymax = mean_eff + (1.96 * sd)), color = "gray5", linewidth = 0.6) + 
    geom_point(aes(x = groupID, y = mean_eff), color = "gray5", shape = 16, size = 4) +
    scale_color_manual(values = c("#009E73","#CC79A7","#E69F00", "#D55E00", "#56B4E9", "#006633", "gray50")) +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "gray5", linetype = "dashed") +
    theme_bw() + 
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 11),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 11),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 11), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank()) + 
    guides(fill = "none", shape = "none", color = guide_legend(override.aes = list(alpha = 1, size = 3), ncol = 2)) + 
    annotate(geom = "text", x = "NORTHERN FORESTS", y = ann_y_1, label = paste ("n =", n), hjust = 0, size = 4) + 
    annotate(geom = "text", x = "NORTHERN FORESTS", y = ann_y_2, label = paste ("p =", p_val), hjust = 0, size = 4) 
  
  return(plot_out)
}  

# ------------------------------------------------ trait splits ----------------------------------------------------

CWM_80_splits <- split(CWM_plot_info_80, CWM_plot_info_80$Trait)
CWM_80_splits_2per <- split(CWM_plot_info_80_2per, CWM_plot_info_80_2per$Trait)
inv_nat_cwm_er_splits <- split(inv_nat_cwm_er, inv_nat_cwm_er$Trait)

# ----------------------------------------------- ANALYSES -------------------------------------------------------
# ----------------------------- abundance gradient output - full dataset ------------------------------------------

SLA_mod <- mod_function(CWM_80_splits$`Specific leaf area`)
LDMC_mod <- mod_function(CWM_80_splits$`Leaf dry matter content`)
leaf_N_mod <- mod_function(CWM_80_splits$`Leaf nitrogen concentration`)
leaf_P_mod <- mod_function(CWM_80_splits$`Leaf phosphorus concentration`)
root_diam_mod <- mod_function(CWM_80_splits$`Root diameter`)
RTD_mod <- mod_function(CWM_80_splits$`Root tissue density`) 
root_N_mod <- mod_function(CWM_80_splits$`Root nitrogen concentration`)
SRL_mod <- mod_function(CWM_80_splits$`Specific root length`)
max_height_mod <- mod_function(CWM_80_splits$`Maximum height`)
root_maxdepth_mod <- mod_function(CWM_80_splits$`Maximum rooting depth`)

# ------------------------- abundance gradient output - > 2% introduced dataset ------------------------------------------

SLA_mod_2per <- mod_function(CWM_80_splits_2per$`Specific leaf area`)
LDMC_mod_2per <- mod_function(CWM_80_splits_2per$`Leaf dry matter content`)
leaf_N_mod_2per <- mod_function(CWM_80_splits_2per$`Leaf nitrogen concentration`)
leaf_P_mod_2per <- mod_function(CWM_80_splits_2per$`Leaf phosphorus concentration`)
root_diam_mod_2per <- mod_function(CWM_80_splits_2per$`Root diameter`)
RTD_mod_2per <- mod_function(CWM_80_splits_2per$`Root tissue density`) 
root_N_mod_2per <- mod_function(CWM_80_splits_2per$`Root nitrogen concentration`)
SRL_mod_2per <- mod_function(CWM_80_splits_2per$`Specific root length`)
max_height_mod_2per <- mod_function(CWM_80_splits_2per$`Maximum height`)
root_maxdepth_mod_2per <- mod_function(CWM_80_splits_2per$`Maximum rooting depth`)

# ----------------------------------- co-occurrence output  -----------------------------------------

SLA_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Specific leaf area`)
LDMC_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Leaf dry matter content`)
LeafN_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Leaf nitrogen concentration`)
LeafP_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Leaf phosphorus concentration`)
RootN_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Root nitrogen concentration`)
RTD_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Root tissue density`)
SRL_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Specific root length`)
RootDiam_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Root diameter`)
height_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Maximum height`)
rootdepth_cooccurrence_mod <- diff_mod_function(inv_nat_cwm_er_splits$`Maximum rooting depth`)

# ------------------------------------ abundance gradient plots  -------------------------------------------------------

maxheight_mod_plot <- mod_plot_function(max_height_mod, 0.32, 0.225, 0.15, 0.11, expression(paste("CWM max height (m"^{},")"))) +
  scale_y_continuous(trans = "log2", limits = c(0.11, 46))

maxroot_mod_plot <- mod_plot_function(root_maxdepth_mod, 0.32, 0.225, 0.15, 0.11, "CWM max rooting depth (m)") + 
  scale_y_continuous(trans = "log2", limits = c(0.11, 46)) 

SLA_mod_plot <- mod_plot_function(SLA_mod, 6, 4.8, 3.76, 3.1,  expression(paste("CWM SLA (m" ^{2}," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(3, 60)) 

LDMC_mod_plot <- mod_plot_function(LDMC_mod, 0.195, 0.165, 0.133, 0.105, expression(paste("CWM LDMC (g" ," g" ^{-1},")"))) + 
  ylim(0.1, 0.52) 

LeafN_mod_plot <- mod_plot_function(leaf_N_mod, 12,10.5,9,8, expression(paste("CWM Leaf N (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(8, 50)) 

LeafP_mod_plot <- mod_plot_function(leaf_P_mod, 0.82, 0.695,0.575,0.5, expression(paste("CWM Leaf P (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(0.5, 5)) 

RootN_mod_plot <- mod_plot_function(root_N_mod, 6.9,6,5.175,4.6, expression(paste("CWM Root N (mg" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(4.5, 32))

RTD_mod_plot <- mod_plot_function(RTD_mod, 0.135, 0.1, 0.06, 0.028, expression(paste("CWM RTD (g" ," cm" ^{-3},")"))) + 
  ylim(0.025, 0.52) 

RootDiam_mod_plot <- mod_plot_function(root_diam_mod, 0.145, 0.127,0.11,0.1, "CWM root diameter (mm)") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 0.8))

SRL_mod_plot <- mod_plot_function(SLA_mod, 5.55, 4.25, 3.225, 2.51,   expression(paste("CWM SRL (m" ," g" ^{-1},")"))) + 
  scale_y_continuous(trans = "log2", limits = c(2.5, 130)) 

# put abundance gradient figures together 

LES_title <- ggdraw() +
  draw_label("Leaf conservation gradient")

RootCons_title <- ggdraw() +
  draw_label("Root conservation gradient")

Size_title <- ggdraw() +
  draw_label("Size gradients")

Collab_title <- ggdraw() +
  draw_label("Root collaboration gradient")

LES_grid <- plot_grid(LES_title, plot_grid(SLA_mod_plot, LDMC_mod_plot, LeafN_mod_plot, LeafP_mod_plot,
                                           ncol = 2, rel_widths = c(1,1,1,1), labels = c("a", "b", "c", "d"), label_x = 0.1, scale = 0.95),
                      ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCons_grid <- plot_grid(RootCons_title, plot_grid(RootN_mod_plot, RTD_mod_plot, 
                                                     ncol = 1, rel_heights = c(1,1), labels = c("e", "f"), label_x = 0.1, scale = 0.95),
                           ncol = 1, rel_heights = c(0.5,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

Collab_grid <- plot_grid(Collab_title, plot_grid(SRL_mod_plot, RootDiam_mod_plot,
                                                 ncol = 1, rel_heights = c(1,1), labels = c("a", "b"), label_x = 0.1, scale = 0.95),
                         ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

Size_grid <- plot_grid(Size_title, plot_grid(maxheight_mod_plot, maxroot_mod_plot,
                                             ncol = 1, rel_heights = c(1,1), labels = c("c", "d"), label_x = 0.1, scale = 0.95),
                       ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))


# ------------------------------------ co-occurance plots  -------------------------------------------------------


legend <- get_legend(SLA_plot)
legend_height <- get_legend(Height_plot)

SLA_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Specific leaf area`, -1.65, -2) + 
  labs(y = "ln(SLA Introduced) - ln(SLA Native)",
       title = expression(paste("SLA (m" ^{2}," g" ^{-1},")"))) + ylim(-2, 2.5)
LDMC_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Leaf dry matter content`, -0.36, -0.4) + 
  labs(y = "LDMC Introduced -  LDMC Native",
       title = expression(paste("LDMC (g" ," g" ^{-1},")"))) + ylim(-0.4, 0.2)
leafN_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Leaf nitrogen concentration`, -1.05, -1.25) + 
  labs(y = "ln(Leaf N Introduced) -  ln(Leaf N Native)",
       title = expression(paste("Leaf N (mg" ," g" ^{-1},")")))
leafP_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Leaf phosphorus concentration`, -1.725, -2) + 
  labs(y = "ln(Leaf P Introduced) -  ln(Leaf P Native)", 
       title = expression(paste("Leaf P (mg" ," g" ^{-1},")"))) + ylim(-2, 1.75)
rootN_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Root nitrogen concentration`, -1.075, -1.25) + 
  labs(y = "ln(Root N Introduced) - ln(Root N Native)",
       title = expression(paste("Root N (mg" ," g" ^{-1},")"))) + ylim(-1.25, 1.25)
RTD_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Root tissue density`,-0.42, -0.5) + 
  labs(y = "RTD Introduced -  RTD Native", 
       title = expression(paste("RTD (g" ," cm" ^{-3},")"))) + ylim(-0.5, 0.5)
SRL_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Specific root length`, 3.5, 3) + 
  labs(y = "ln(SRL Introduced) - ln(SRL Native)", 
       title = expression(paste("CWM SRL (m" ," g" ^{-1},")"))) + ylim(-3, 4)
RootDiam_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Root diameter`, 1.8, 1.6) + 
  labs(y = "Root diameter (mm)", 
       title = "Root diameter (mm)") + ylim(-1.5, 2)
Height_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Maximum height`, 3.7, 3) + 
  labs(y = "ln(Height Introduced) - ln(Height Native)", 
       title = expression(paste("CWM max height (m"^{},")"))) + ylim(-6, 5)
RootDepth_plot <- diff_mod_function_plot(inv_nat_cwm_er_splits$`Maximum rooting depth`, 3.5, 3) + 
  labs(y = "ln(Depth Introduced) - ln(Depth Native)", 
       title = "Max Rooting Depth (m)") + ylim(-4, 4)

# Grids 
LES_within_comm_grid <- plot_grid(LES_title, plot_grid(SLA_plot, LDMC_plot, leafN_plot, leafP_plot,
                                                       ncol = 2, rel_widths = c(1,1,1,1), labels = c("g", "h", "i", "j"), label_x = 0.1, scale = 0.95),
                                  ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCons_within_comm_grid <- plot_grid(RootCons_title, plot_grid(rootN_plot, RTD_plot, 
                                                                 ncol = 1, rel_heights = c(1,1), labels = c("k", "l"), label_x = 0.1, scale = 0.95),
                                       ncol = 1, rel_heights = c(0.5,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

Collab_within_comm_grid <- plot_grid(Collab_title, plot_grid(SRL_plot, RootDiam_plot,
                                                             ncol = 1, rel_heights = c(1,1), labels = c("e", "f"), label_x = 0.1, scale = 0.95),
                                     ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

Size_within_comm_grid <- plot_grid(Size_title, plot_grid(Height_plot, RootDepth_plot,
                                                         ncol = 1, rel_heights = c(1,1), labels = c("g", "h"), label_x = 0.1, scale = 0.95),
                                   ncol = 1, rel_heights = c(0.5,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

# All figures together 
across_title <- ggdraw() +
  draw_label("TRAIT SHIFTS ALONG GRADIENTS OF INTRODUCED SPECIES ABUNDANCE")

within_title <- ggdraw() +
  draw_label("TRAIT DIFFERENCES BETWEEN CO-OCCURRING NATIVE AND INTRODUCED SPECIES")

# Conservation gradients 
conservation_across <- plot_grid(LES_grid, RootCons_grid, nrow = 1, rel_widths = c(2,1))
conservation_within <- plot_grid(LES_within_comm_grid, RootCons_within_comm_grid, nrow = 1, rel_widths = c(2,1))

# Collab & size 
collab_size_across <- plot_grid(Collab_grid, Size_grid, nrow = 1, rel_widths = c(1,1), align = "v")
collab_size_within <- plot_grid(Collab_within_comm_grid, Size_within_comm_grid, nrow = 1, rel_widths = c(1,1), align = "v")

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/across_within_community_plots_conservationgradients_1.pdf", height = 17, width = 12)
plot_grid(across_title, conservation_across, within_title, conservation_within, legend, ncol = 1, rel_heights = c(0.5,5,0.5,5,1))
dev.off()

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/across_within_community_plots_collabsize_1.pdf", height = 19, width = 9)
plot_grid(across_title, collab_size_across, within_title, collab_size_within, legend_height, ncol = 1, rel_heights = c(0.5,5,0.5,5,1))
dev.off()
