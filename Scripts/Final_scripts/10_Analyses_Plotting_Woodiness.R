# project: Trait-Shifts 
# objective: woody vs herbacious analyses and figures (for revisions)
# author: Magda Garbowski 
# date: May 3, 2024

library(data.table)
library(lme4)
library(lmerTest)
library(MuMIn)
library(merTools)
library(cowplot)
library(gridGraphics)
library(ggplot2)

# --------------------------------- data ---------------------------------------------

woodiness_dat <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/Woodiness_by_plot.csv")
woodiness_CWM_dat <-read.csv( "/Users/magdagarbowski/TraitShifts/Submitted_datasets/Woodiness_CWM.csv")

cb_palette <- c("#009E73","#CC79A7","#E69F00", "#D55E00", "#56B4E9", "#006633", "black")

# --------------------------------- models ---------------------------------------------

# abundance gradients 
mod_function <-function(dat){
  dat <- dat[!is.na(dat$EcoRegionLevelI),]
  dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
  dat$unscaled <- dat$I_relcov_scaled * attr(dat$I_relcov_scaled, "scaled:scale") + attr(dat$I_relcov_scaled, "scaled:center")

  mod <- lmer(rel_woody ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
  
  summary_mod <- summary(mod)
  coeffs_out <- coef(mod)
  pvalue <- summary_mod$coefficients[[10]]
  R2_marginal <- as.numeric(round(r.squaredGLMM(mod)[1],3))
  R2_conditional <- as.numeric(round(r.squaredGLMM(mod)[2],3))
  
  ls <- list(dat, R2_marginal, R2_conditional, mod, pvalue, summary_mod, coeffs_out)
  return(ls)
}

# --------------------------------- plotting ---------------------------------------------

mod_plot_function <- function(ls, y1, y2, y3, y4, title_lab){
  # data for plot 
  dat = ls[[1]]
  R2_marginal <- ls[[2]]
  R2_conditional <- ls[[3]]
  mod = ls[[4]]
  p = ls[[5]]
  p_val = paste(round(ls[[5]],4))
  trait = "Woodiness"
  
  scale_att_1 <- attr(dat$I_relcov_scaled, "scaled:scale")
  scale_center <- attr(dat$I_relcov_scaled, "scaled:center")
  dat$unscaled <- dat$I_relcov_scaled * attr(dat$I_relcov_scaled, "scaled:scale") + attr(dat$I_relcov_scaled, "scaled:center")
  xmin <- min(dat$I_relcov_scaled)
  xmax <- max(dat$I_relcov_scaled)
  
  R2_M_label = paste("R [M]^2 ==",R2_marginal)
  R2_C_label = paste("R [C]^2 ==",R2_conditional)
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
  
  out_all <- cbind(out_all[,c(1,2,3,7)],out_all[4:6])

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
  
  main_all <- cbind(main_all[,c(1,2,3,7)], main_all[4:6])

  # get min and max for x-axis 
  xmin_unscaled = min(dat$unscaled)
  xmax_unscaled = max(dat$unscaled)
  
  # add "overall" category to dat
  overall_dat <- dat[c(1),]
  overall_dat[c(1),] <- NA
  overall_dat$EcoRegionLevelI <- "OVERALL"
  
  dat <- rbind(dat, overall_dat)
  
  plot_2 <- ggplot() + 
    geom_point(data = dat, aes(x = unscaled, y = rel_woody, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.8, alpha = 0.2) + 
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
    scale_y_continuous( expand = c(0,0))+
    scale_color_manual(values = cb_palette) + 
    scale_fill_manual(values = c(cb_palette, "black")) +
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
          plot.margin = margin(10,10,10,10),
          plot.title = element_text(size = 12, hjust = 0.5)) + 
    guides(fill = "none", shape = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4), ncol = 2)) +
    annotate(geom = "text", x = 75, y = y1, label = R2_M_label, parse = TRUE, hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = 75, y = y2, label = R2_C_label, parse = TRUE, hjust = 0, size = 3.5) +
    annotate(geom = "text", x = 75, y = y3, label = p_label, parse = TRUE, hjust = 0, size = 3.5) + 
    annotate(geom = "text", x = 75, y = y4, label = n_label, parse = TRUE, hjust = 0, size = 3.5) 
  return(plot_2)
}

woody_props <- mod_plot_function(GH_pros_mod, 0.90, 0.85, 0.80, 0.75, expression(paste("PROPORTION WOODY"))) + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0,100))

props_legend <- cowplot::get_plot_component(woody_props, "guide-box-bottom", return_all = TRUE)


# ----------------------------- community weighted means  --------------------------------------

mod_function_2 <-function(dat, x_1, x_2){
  dat_ss <- dat[!is.na(dat$EcoRegionLevelI),]
  dat_ss$I_relcov_scaled <- scale(dat_ss$I_relcov, center = TRUE, scale = TRUE)
  dat_ss$unscaled <- dat_ss$I_relcov_scaled * attr(dat_ss$I_relcov_scaled, "scaled:scale") + attr(dat_ss$I_relcov_scaled, "scaled:center")
  trait <- dat_ss$Trait[1]
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    mod_x1 <- lmer(log(x_1) ~ I_relcov_scaled + (1|EcoRegionLevelI) + (1|US_L4CODE), data = dat_ss)
    mod_x2 <- lmer(log(x_2) ~ I_relcov_scaled + (1|EcoRegionLevelI) + (1|US_L4CODE), data = dat_ss)

  }
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    mod_x1 <- lmer(x_1 ~ I_relcov_scaled + (1|EcoRegionLevelI) + (1|US_L4CODE), data = dat_ss)
    mod_x2 <- lmer(x_2 ~ I_relcov_scaled + (1|EcoRegionLevelI) + (1|US_L4CODE), data = dat_ss)
    
  }
  summary_mod_x1 <- summary(mod_x1)
  coeffs_out_x1 <- coef(mod_x1)
  pvalue_x1 <- summary_mod_x1$coefficients[[10]]
  R2_marginal_x1 <- as.numeric(round(r.squaredGLMM(mod_x1)[1],3))
  R2_conditional_x1 <- as.numeric(round(r.squaredGLMM(mod_x1)[2],3))
  
  summary_mod_x2 <- summary(mod_x2)
  coeffs_out_x2 <- coef(mod_x2)
  pvalue_x2 <- summary_mod_x2$coefficients[[10]]
  R2_marginal_x2 <- as.numeric(round(r.squaredGLMM(mod_x2)[1],3))
  R2_conditional_x2 <- as.numeric(round(r.squaredGLMM(mod_x2)[2],3))
  
  ls <- list(dat_ss, 
             R2_marginal_x1, R2_conditional_x1, mod_x1, pvalue_x1, summary_mod_x1, coeffs_out_x1,
             R2_marginal_x2, R2_conditional_x2, mod_x2, pvalue_x2, summary_mod_x2, coeffs_out_x2)
  return(ls)
}

# plotting function - abundance gradients - woody and herbacious ~ indroduced sps abundance 
mod_plot_function_2 <- function(ls, y1, y2, y3, y4, title_lab, color_1, color_2){
  # data for plot 
  dat = ls[[1]]
  x1 = dat[[5]]
  x2 = dat[[6]]
  trait = dat$Trait[[1]]
  
  # Herbacious
  R2_marginal_x1 <- ls[[2]]
  R2_conditional_x1 <- ls[[3]]
  mod_x1 = ls[[4]]
  p_x1 = ls[[5]]
  p_val_x1 = paste(round(ls[[5]],4))
  
  # Woody
  R2_marginal_x2 <- ls[[8]]
  R2_conditional_x2 <- ls[[9]]
  mod_x2 = ls[[10]]
  p_x2 = ls[[11]]
  p_val_x2 = paste(round(ls[[11]],4))
  
  scale_att_1 <- attr(dat$I_relcov_scaled, "scaled:scale")
  scale_center <- attr(dat$I_relcov_scaled, "scaled:center")
  dat$unscaled <- dat$I_relcov_scaled * attr(dat$I_relcov_scaled, "scaled:scale") + attr(dat$I_relcov_scaled, "scaled:center")
  xmin <- min(dat$I_relcov_scaled)
  xmax <- max(dat$I_relcov_scaled)
  
  R2_M_label_x1 = paste("R [M]^2 ==",R2_marginal_x1)
  R2_C_label_x1 = paste("R [C]^2 ==",R2_conditional_x1)
  R2_M_label_x2 = paste("R [M]^2 ==",R2_marginal_x2)
  R2_C_label_x2 = paste("R [C]^2 ==",R2_conditional_x2)  
  
  n_plot = nrow(dat)
  n_label = paste("n ==", n_plot)
  
  p_label_x1 = ifelse(p_x1 < 0.0001, paste("p < 0.0001"), paste("p ==", p_val_x1))
  p_label_x2 = ifelse(p_x2 < 0.0001, paste("p < 0.0001"), paste("p ==", p_val_x2))
  
  # simulate data, set random-effect of ER4 to mean, get predictions 
  # introduced species 
  coefs_mod_x1 <- coef(mod_x1)
  avg_ER4_mod_x1 <- mean(coefs_mod_x1$US_L4CODE$`(Intercept)`)
  L4_coeffs_mod_x1 <- as.data.frame(coefs_mod_x1$US_L4CODE)
  
  target.index_mod_x1 <- which(abs(L4_coeffs_mod_x1$`(Intercept)` - avg_ER4_mod_x1) == min(abs(L4_coeffs_mod_x1$`(Intercept)`  - avg_ER4_mod_x1)))
  ER4_avg_value_mod_x1 <- L4_coeffs_mod_x1$`(Intercept)`[target.index_mod_x1]
  ER4_avg_L4_code_mod_x1 <- rownames(L4_coeffs_mod_x1)[L4_coeffs_mod_x1$`(Intercept)` == ER4_avg_value_mod_x1]
  
  # simulate data, set random-effect of ER4 to mean, get predictions 
  # native species 
  coefs_mod_x2 <- coef(mod_x2)
  avg_ER4_mod_x2 <- mean(coefs_mod_x2$US_L4CODE$`(Intercept)`)
  L4_coeffs_mod_x2 <- as.data.frame(coefs_mod_x2$US_L4CODE)
  
  target.index_mod_x2 <- which(abs(L4_coeffs_mod_x2$`(Intercept)` - avg_ER4_mod_x2) == min(abs(L4_coeffs_mod_x2$`(Intercept)`  - avg_ER4_mod_x2)))
  ER4_avg_value_mod_x2 <- L4_coeffs_mod_x2$`(Intercept)`[target.index_mod_x2]
  ER4_avg_L4_code_mod_x2 <- rownames(L4_coeffs_mod_x2)[L4_coeffs_mod_x2$`(Intercept)` == ER4_avg_value_mod_x2]
  
  # get overall trend line  
  # introduced species 
  overall_mod_x1ntercept_mod_x1 <- mean(coefs_mod_x1$EcoRegionLevelI$`(Intercept)`)
  overall_slope_mod_x1 <- mean(coefs_mod_x1$EcoRegionLevelI$I_relcov_scaled)
  
  new_main_dat_mod_x1 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                                    EcoRegionLevelI = "overall",
                                    US_L4CODE = ER4_avg_L4_code_mod_x1)
  
  predict_mod_x1nt_main_mod_x1 <- predictInterval(mod_x1, new_main_dat_mod_x1, n.sims = 2000,
                                                stat = "mean",
                                                type = "linear.prediction",
                                                include.resid.var = FALSE)
  
  main_all_mod_x1 <- cbind(new_main_dat_mod_x1, predict_mod_x1nt_main_mod_x1)
  main_all_mod_x1$unscaled <- main_all_mod_x1$I_relcov_scaled * scale_att_1 + scale_center
  main_all_mod_x1$unscaled <- round(main_all_mod_x1$unscaled, 6)
  
  # get overall trend line  
  overall_mod_x2ntercept_mod_x2 <- mean(coefs_mod_x2$EcoRegionLevelI$`(Intercept)`)
  overall_slope_mod_x2 <- mean(coefs_mod_x2$EcoRegionLevelI$I_relcov_scaled)
  
  new_main_dat_mod_x2 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                                    EcoRegionLevelI = "overall",
                                    US_L4CODE = ER4_avg_L4_code_mod_x2)
  
  predict_mod_x2nt_main_mod_x2 <- predictInterval(mod_x2, new_main_dat_mod_x2, n.sims = 2000,
                                                stat = "mean",
                                                type = "linear.prediction",
                                                include.resid.var = FALSE)
  
  main_all_mod_x2 <- cbind(new_main_dat_mod_x2, predict_mod_x2nt_main_mod_x2)
  main_all_mod_x2$unscaled <- main_all_mod_x2$I_relcov_scaled * scale_att_1 + scale_center
  main_all_mod_x2$unscaled <- round(main_all_mod_x2$unscaled, 6)
  
  # exponentiate overall trend line for correct traits 
  if(trait == "Specific leaf area"| trait == "Maximum height" | trait == "Leaf phosphorus concentration" | 
     trait == "Leaf nitrogen concentration" | trait == "Maximum rooting depth"| trait == "Root diameter" | 
     trait == "Root nitrogen concentration"| trait =="Specific root length"){
    main_all_mod_x1 <- cbind(main_all_mod_x1[,c(1,2,3,7)], exp(main_all_mod_x1[4:6]))
    main_all_mod_x2 <- cbind(main_all_mod_x2[,c(1,2,3,7)], exp(main_all_mod_x2[4:6]))
  }
  
  if(trait != "Specific leaf area" & trait != "Maximum height" & trait != "Leaf phosphorus concentration" & 
     trait != "Leaf nitrogen concentration" & trait != "Maximum rooting depth"& 
     trait != "Root diameter" & trait != "Root nitrogen concentration"& trait !="Specific root length"){
    main_all_mod_x1 <- cbind(main_all_mod_x1[,c(1,2,3,7)], main_all_mod_x1[4:6])
    main_allwoodiness_mod_x2 <- cbind(main_all_mod_x2[,c(1,2,3,7)], main_all_mod_x2[4:6])
    
  }
  
  # get min and max for x-axis 
  xmin_unscaled = min(dat$unscaled)
  xmax_unscaled = max(dat$unscaled)
  
  plot_2 <- ggplot() + 
    geom_point(data = dat, aes(x = unscaled, y = x1, fill = EcoRegionLevelI, color = color_1), size = 0.8, alpha = 0.1) +
    geom_point(data = dat, aes(x = unscaled, y = x2, fill = EcoRegionLevelI, color = color_2), size = 0.8, alpha = 0.1) +
    geom_point(data = main_all_mod_x1, aes(x = unscaled, y = fit), color = color_1, size = 1) + 
    geom_ribbon(data = main_all_mod_x1, aes(x = unscaled, ymin = lwr, ymax = upr), alpha = 0.3) + 
    geom_point(data = main_all_mod_x2, aes(x = unscaled, y = fit), color = color_2, size = 1) + 
    geom_ribbon(data = main_all_mod_x2, aes(x = unscaled, ymin = lwr, ymax = upr), alpha = 0.3) + 
    theme_bw() + 
    scale_x_continuous(limits = c(0, 100), expand = c(0,0))+
    scale_color_manual(labels = c("WOODY", "HERBACIOUS"),
                       values = c(color_2, color_1)) + 
    labs(title = title_lab,
         x = "Relative abundance of introduced species",
         y = title_lab) +
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.text = element_text(size = 10), 
          legend.spacing.y = unit(0.15, "cm"),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(12,12,12,12),
          plot.title = element_text(size = 12, hjust = 0.5)) + 
    guides(fill = "none", shape = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4), ncol = 2)) +
    annotate(geom = "text", x = 35, y = y1, label = R2_M_label_x1, parse = TRUE, hjust = 0, size = 3, color = color_1) + 
    annotate(geom = "text", x = 35, y = y2, label = R2_C_label_x1, parse = TRUE, hjust = 0, size = 3, color = color_1) +
    annotate(geom = "text", x = 35, y = y3, label = p_label_x1, parse = TRUE, hjust = 0, size = 3, color = color_1) + 
    annotate(geom = "text", x = 35, y = y4, label = n_label, parse = TRUE, hjust = 0, size = 3, color = color_1) + 
    annotate(geom = "text", x = 75, y = y1, label = R2_M_label_x2, parse = TRUE, hjust = 0, size = 3, color = color_2) + 
    annotate(geom = "text", x = 75, y = y2, label = R2_C_label_x2, parse = TRUE, hjust = 0, size = 3, color = color_2) +
    annotate(geom = "text", x = 75, y = y3, label = p_label_x2, parse = TRUE, hjust = 0, size = 3, color = color_2) + 
    annotate(geom = "text", x = 75, y = y4, label = n_label, parse = TRUE, hjust = 0, size = 3, color = color_2) 
  return(plot_2)
}


# ----------------------------------------------- ANALYSES -------------------------------------------------------

# CWM of woody vs. herbacious 
CWM_H_W_splits <- split(woodiness_CWM_dat, woodiness_CWM_dat$Trait)

SLA_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Specific leaf area`, 
                                    CWM_H_W_splits$`Specific leaf area`$CWM_H, CWM_H_W_splits$`Specific leaf area`$CWM_W)

LDMC_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Leaf dry matter content`,
                                     CWM_H_W_splits$`Leaf dry matter content`$CWM_H, CWM_H_W_splits$`Leaf dry matter content`$CWM_W)

leaf_N_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Leaf nitrogen concentration`, 
                                       CWM_H_W_splits$`Leaf nitrogen concentration`$CWM_H, CWM_H_W_splits$`Leaf nitrogen concentration`$CWM_W)

leaf_P_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Leaf phosphorus concentration`, 
                                       CWM_H_W_splits$`Leaf phosphorus concentration`$CWM_H, CWM_H_W_splits$`Leaf phosphorus concentration`$CWM_W)

root_diam_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Root diameter`, 
                                          CWM_H_W_splits$`Root diameter`$CWM_H, CWM_H_W_splits$`Root diameter`$CWM_W)

RTD_woodiness_mod<- mod_function_2(CWM_H_W_splits$`Root tissue density`, 
                                   CWM_H_W_splits$`Root tissue density`$CWM_H, CWM_H_W_splits$`Root tissue density`$CWM_W) 

root_N_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Root nitrogen concentration`,
                                       CWM_H_W_splits$`Root nitrogen concentration`$CWM_H, CWM_H_W_splits$`Root nitrogen concentration`$CWM_W) 

SRL_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Specific root length`, 
                                    CWM_H_W_splits$`Specific root length`$CWM_H, CWM_H_W_splits$`Specific root length`$CWM_W) 

max_height_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Maximum height`,
                                           CWM_H_W_splits$`Maximum height`$CWM_H, CWM_H_W_splits$`Maximum height`$CWM_W)

root_maxdepth_woodiness_mod <- mod_function_2(CWM_H_W_splits$`Maximum rooting depth`, 
                                              CWM_H_W_splits$`Maximum rooting depth`$CWM_H, CWM_H_W_splits$`Maximum rooting depth`$CWM_W)

# ------------------------------------------------ FIGURES ----------------------------------------------

maxheight_woodiness_mod_plot <- mod_plot_function_2(max_height_woodiness_mod, 0.26, 0.17, 0.115, 0.085, expression(paste("CWM max height (m"^{},")")), "palegreen4", "burlywood4") +
  scale_y_continuous(trans = "log2", limits = c(0.07, 60),  expand = c(0,0))

maxroot_woodiness_mod_plot <- mod_plot_function_2(root_maxdepth_woodiness_mod, 0.26, 0.17, 0.115, 0.085,"CWM max rooting depth (m)", "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(0.07, 60),  expand = c(0,0)) 

SLA_woodiness_mod_plot <- mod_plot_function_2(SLA_woodiness_mod, 6, 4.7, 3.9, 3.3,  expression(paste("CWM SLA (m" ^{2}," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(3, 60),  expand = c(0,0)) 

LDMC_woodiness_mod_plot <- mod_plot_function_2(LDMC_woodiness_mod, 0.175, 0.148, 0.128, 0.11, expression(paste("CWM LDMC (g" ," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(limits = c(0.1, 0.52),  expand = c(0,0))

LeafN_woodiness_mod_plot <- mod_plot_function_2(leaf_N_woodiness_mod, 12,10.6,9.4,8.5, expression(paste("CWM Leaf N (mg" ," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(8, 50),  expand = c(0,0)) 

LeafP_woodiness_mod_plot <- mod_plot_function_2(leaf_P_woodiness_mod, 0.78, 0.68,0.6,0.53, expression(paste("CWM Leaf P (mg" ," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(0.5, 5),  expand = c(0,0)) 

RootN_woodiness_mod_plot <- mod_plot_function_2(root_N_woodiness_mod, 7.1,6.1,5.4,4.8, expression(paste("CWM Root N (mg" ," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(4.5, 32),  expand = c(0,0))

RTD_woodiness_mod_plot <- mod_plot_function_2(RTD_woodiness_mod, 0.11, 0.08, 0.055, 0.035, expression(paste("CWM RTD (g" ," cm" ^{-3},")")), "palegreen4", "burlywood4") + 
   scale_y_continuous(limits = c(0.025, 0.52),  expand = c(0,0))

RootDiam_woodiness_mod_plot <- mod_plot_function_2(root_diam_woodiness_mod, 0.14, 0.125, 0.115,0.105, "CWM root diameter (mm)", "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(0.1, 0.8),  expand = c(0,0))

SRL_woodiness_mod_plot <- mod_plot_function_2(SLA_woodiness_mod, 5.85, 4.4, 3.5, 2.81,   expression(paste("CWM SRL (m" ," g" ^{-1},")")), "palegreen4", "burlywood4") + 
  scale_y_continuous(trans = "log2", limits = c(2.5, 130),  expand = c(0,0)) 

woodiness_legend <- cowplot::get_plot_component(maxheight_woodiness_mod_plot, "guide-box-bottom", return_all = TRUE)


# ----------------------------------------- FIGURE GRIDS  -------------------------------------------------------
props_title <- ggdraw() + draw_label(" ")

woodiness_cwm_title <- ggdraw() + draw_label("CWM TRAITS OF WOODY AND HERBACIOUS SPECIES")

woodiness_grid <- plot_grid(plot_grid(SLA_woodiness_mod_plot,  LeafN_woodiness_mod_plot, 
                             RootN_woodiness_mod_plot, SRL_woodiness_mod_plot, 
                             maxheight_woodiness_mod_plot, maxroot_woodiness_mod_plot, 
                             ncol = 3, labels = c("b", "c", "d", "e", "f", "g")), 
                            woodiness_legend, ncol = 1, rel_heights = c(5, 0.25))

woodiness_grid_2 <- plot_grid(woodiness_cwm_title, woodiness_grid, ncol = 1, rel_heights = c(0.5,6))

props_grid <- plot_grid(woody_props, props_legend, ncol = 1, rel_heights = c(3.5, 1.5), labels = c("a", ""))

props_grid_2 <- plot_grid(props_title, props_grid, ncol = 1, rel_heights = c(0.5, 5))

woodiness_SI <- plot_grid(plot_grid(LDMC_woodiness_mod_plot,  LeafP_woodiness_mod_plot, 
                                    RTD_woodiness_mod_plot, RootDiam_woodiness_mod_plot,
                                    ncol = 2, labels = c("a", "b", "c", "d")), 
                          woodiness_legend, ncol = 1, rel_heights = c(5, 0.25))

# ----------------------------------------- OUTPUT PLOTS  -------------------------------------------------------

pdf(file = "/Users/magdagarbowski/Desktop/combined_woodiness.pdf", height = 8.2, width = 18)
plot_grid(props_grid_2, woodiness_grid_2, ncol = 2, rel_widths = c(1.5, 3))
dev.off()

pdf(file = "/Users/magdagarbowski/Desktop/woodiness_SI.pdf", height = 8.5, width = 9)
woodiness_SI
dev.off()
