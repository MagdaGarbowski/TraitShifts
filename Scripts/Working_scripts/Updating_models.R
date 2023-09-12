
# updating traitshifts models

library(sjPlot)

dat <- CWM_80_splits$max_975_heightveg_m

dat <- dat[!is.na(dat$EcoRegionLevelI),]
dat$I_relcov_scaled <- scale(dat$I_relcov, center = TRUE, scale = TRUE)
trait <- dat$Trait[1]

mod_2 <- lmer(log(CWM) ~ I_relcov_scaled + (I_relcov_scaled|EcoRegionLevelI) + (1|US_L4CODE), data = dat)
summary(mod_2)
R2_marginal <- round(r.squaredGLMM(mod_2)[1],3)
R2_conditional <- round(r.squaredGLMM(mod_2)[2],3)

# prediction intervals (mod 1)
xmin = min(dat$I_relcov_scaled)
xmax = max(dat$I_relcov_scaled)

new_dat_full <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                            EcoRegionLevelI = c(levels(as.factor(dat$EcoRegionLevelI))))

predict_int_full <- predictInterval(mod_1, new_dat_full, n.sims = 2000,
                                    stat = "median",
                                    type = "linear.prediction",
                                    include.resid.var = TRUE)


# prediction intervals (mod 2)

new_dat_full_2 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                              EcoRegionLevelI = c(levels(as.factor(dat$EcoRegionLevelI))))

predict_int_full_2 <- predictInterval(mod_2, new_dat_full_2, n.sims = 2000,
                                    stat = "median",
                                    type = "linear.prediction",
                                    include.resid.var = TRUE)
coefs <- coef(mod_2)

# find average intercept for ecoregion 4 
avg_ER4 <- mean(coefs$US_L4CODE$`(Intercept)`)
L4_coeffs <- as.data.frame(coefs$US_L4CODE)

target.index <- which(abs(L4_coeffs$`(Intercept)` - avg_ER4) == min(abs(L4_coeffs$`(Intercept)`  - avg_ER4)))
ER4_avg_value <- vector_intercept[target.index]
ER4_avg_L4_code <- rownames(L4_coeffs)[L4_coeffs$`(Intercept)` == ER4_avg_value]

new_dat_full_3 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                              EcoRegionLevelI = c(levels(as.factor(dat$EcoRegionLevelI))),
                              US_L4CODE = "17k")

predict_int_full_3 <- predictInterval(mod_2, new_dat_full_3, n.sims = 2000,
                                      stat = "median",
                                      type = "linear.prediction",
                                      include.resid.var = TRUE)

out_all <- cbind(new_dat_full_3, predict_int_full_3)
out_all_exp <- cbind(out_all[1:3],exp(out_all[4:6]))

# overall trend line 
overall_intercept <- mean(coefs$EcoRegionLevelI$`(Intercept)`)
overall_slope <- mean(coefs$EcoRegionLevelI$I_relcov_scaled)

new_main_dat_3 <- expand.grid(I_relcov_scaled = seq(xmin,xmax,0.005), 
                              EcoRegionLevelI = "overall",
                              US_L4CODE = "17k")

predict_int_main_3 <- predictInterval(mod_2, new_main_dat_3, n.sims = 5000,
                                      stat = "median",
                                      type = "linear.prediction",
                                      include.resid.var = TRUE)

main_all <- cbind(new_main_dat_3, predict_int_main_3)
main_all_exp <- cbind(main_all[1:3],exp(main_all[4:6]))


ggplot() +
  geom_point(data = dat, aes(x = I_relcov_scaled, y = CWM, fill = EcoRegionLevelI, color = EcoRegionLevelI), size = 0.8, alpha = 0.4) +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI == "EASTERN TEMPERATE FORESTS", ], aes(x = I_relcov_scaled, y = fit), color = "#006633") +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI =="GREAT PLAINS",], aes(x = I_relcov_scaled, y = fit), color = "#E69F00", linewidth = 1)  +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI =="MEDITERRANEAN CALIFORNIA",],aes(x = I_relcov_scaled, y = fit), color = "#F0E442", linewidth = 1)  +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI =="NORTH AMERICAN DESERTS",],aes(x = I_relcov_scaled, y = fit), color = "#D55E00", linewidth = 1)  +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI =="NORTHERN FORESTS",], aes(x = I_relcov_scaled, y = fit), color = "#009E73", linewidth = 1)  +
  geom_smooth(data = out_all_exp[out_all_exp$EcoRegionLevelI =="NORTHWESTERN FORESTED MOUNTAINS",], aes(x = I_relcov_scaled, y = fit), color = "#56B4E9", linewidth = 1) + 
  geom_smooth(data = main_all_exp, aes(x = I_relcov_scaled, y = fit), color = "black", linewidth = 1) + 
  geom_ribbon(data = main_all_exp, aes(x = I_relcov_scaled, ymin = lwr, ymax = upr), alpha = 0.25) + 
  labs(title = paste(dat$Trait[1], "\n", "R2_m =", R2_marginal,"; R2_c =", R2_conditional),
       x = "Scaled invasion level", y = paste(dat$Trait[1], "(CWM)")) +
  theme_bw() + 
  scale_x_continuous(limits = c(xmin, xmax), expand = c(0,0))+
  scale_color_manual(values = cb_palette) + 
  scale_fill_manual(values = cb_palette) +
  theme(legend.position = "none", 
        legend.text = element_text(size = 8), 
        legend.spacing.y = unit(0.15, "cm"),
        legend.title = element_blank()) + 
  guides(fill = guide_legend(nrow = 3, byrow = TRUE), 
         color = guide_legend(override.aes = list(alpha = 1))) + 
  scale_y_continuous(trans = "log2", limits = c(0.04, 46))

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
############ Within community analyses ##############################

SLA_dat <- inv_nat_cwm_er_splits$`SLA_mm2/mg`
LDMC_dat <- inv_nat_cwm_er_splits$`LDMC_g/g`
LeafN_dat <- inv_nat_cwm_er_splits$`leafN_mg/g`
LeafP_dat <- inv_nat_cwm_er_splits$`leafP_mg/g`
RootN_dat <- inv_nat_cwm_er_splits$`Root_N_concentration`
RTD_dat <- inv_nat_cwm_er_splits$`Root_tissue_density`

# original model # 
mod_1_function <- function(df){
  TraitNameAbr <- df$TraitNameAbr[1]
  dat_long <- reshape(df, 
                      varying = c("CWM_I", "CWM_N"),
                      v.names = "CWM",
                      idvar = c("Plot"),
                      timevar = "NativeStatus",
                      times = c("I", "N"),
                      direction = "long")
  
  mod_out <- lmer(log(CWM) ~ NativeStatus + (1|EcoRegionLevelI/US_L4NAME) + (1|Plot), data = dat_long)
  sum_mod <- summary(mod_out)
  return(sum_mod)
  
}

mod_1_function(SLA_dat)
mod_1_function(LDMC_dat)
mod_1_function(LeafN_dat)
mod_1_function(LeafP_dat)
mod_1_function(RootN_dat)
mod_1_function(RTD_dat)


# difference model # 

mod_2_function <- function(df){
  TraitNameAbr <- df$TraitNameAbr[1]
  
  if(TraitNameAbr == "SLA_mm2/mg"| TraitNameAbr == "heightveg_m" | TraitNameAbr == "leafP_mg/g" | 
     TraitNameAbr == "leafarea_mm2" | TraitNameAbr == "leafN_mg/g" | TraitNameAbr == "max_rooting_depth_m"| 
     TraitNameAbr == "Mean_Root_diameter" | TraitNameAbr == "Root_N_concentration"| TraitNameAbr =="Specific_root_length"){
    df$difference = log(df$CWM_I) - log(df$CWM_N)
  }
  
  if(TraitNameAbr != "SLA_mm2/mg"& TraitNameAbr != "heightveg_m" & TraitNameAbr != "leafP_mg/g" & 
     TraitNameAbr != "leafarea_mm2" & TraitNameAbr != "leafN_mg/g" & TraitNameAbr != "max_rooting_depth_m"& 
     TraitNameAbr != "Mean_Root_diameter" & TraitNameAbr != "Root_N_concentration" & TraitNameAbr !="Specific_root_length"){
    df$difference = df$CWM_I - df$CWM_N
  }
  
  mod_out <- lmer(difference ~ (1|EcoRegionLevelI/US_L4NAME), data = df)
  sum_mod <- summary(mod_out)
  return(sum_mod)
}

mod_2_function(SLA_dat)
mod_2_function(RootN_dat)

mod_2_function_plot <- function(df, ann_y_1, ann_y_2){
  df = SLA_dat

  TraitNameAbr <- df$TraitNameAbr[1]
  
  if(TraitNameAbr == "SLA_mm2/mg"| TraitNameAbr == "heightveg_m" | TraitNameAbr == "leafP_mg/g" | 
     TraitNameAbr == "leafarea_mm2" | TraitNameAbr == "leafN_mg/g" | TraitNameAbr == "max_rooting_depth_m"| 
     TraitNameAbr == "Mean_Root_diameter" | TraitNameAbr == "Root_N_concentration"| TraitNameAbr =="Specific_root_length"){
    df$difference = log(df$CWM_I) - log(df$CWM_N)
  }
  
  if(TraitNameAbr != "SLA_mm2/mg"& TraitNameAbr != "heightveg_m" & TraitNameAbr != "leafP_mg/g" & 
     TraitNameAbr != "leafarea_mm2" & TraitNameAbr != "leafN_mg/g" & TraitNameAbr != "max_rooting_depth_m"& 
     TraitNameAbr != "Mean_Root_diameter" & TraitNameAbr != "Root_N_concentration"& TraitNameAbr !="Specific_root_length"){
    df$difference = df$CWM_I - df$CWM_N
  }

  mod_out <- lmer(difference ~ (1|EcoRegionLevelI/US_L4NAME), data = df)
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
    geom_point(data = df, aes(x = EcoRegionLevelI, y = difference, color = EcoRegionLevelI), alpha = 0.1, position = position_jitter(width = 0.25)) + 
    geom_linerange(aes(ymin = mean_eff - (1.96 * sd), 
                      ymax = mean_eff + (1.96 * sd)), color = "black", size = 0.6) + 
    geom_point(aes(x = groupID, y = mean_eff, shape = groupID, size = groupID, fill = groupID), color = "black") +
    scale_shape_manual(values = c(22,22,22,22,22,22,22)) +
    scale_color_manual(values = c("#009E73","#CC79A7","#E69F00", "#D55E00", "#56B4E9", "#006633", "black")) +
    scale_fill_manual(values = c("#009E73","#CC79A7","#E69F00", "#D55E00", "#56B4E9", "#006633", "black")) +
    scale_size_manual(values = c(4,4,4,4,4,4,4)) + 
    geom_hline(yintercept = 0, size = 0.6, color = "gray25", linetype = "dashed") +
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
   guides(color = "none") + 
   annotate(geom = "text", x = "NORTHERN FORESTS", y = ann_y_1, label = paste ("n =", n), hjust = 0, size = 3) + 
   annotate(geom = "text", x = "NORTHERN FORESTS", y = ann_y_2, label = paste ("p =", p_val), hjust = 0, size = 3) 

 return(plot_out)
}  

legend <- get_legend(SLA_plot_2)

SLA_plot <- mod_2_function_plot(SLA_dat, 0.85, 0.77) + labs(y = "ln(SLA Introduced) - ln(SLA Native)", 
                                    title = expression(paste("SLA (m" ^{2}," g" ^{-1},")")))
LDMC_plot <- mod_2_function_plot(LDMC_dat, -0.1, -0.14) + labs(y = "LDMC Introduced -  LDMC Native", 
                                                 title = expression(paste("LDMC (g" ," g" ^{-1},")")))
leafN_plot <- mod_2_function_plot(LeafN_dat, 0.39, 0.36) + labs(y = "ln(Leaf N Introduced) -  ln(Leaf N Native)", 
                                                    title = expression(paste("Leaf N (mg" ," g" ^{-1},")")))
leafP_plot <- mod_2_function_plot(LeafP_dat, 0.7, 0.65) + labs(y = "ln(Leaf P Introduced) -  ln(Leaf P Native)", 
                                                    title = expression(paste("Leaf P (mg" ," g" ^{-1},")")))
rootN_plot <- mod_2_function_plot(RootN_dat, 0.185, 0.155) + labs(y = "ln(Root N Introduced) - ln(Root N Native)", 
                                                    title = expression(paste("Root N (mg" ," g" ^{-1},")")))
RTD_plot <- mod_2_function_plot(RTD_dat, 0.325, 0.26) + labs(y = "ln(RTD Introduced) -  ln(RTD Native)", 
                                                    title = expression(paste("RTD (g" ," cm" ^{-3},")")))

# Option 2 
SLA_plot_2 <- mod_2_function_plot(SLA_dat, 2, 1.8) + labs(y = "ln(SLA Introduced) - ln(SLA Native)", 
                                                            title = expression(paste("SLA (m" ^{2}," g" ^{-1},")"))) + ylim(-2, 2.5)
LDMC_plot_2 <- mod_2_function_plot(LDMC_dat, .18, 0.14) + labs(y = "LDMC Introduced -  LDMC Native", 
                                                               title = expression(paste("LDMC (g" ," g" ^{-1},")"))) + ylim(-0.4, 0.2)
leafN_plot_2 <- mod_2_function_plot(LeafN_dat, 1.5, 1.35) + labs(y = "ln(Leaf N Introduced) -  ln(Leaf N Native)", 
                                                                title = expression(paste("Leaf N (mg" ," g" ^{-1},")"))) 
leafP_plot_2 <- mod_2_function_plot(LeafP_dat, 1.5, 1.35) + labs(y = "ln(Leaf P Introduced) -  ln(Leaf P Native)", 
                                                               title = expression(paste("Leaf P (mg" ," g" ^{-1},")"))) + ylim(-2, 1.75)
rootN_plot_2 <- mod_2_function_plot(RootN_dat, 0.85, 0.75) + labs(y = "ln(Root N Introduced) - ln(Root N Native)", 
                                                                  title = expression(paste("Root N (mg" ," g" ^{-1},")"))) + ylim(-1.25, 1.25)
RTD_plot_2 <- mod_2_function_plot(RTD_dat, 0.325, 0.26) + labs(y = "RTD Introduced -  RTD Native", 
                                                             title = expression(paste("RTD (g" ," cm" ^{-3},")"))) + ylim(-0.5, 0.5)


LES_title <- ggdraw() +
  draw_label("Leaf conservation gradient")

RootCons_title <- ggdraw() +
  draw_label("Root conservation gradient")

LES_grid <- plot_grid(LES_title, plot_grid(SLA_plot_2, LDMC_plot_2, leafN_plot_2, leafP_plot_2, ncol = 1, 
                                           rel_widths = c(1,1,1,1), labels = c("a", "b", "c", "d"), scale = 0.95),
                      ncol = 1, rel_heights = c(1,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCons_grid <- plot_grid(RootCons_title, plot_grid(rootN_plot_2, RTD_plot_2, ncol = 1, rel_heights = c(1,1), labels = c("e", "f"), scale = 0.95),
                           ncol = 1, rel_heights = c(1,4))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))


pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/within_community_plots_differences_consgrad_1.pdf", height = 14, width = 5)
plot_grid(plot_grid(LES_grid, RootCons_grid, legend, ncol = 1, rel_heights = c(2,1, 0.25)))
dev.off()

LES_grid_2 <- plot_grid(LES_title, plot_grid(SLA_plot_2, LDMC_plot_2, leafN_plot_2, leafP_plot_2, ncol = 2, rel_widths = c(1,1,1,1), labels = c("a", "b", "c", "d"), scale = 0.95),
                      ncol = 1, rel_heights = c(1,8)) + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

RootCons_grid_2 <- plot_grid(RootCons_title, plot_grid(rootN_plot_2, RTD_plot_2, ncol = 1, rel_heights = c(1,1), labels = c("e", "f"), scale = 0.95),
                           ncol = 1, rel_heights = c(1,8))  + 
  theme(plot.background = element_rect(color = "black", linewidth = 1))

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/within_community_plots_conservationgradients_3.pdf", height = 8, width = 10.5)
plot_grid(plot_grid(LES_grid_2, RootCons_grid_2, nrow = 1, rel_widths = c(2,1)), legend, ncol = 1, rel_heights = c(10,1))
dev.off()
