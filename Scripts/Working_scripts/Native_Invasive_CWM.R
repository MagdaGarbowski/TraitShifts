# Trait-Shifts 
# goal: get CWM for each plot by trait combo 

library(data.table)
library(ggplot2)
library(lme4)
library(emmeans)
library(cowplot)

cb_palette <- c("#006633",  "#E69F00","#D55E00", "#F0E442", "#56B4E9","#009E73",  "#0072B2",  "#CC79A7","#999999" )

# ------------------------------- data ---------------------------------------------

dat <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_10272022_wtraits.csv", 
                           select = c("AcceptedTaxonName", "Plot", "Year", "NativeStatus", "PctCov_100", "rel_100", "TraitNameAbr", "mean")))

coverage <- as.data.frame(fread("/Users/MagdaGarbowski 1/TraitShifts/Generated_Data/SPCIS_TRY_coverage.csv"))

# ------------------------------- prep data  ---------------------------------------------

# select out plots with >80 coverage for traits 
coverage_80 <- coverage[coverage$trait_cov > 80,]

# drop EcoRegions without good coverage 
u_er_plots <- unique(coverage_80[c("EcoRegionLevelI", "Plot")])
table(u_er_plots$EcoRegionLevelI)

coverage_80 <- coverage_80[coverage_80$EcoRegionLevelI %in% c("EASTERN TEMPERATE FORESTS", "GREAT PLAINS", "NORTHERN FORESTS",
                                                              "MEDITERRANEAN CALIFORNIA", "NORTH AMERICAN DESERTS", "NORTHWESTERN FORESTED MOUNTAINS"),]
        
coverage_80_ss <- coverage_80[c("Plot", "Trait", "Year", "trait_cov")]

# merge datasets to select plots from dat that have >80 for specific traits 
dat_80_cov <- merge(dat, coverage_80_ss, by.x = c("Plot", "TraitNameAbr", "Year"), 
                    by.y = c("Plot", "Trait", "Year"), all.y = TRUE)

dat_80_cov_noNI <- dat_80_cov[!dat_80_cov$NativeStatus == "NI",]
dat_80_cov_noNI <- dat_80_cov_noNI[!dat_80_cov_noNI$TraitNameAbr %in% c("Duration", "Growth.Habit", "Mycorrhizal.type"),]






# ------------  recalculate relative cover for natives and invasives alone -------- 

dat_80_splits <- split(dat_80_cov_noNI, list(dat_80_cov_noNI$Plot, dat_80_cov_noNI$Year), drop = TRUE)

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

inv_nat_cwm <- do.call(rbind, lapply( N_I_plots, CWM_function_N_I))

# merge with ecoregion 

eco_r <- unique(coverage[c("Plot", "Year", "EcoRegionLevelI")])

inv_nat_cwm_er <- merge(inv_nat_cwm, eco_r, by = c("Plot", "Year"))

inv_nat_cwm_er_splits <- split(inv_nat_cwm_er, inv_nat_cwm_er$TraitNameAbr)


# need averages for the different ecoregions? 

plot_dat <- function(df){
  dat_long <- reshape(df, 
                      varying = c("CWM_I", "CWM_N"),
                      v.names = "CWM",
                      idvar = c("Plot", "Year"),
                      timevar = "NativeStatus",
                      times = c("I", "N"),
                      direction = "long")
  
  dat_long <- dat_long[!dat_long$CWM == -Inf,]
  
  mod <- lmer(CWM ~ NativeStatus*EcoRegionLevelI + (1|Plot),
              data = dat_long)
  
  est_means_out <- as.data.frame(emmeans(mod, ~EcoRegionLevelI+NativeStatus))
  
  est_means_out_wide <- reshape(est_means_out,  idvar = "EcoRegionLevelI",
                                timevar = "NativeStatus",
                                direction = "wide")
  return(est_means_out_wide)
}

height <- plot_dat(inv_nat_cwm_er_splits$heightveg_m)
SLA <- plot_dat(inv_nat_cwm_er_splits$`SLA_mm2/mg`)
LDMC <- plot_dat(inv_nat_cwm_er_splits$`LDMC_g/g`)
leafN <- plot_dat(inv_nat_cwm_er_splits$`leafN_mg/g`)
root_depth <- plot_dat(inv_nat_cwm_er_splits$max_rooting_depth_m)
RTD <- plot_dat(inv_nat_cwm_er_splits$Root_tissue_density)
SRL <- plot_dat(inv_nat_cwm_er_splits$Specific_root_length)

plot_function <- function(df, means_df, min, max, trait){
  plot <-   ggplot(df, aes (x = CWM_N, y = CWM_I)) + 
    geom_point(aes(color = EcoRegionLevelI), size = 1, alpha = 0.3) + 
    ylim(min, max) +
    xlim(min, max) + 
    geom_point(data = means_df,
               aes(x = emmean.N, y = emmean.I, fill = EcoRegionLevelI),
               shape = 21, size = 4, alpha = 0.8) + 
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    scale_color_manual(values = cb_palette) +
    scale_fill_manual(values = cb_palette) +
    ylab(paste("Introducted CWM", trait)) + 
    xlab(paste("Native CWM", trait)) + 
    theme_bw() + 
    theme(legend.position = "none")
  return(plot)
}

height_plot <- plot_function(inv_nat_cwm_er_splits$heightveg_m, height, 0, 30, "Height")
SLA_plot <- plot_function(inv_nat_cwm_er_splits$`SLA_mm2/mg`, SLA, 0, 65, "SLA_mm2/mg")
LDMC_plot <- plot_function(inv_nat_cwm_er_splits$`LDMC_g/g`, LDMC, 0, 0.65, "LDMC_g/g")
leafN_plot <- plot_function(inv_nat_cwm_er_splits$`leafN_mg/g`, leafN, 0, 50, "leafN_mg/g")
root_depth_plot <- plot_function(inv_nat_cwm_er_splits$max_rooting_depth_m, root_depth, 0, 50, "rooting_depth")
RTD_plot <- plot_function(inv_nat_cwm_er_splits$Root_tissue_density, RTD, 0, 0.9, "RTD")
SRL_plot <- plot_function(inv_nat_cwm_er_splits$Specific_root_length, SRL, 0, 660, "SRL")

legend <- get_legend(height_plot)

pdf(file = "/Users/MagdaGarbowski 1/TraitShifts/Output_Figures/CWM_plots.pdf", height = 12, width = 8)
plot_grid(height_plot, SLA_plot, LDMC_plot, leafN_plot,
          root_depth_plot, RTD_plot, SRL_plot, legend, ncol = 2)
dev.off()

Groot_full <- as.data.frame(fread("/Users/MagdaGarbowski 1/Downloads/GRooTFullVersion.csv", 
                            select = c("referencesAbbreviated", "genus", "species", "traitName","traitValue")))



