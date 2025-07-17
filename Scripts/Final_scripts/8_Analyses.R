
# project: Trait-Shifts 
# objective: analyses
# author: Magda Garbowski 
# date: November 22, 2023

library(lme4)
library(lmerTest)
library(MuMIn)
library(merTools)

CWM_plot_info_80 <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80.csv")
CWM_plot_info_80_2per <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80_2percent.csv")
inv_nat_cwm_er <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_cooccurring_weighted_means.csv")

# -------------------------------------- trait splits for models ----------------------------------------

CWM_80_splits <- split(CWM_plot_info_80, CWM_plot_info_80$Trait)
CWM_80_splits_2per <- split(CWM_plot_info_80_2per, CWM_plot_info_80_2per$Trait)
inv_nat_cwm_er_splits <- split(inv_nat_cwm_er, inv_nat_cwm_er$Trait)

# --------------------------------------- abundance model function -------------------------------------

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

# ----------------------------------- co-occurrance model function -------------------------------------

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

# ----------------------- abundance gradient output - full dataset  ---------------------------------

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

# --------------------------- abundance gradient output - 2% dataset  -------------------------------

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

