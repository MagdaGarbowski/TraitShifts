
# project: Trait-Shifts 
# objective: make maps associated with analyses 
# author: Magda Garbowski 
# date: May 8, 2024

rm(list = ls())

library(ggplot2)
library(data.table)
library(ggplot2)
library(usmap)
library(mapdata)
library(tidyr)
library(sf)

theme_set(theme_void())

CWM_plot_info_80 <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_CWM_80.csv")
inv_nat_cwm_er <- read.csv("/Users/magdagarbowski/TraitShifts/Submitted_datasets/SPCIS_cooccurring_weighted_means.csv")
plot_info <- read.csv("/Users/magdagarbowski/TraitShifts/Data/SPCIS_ecoreg_03182023.csv")

CWM_plot_info_80_wlatlong <- merge(CWM_plot_info_80, plot_info[c("Plot", "Lat", "Long")], by = c("Plot", "Lat", "Long"))

state <- map_data("state")

# create co-occurance column 
CWM_plot_info_80_wlatlong$cooccurance <- ifelse(CWM_plot_info_80_wlatlong$N_relcov > 0 & CWM_plot_info_80_wlatlong$I_relcov > 0, "Yes", "No")

# split by trait for plotting 
splits <- split(CWM_plot_info_80_wlatlong, CWM_plot_info_80_wlatlong$Trait)

map_plot_function <- function(df){
  df$EcoRegionLevelI <- as.factor(as.character(df$EcoRegionLevelI))
  
p_title <- df$Trait[[1]]
n = nrow(df)
n_coocur = nrow(df[df$cooccurance == "Yes",])

plot <- ggplot(data=state, aes(x=long, y=lat, group=group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color="black", fill="white") + 
  ggtitle(p_title)  +
  geom_point(aes(x = Long, y = Lat, group = EcoRegionLevelI, color = EcoRegionLevelI),
             size = 0.5,
             shape = 16,
             alpha = 0.5,
             data = df) +
  scale_color_manual(values = c("#006633",  "#E69F00",  "#F0E442","#D55E00","#009E73",  "#56B4E9")) + 
  labs(x = "Latitude", y = "Longitude") +  
    theme_bw() + 
    theme(legend.position = "none",
          plot.margin = margin(-25,10,-25,10),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          title = element_text(size = 18),
          legend.title = element_blank()) +
  guides(fill = "none", shape = "none", color = guide_legend(override.aes = list(alpha = 1, size = 4), ncol = 3)) +
  annotate(geom = "text", x = -120, y = 28, label = paste ("n =", n), hjust = 0, size = 5) + 
  annotate(geom = "text", x = -120, y = 25.5, label = paste ("n cooccurance =", n_coocur), hjust = 0, size = 5) 
return(plot)
}

SLA <- map_plot_function(splits$`Specific leaf area`)
LDMC <- map_plot_function(splits$`Leaf dry matter content`)
LeafN <- map_plot_function(splits$`Leaf nitrogen concentration`)
LeafP <- map_plot_function(splits$`Leaf phosphorus concentration`)

SRL <- map_plot_function(splits$`Specific root length`)
Diam <- map_plot_function(splits$`Root diameter`)
RootN <- map_plot_function(splits$`Root nitrogen concentration`)
RTD <- map_plot_function(splits$`Root tissue density`)

Height <- map_plot_function(splits$`Maximum height`)
Depth <- map_plot_function(splits$`Maximum rooting depth`)

legend <- cowplot::get_plot_component(SLA, "guide-box-bottom", return_all = TRUE)

pdf(file = "/Users/magdagarbowski/Desktop/maps_les_res.pdf", height = 17, width = 14)
plot_grid(plot_grid(SLA, LDMC, LeafN, LeafP, RootN, RTD, ncol = 2,
                    rel_heights = c(1,1,1,1,1,1), labels = c("a", "b", "c", "d", "e", "f")),
          legend, rel_heights = c(8,1), ncol = 1)
dev.off()

pdf(file = "/Users/magdagarbowski/Desktop/maps_roots.pdf", height = 10, width = 14)
plot_grid(plot_grid(SRL, Diam, RootN, RTD, ncol = 2, rel_heights = c(1,1,1,1), labels = c("a", "b", "c", "d")),
          legend, rel_heights = c(6,1), ncol = 1)
dev.off()

pdf(file = "/Users/magdagarbowski/Desktop/maps_collab_size.pdf", height = 10, width = 14)
plot_grid(plot_grid(SRL, Diam, Height, Depth, ncol = 2, labels = c("a", "b", "c", "d")), 
          legend, rel_heights = c(6,1), ncol = 1)
dev.off()







