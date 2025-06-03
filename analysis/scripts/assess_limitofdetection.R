# Density histogram of mutational cancer cell fraction for simulated data

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(scales)

# PATHS -------------------------------------------------------------------
run = "run001"
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
input_dir = file.path(BASE, "data/simulations", run, "out/")


# LOAD DATA ---------------------------------------------------------------
if(file.exists(input_dir)){
  ccf_comparisons_rep = read.delim(paste0(input_dir, "ccf_comparisons_data_repseq.txt"))
  ccf_comparisons_sr = read.delim(paste0(input_dir, "ccf_comparisons_data_sr.txt"))
}


# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"



# GET PROPORTION AND SUM VARIANTS CAPTURED --------------------------------
thresholds = 10^seq(-6, 0, by = .01) # thresholds between 1 and 10^-6

# proportion and sum of variants per single region after 1000 runs
proportions_sr  = sapply(thresholds, function(threshold) mean(ccf_comparisons_sr$ccfs > threshold))
sum_sr  = sapply(thresholds, function(threshold) sum(ccf_comparisons_sr$ccfs > threshold))

# proportion and sum of variants per RepSeq sample after 1000 runs
proportions  = sapply(thresholds, function(threshold) mean(ccf_comparisons_rep$ccfs > threshold))
sum  = sapply(thresholds, function(threshold) sum(ccf_comparisons_rep$ccfs > threshold))

# DRAW PLOT ---------------------------------------------------------------
plot_grid(
  as_tibble(rbind(cbind(proportions_sr, thresholds, "SingReg"), cbind(proportions, thresholds, "RepSamp")))%>%
    mutate(Sample = V3)%>%
    ggplot(aes(as.numeric(thresholds), as.numeric(proportions_sr), col = Sample))+
    geom_smooth()+
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)),
      
    )+
    ylim(0,1)+
    theme_classic(base_size = 30)+
    guides(col = guide_legend(override.aes = list(shape = 16))) +
    scale_color_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
    theme(axis.title.x = element_blank(),
          legend.position = "none")+
    ylab("Subclones/Sample"),
  
  as_tibble(rbind(cbind(sum_sr, thresholds, "SingReg"), cbind(sum, thresholds, "RepSamp")))%>%
    mutate(Sample = V3)%>%
    ggplot(aes(as.numeric(thresholds), as.numeric(sum_sr), fill = Sample))+
    scale_x_log10(    breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))+
    geom_col()+
    theme_classic(base_size = 30)+
    scale_fill_manual(values = c("RepSamp" = rs_col, "SingReg" = sr_col))+
    xlab("Limit of Detection")+
    ylab("Subclones (total)")+
    theme(legend.position = "bottom"),
  ncol = 1,rel_heights = c(1,1.3)
)


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Figure2B_LODplots.png"))     
