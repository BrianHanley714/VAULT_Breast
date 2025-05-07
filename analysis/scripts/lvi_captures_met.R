# Create a plot looking at the proportion of subclonal variants per driver type

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(epitools)


# ORGANISE PATHS ----------------------------------------------------------
BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast/"
OUT_DIR = file.path(BASE, "analysis", "figures")


# BACKGROUND --------------------------------------------------------------
# 6/8 cases which captured the metastasising subclone contained LVI, while this was 1/8 for cases which did not capture it.
# hence, the odds ratio was calculate as follows:
oddsratio(c(6,2,1,7))


# GRAPHICAL PARAMETERS ----------------------------------------------------
rs_col = "#54278F"
sr_col = "#33A02Cb7"

# CREATE DATAFRAME --------------------------------------------------------
df = as.data.frame(cbind(c(6,2,1,7), 
      c("Present", "Absent", "Present", "Absent"),
      c("Met Captured", "Met Captured", "Met Not Captured", "Met Not Captured")))

names(df) = c("n", "LVI", "met")
df$n = as.numeric(df$n)
df$LVI = factor(df$LVI, levels = c("Absent", "Present"))



# DRAW PLOT ---------------------------------------------------------------
df%>%
  ggplot(aes(met, n, fill = LVI))+
  geom_col()+
  theme_classic(base_size = 30)+
  theme(legend.position = "top",
        axis.title.x = element_blank())+
  scale_fill_manual(values = c(sr_col, rs_col))+
  ylab("SingReg Count")


# WRITE PLOT --------------------------------------------------------------
ggsave(file.path(OUT_DIR, "Extended_Data_Figure7_LVIregions_met.png"))


