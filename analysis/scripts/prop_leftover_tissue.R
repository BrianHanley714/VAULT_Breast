# Compare leftover tissue to other sample types

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(highcharter)
library(ggpubr)
library(htmlwidgets)

# PATHS -------------------------------------------------------------------
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
TUM_COUNTS = file.path(BASE, "data", "image_analysis", "tumour_cell_counts.tsv")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))


# LOAD DATA ---------------------------------------------------------------
annotations_tumour = read.delim(TUM_COUNTS)



# Graphical Parameters ----------------------------------------------------
thickness_ffpe = 4000 # in microns
col_grams = "#54278F"
rs_col = "#54278F"
col_no_tumourcells = "#33A02Cb7"
sr_col = "#33A02Cb7"



# Organise data -----------------------------------------------------------
annotations_tumour = annotations_tumour%>%
  mutate(numtumourcells_ERWSI = num_tumourcells,
         numtumourcells_50um = tumcell_perum3 * Area.µm.2*50,
         numtumourcells_FFPEmol = tumcell_perum3*vol_ERblock,
         numtumourcells_FFPEother = tumcell_perum3*blocks_vol,
         leftover_tumour_vol = Tumour_Weight*(10^12),
         numtumourcells_leftover = totalcell_perum3* leftover_tumour_vol*purity_homogenate,
         mass_ERWSI = Area.µm.2*3/(10^12),
         mass_50um = Area.µm.2*50/(10^12),
         mass_FFPEmol = vol_ERblock/(10^12),
         mass_FFPEother = blocks_vol/(10^12),
         mass_leftover = Tumour_Weight
  )


plotting_dataframe = (annotations_tumour%>%
                        dplyr::select(starts_with(c("Trial_id", "numtumour", "mass")))%>%
                        pivot_longer(cols = starts_with(c("numtumour", "mass")), names_to ="sample_type", values_to =  "number_tumourcells" )%>%
                        mutate(sample = sub(".*_", "", sample_type),
                               col = factor(sub("_.*", "", sample_type), levels = c("numtumourcells", "mass")),
                               sample = case_when(sample == "50um" ~ "Paraffin curls",
                                                  sample == "ERWSI" ~ "IHC slide",
                                                  sample == "FFPEmol" ~ "FFPE for molecular",
                                                  sample == "FFPEother" ~ "All other FFPEs",
                                                  sample == "leftover" ~ "Leftover tissue",
                               ),
                               sample = factor(sample, levels = c("IHC slide", 
                                                                  "Paraffin curls", 
                                                                  "FFPE for molecular",
                                                                  "All other FFPEs",
                                                                  "Leftover tissue"
                               )))%>%
                        mutate(number_tumourcells = if_else(col == "mass", number_tumourcells*10^8, number_tumourcells)))


# Output plots ------------------------------------------------------------
p = plotting_dataframe%>%
  ggplot(aes(sample, number_tumourcells, fill = col))+
  geom_boxplot()+
  scale_y_log10(
    labels = scales::label_math(expr = 10^.x, format = log10),
    sec.axis = sec_axis(~ . / 10^8, name = "Tissue Mass (g)", 
                        labels = scales::label_math(expr = 10^.x, format = log10)))+
  scale_fill_manual(values = c("mass" = col_grams, "numtumourcells" = col_no_tumourcells))+
  theme_classic(base_size = 40)+
  theme(axis.title.x = element_blank())+
  rot.lab()+
  theme(
    axis.title.y = element_text(color = col_no_tumourcells, size = 30, face = "bold"),
    axis.title.y.right = element_text(color = col_grams, size = 30, face = "bold"),
    axis.text.x = element_text(color = "black", size = 20, face = "bold"),
    legend.position = "none"
  )+
  ylab("No. Tumour Cells")



p+
  stat_compare_means(data = subset(plotting_dataframe, col == "mass"),
                     aes(group = sample),
                     #method = "t.test",  # Can also use "wilcox.test"
                     label = "p.signif",
                     comparisons = list(c("IHC slide", "Leftover tissue"),
                                        c("Paraffin curls", "Leftover tissue"),
                                        c("FFPE for molecular", "Leftover tissue"),
                                        c("All other FFPEs", "Leftover tissue")), 
                     size = 6,  # Increase text size
                     bracket.size = 1, 
                     label.y.npc = 0.95 )

ggsave(file.path(OUT_DIR, "Figure1E_boxplot.png"))

pie_df = plotting_dataframe%>%
  mutate(sample = sample_type,
         sample = case_when(sample == "50um" ~ "Paraffin curls",
                            sample == "mass_ERWSI" ~ "IHC slide",
                            sample == "mass_FFPEmol" ~ "FFPE for molecular",
                            sample == "mass_FFPEother" ~ "All other FFPEs",
                            sample == "mass_leftover" ~ "Leftover tissue"))%>%
  filter(col == "mass")%>%
  dplyr::select(Trial_id, sample_type, sample, number_tumourcells)%>%
  filter(sample_type %in% c("mass_ERWSI", "mass_FFPEmol", "mass_FFPEother", "mass_leftover"))%>%
  group_by(sample)%>%
  reframe(sum = sum(number_tumourcells))%>%
  mutate(xplode = if_else(sample == "IHC slide", 0.2, 0))



hc = highchart() %>%
  hc_chart(type = "pie") %>%
  hc_plotOptions(pie = list(
    dataLabels = list(
      enabled = TRUE,
      style = list(
        fontSize = "20px",  
        fontWeight = "bold"
      )
    )
  )) %>%
  hc_series(
    list(
      data = list(
        list(name = pie_df$sample[3], y = pie_df$sum[3], sliced = TRUE, selected = TRUE, color = "#E31A1CB7"),  # Explode
        list(name = pie_df$sample[1], y = pie_df$sum[1], color = sr_col),
        list(name = pie_df$sample[2], y = pie_df$sum[2], color = "#1F78B4FF"),
        list(name = pie_df$sample[4], y = pie_df$sum[4], color = rs_col)
      )
    )
  )
saveWidget(hc, file = file.path(OUT_DIR, "Figure1E_piechart.html"), selfcontained = TRUE)




