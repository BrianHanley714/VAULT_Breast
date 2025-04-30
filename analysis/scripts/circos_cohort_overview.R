# Generate a Circos overview plot of the VAULT breast cancer samples

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(geomtextpath)
library(tidyverse)

# PATHS -------------------------------------------------------------------


BASE = here::here()
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
OUT_DIR = file.path(BASE, "analysis", "figures")
CLINDATA = file.path(BASE, "data", "metadata", "clinical_data.txt")
ENDPOINT = file.path(BASE, "data", "metadata", "VAULT_endpoint_data.txt")
TUM_COUNTS = file.path(BASE, "data", "image_analysis", "tumour_cell_counts.tsv")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_features.R"))


# LOAD DATA ---------------------------------------------------------------
clinical_data = read.delim(CLINDATA)
endpoint_data = read.delim(ENDPOINT)
annotations = read.delim(TUM_COUNTS)

# get clinical data
clinical_data = left_join(clinical_data, endpoint_data, by = "Trial.ID")

# annotate recurrences
clinical_data$Trial.ID = if_else(clinical_data$recurrence_current_sample_is_recurence == "Yes", paste0("*", clinical_data$Trial.ID), clinical_data$Trial.ID)

# rename treatment classes for clarity
clinical_data$neo_adj_treatment_class = if_else(grepl("SACT", clinical_data$neo_adj_treatment_class), "Neoadj_other", clinical_data$neo_adj_treatment_class)
clinical_data$neo_adj_treatment_class = if_else(grepl("NET", clinical_data$neo_adj_treatment_class), "NET_only", clinical_data$neo_adj_treatment_class)
clinical_data$neo_adj_treatment_class = if_else(is.na(clinical_data$neo_adj_treatment_class), "Rx_naive", clinical_data$neo_adj_treatment_class)

# simplify names
clinical_data$Tumour_Weight = as.numeric(clinical_data$Tumour.leftover.tissue.weight..grams.)
clinical_data$Normal_Weight = as.numeric(clinical_data$Normal.leftover.tissue.weight..grams.)
clinical_data$Tumadj_Weight = as.numeric(clinical_data$Tumour.leftover.leftover.tissue.weight..grams.)

# order dataframe
clinical_data$Receptor.Subtype = factor(clinical_data$Receptor.Subtype, levels = c("HER2+", "ER+/HER2-","TN"))
clinical_data = clinical_data%>%
  arrange(Receptor.Subtype, Histological.Subtype, neo_adj_treatment_class, Tumour_Weight)

# get the AVERAGE volume of pathologist chosen block
median_wt_molecular_block = annotations%>%
  mutate(wt_block = vol_ERblock/(10^12))%>%
  pull(wt_block)%>%median()





#graphical parameters
out_width = 0.5
col_leftover = "#EFEDF5"
col_block = "#54278F"
aes_width = 0.05
ER_col = nice_cols[5] 
HER2_col = nice_cols[6]
TN_col = nice_cols[7]
NST_col = "#DECBE4"
Muc_col = "#FBB4AE"
Mix_col = "#B3CDE3"
ILC_col = "#CCEBC5"
vertical_linewidth = 1.5
break_width = 0.05
Naive_col = "light grey"
NET_col = "#9e9e9e"
SACT_col = "#36454F" 
inner_width = 0.1
relsize_tumadj_norm = 0.15


# add graphic data
clinical_data = clinical_data%>%
  mutate(inner = inner_width,
         prop_tumour_lo = Tumour_Weight/(Tumour_Weight+median_wt_molecular_block),#lo indicates leftover
         prop_tumour_block = median_wt_molecular_block/(Tumour_Weight+median_wt_molecular_block),
         
         prop_normal_lo = (Normal_Weight/(Normal_Weight+median_wt_molecular_block)*relsize_tumadj_norm),#lo indicates leftover
         prop_normal_block = (median_wt_molecular_block/(Normal_Weight+median_wt_molecular_block)*relsize_tumadj_norm),
         
         prop_ta_lo = (Tumadj_Weight/(Tumadj_Weight+median_wt_molecular_block)*relsize_tumadj_norm),#lo indicates leftover, ta indicates tumour adjacent
         prop_ta_block = (median_wt_molecular_block/(Tumadj_Weight+median_wt_molecular_block)*relsize_tumadj_norm),
         outer = out_width,
         break_w = break_width, 
         aes = aes_width
  )%>%
  filter(!is.na(prop_tumour_lo))%>%
  mutate(order = row_number())


clinical_data_long = clinical_data%>%
  pivot_longer(cols = c(starts_with("prop"), "outer","inner", "aes" , "break_w"), names_to = "sample_type", values_to = "Proportion")%>%
  dplyr::select(Trial.ID, Proportion, sample_type, order, Receptor.Subtype, Histological.Subtype, neo_adj_treatment_class)%>%
  mutate(sample_type = if_else(sample_type == "aes", paste0(sample_type, "_", Receptor.Subtype), sample_type))

clinical_data_long = bind_rows(clinical_data_long, 
                               clinical_data_long%>%
                                 filter(!duplicated(Trial.ID))%>%
                                 mutate(Proportion= aes_width, 
                                        sample_type = paste0("aes_", Histological.Subtype)),
                               clinical_data_long%>%
                                 filter(!duplicated(Trial.ID))%>%
                                 mutate(Proportion= aes_width, 
                                        sample_type = paste0("aes_", neo_adj_treatment_class))
)





# apply an order to the plot ----------------------------------------------
# note the ordering is backward here relative to the plot

clinical_data_long$sample_type = factor(clinical_data_long$sample_type, levels = c("outer",
                                                                                   "aes_ER+/HER2-",
                                                                                   "aes_HER2+",
                                                                                   "aes_TN",
                                                                                   "aes_IDC_NST",
                                                                                   "aes_MucinousCarcinoma",
                                                                                   "aes_MixedILC/IDC",
                                                                                   "aes_ILC",
                                                                                   "aes_NET_only", 
                                                                                   "aes_Rx_naive",
                                                                                   "aes_Neoadj_other",
                                                                                   "break_w",
                                                                                   "prop_normal_block",
                                                                                   "prop_normal_lo",
                                                                                   "prop_ta_block",
                                                                                   "prop_ta_lo",
                                                                                   "prop_tumour_block",
                                                                                   "prop_tumour_lo",
                                                                                   "inner"
))


#assign colours to plot
colours = c("prop_tumour_lo" = col_leftover, 
            "prop_tumour_block" = col_block, 
            "prop_normal_lo" = col_leftover, 
            "prop_normal_block" = col_block, 
            "prop_ta_lo" = col_leftover, 
            "prop_ta_block" = col_block,
            "inner" = "white",
            "outer" = "white",
            "break_w" = "white",
            "aes_ER+/HER2-" = ER_col,
            "aes_HER2+" = HER2_col,
            "aes_TN" = TN_col,
            "aes_IDC_NST" = NST_col,
            "aes_MucinousCarcinoma" = Muc_col,
            "aes_MixedILC/IDC" = Mix_col,
            "aes_ILC" = ILC_col,
            "aes_NET_only" = NET_col, 
            "aes_Rx_naive"= Naive_col,
            "aes_Neoadj_other" = SACT_col)

# get the plot radius
radius = c(clinical_data_long%>%
             group_by(Trial.ID)%>%reframe(sum = sum(Proportion))%>%pull(sum))[1]

# get the starts for plot segments
HER2_start = sum(clinical_data$Receptor.Subtype == "HER2+") + 0.5
TN_start = nrow(clinical_data) -sum(clinical_data$Receptor.Subtype == "TN") + 0.5




# Generate the plot -------------------------------------------------------
clinical_data_long%>%
  mutate(label = if_else(duplicated(Trial.ID), "", Trial.ID),
         angle = 270-((360/nrow(clinical_data))*order))%>%
  ggplot() +
  geom_bar(aes(x = as.factor(order), y = Proportion, fill = sample_type), stat = "identity", position = "stack", width = 1) +
  coord_curvedpolar(theta = "x") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_labelsegment(aes(x = 0.5, xend = HER2_start, y = radius - out_width, yend = radius - out_width), col = HER2_col, fontface = "bold", label = "HER2+")+
  geom_labelsegment(aes(x = nrow(clinical_data)+0.5, xend = TN_start, y = radius - out_width, yend = radius - out_width), fontface = "bold", col = TN_col, label = "TN")+
  geom_labelsegment(aes(x = HER2_start, xend = TN_start, y = radius - out_width, yend = radius - out_width), col = ER_col, fontface = "bold", label = "ER+ HER2-")+
  geom_labelhline(yintercept = 1+inner_width, label = "Tumour", fontface = "bold", size = 3.5)+
  geom_labelhline(yintercept = 1+relsize_tumadj_norm+inner_width, label = "Tumour Adjacent", fontface = "bold", size = 2)+
  geom_labelhline(yintercept = 1+relsize_tumadj_norm*2+inner_width, label = "Distant Normal",fontface = "bold", size = 2)+

  geom_vline(xintercept = 0.5, col ="white", size = vertical_linewidth)+
  geom_vline(xintercept = HER2_start, col ="white",size = vertical_linewidth)+
  geom_vline(xintercept = TN_start, col ="white",size = vertical_linewidth) +
  scale_fill_manual(values = colours)+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")+
  geom_text(aes(order, radius - out_width +0.15, label = label, angle = angle), fontface = "bold", size = 2.5)

ggsave(file.path(OUT_DIR, "Figure1B_circos_overview.png"))

# Generate the Legends (separate to plot) ---------------------------------
plot.new()
legend("bottomleft", legend = c("Invasive Ductal Carcinoma (NST)", "Invasive Lobular Carcinoma","Mixed Carcinoma", "Mucinous Carcinoma", "NET alone", "Neoadj other", "Treatment Naive"), 
       fill = c(NST_col,ILC_col, Mix_col, Muc_col, NET_col, SACT_col, Naive_col), 
       title = expression(bold("Colour Legend")))

legend("bottomright", legend = c("Leftover Tissue", "One \"Representative\" FFPE Block"), 
       fill = c(col_leftover, col_block), 
       title = expression(bold("Sample Type")))

