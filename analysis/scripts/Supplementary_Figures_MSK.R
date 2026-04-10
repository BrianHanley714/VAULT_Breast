
rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(qs)
library(openxlsx)
library(mutSignatures)
# PATHS -------------------------------------------------------------------

BASE = here::here()
BASE = "/Users/hanleyb/Documents/GitHub/VAULT_Breast"
OUT_DIR = file.path(BASE, "analysis", "figures")
SIGNATURES_PATH = file.path(BASE, "data", "variants", "signatures_all.rdata")
CASES_INCL_PATH = file.path(BASE, "data", "metadata", "cases_included.xlsx")
SIGNATURES_PREFILT_PATH = file.path(BASE, "data", "variants", "signatures_out_pretcgafilt.rdata")
SIGNATURES_INDEL_PATH = file.path(BASE, "data", "variants", "signatures_out_indels.rdata")
VAULT_RAWINDEL_PATH = file.path(BASE, "data", "variants", "raw_indels.rdata")
VAULT_RAWSNV_PATH = file.path(BASE, "data", "variants", "raw_snvs.rdata")
TCGA_RAWINDEL_PATH = file.path(BASE, "data", "variants", "TCGA_indels.rdata")
MSK_RAWINDEL_PATH = file.path(BASE, "data", "variants", "MSK_indels.rdata")
SNV_TYPE_VECTOR_PATH = file.path(BASE, "data", "variants", "type_vector_snv.rdata")
SNV_TYPE_VECTOR_PREFILT_PATH = file.path(BASE, "data", "variants", "type_vector_snv_prefilt.rdata")
MSK_SNV_PATH = file.path(BASE, "data", "variants", "MSK_data_out_anno_snvs.rdata")
TCGA_SNV_PATH = file.path(BASE, "data", "variants", "TCGA_data_out_anno_snvs.rdata")
NIKZAINAL_SIGS_PATH = file.path(BASE, "data", "metadata", "breast_signatures_signal_nikzainal.csv")
FFPESIG_PATH = file.path(BASE, "data", "metadata", "ffpesig.csv")
VAULT_ENDPOINT_PATH = file.path(BASE, "data", "metadata", "VAULT_endpoint_data.txt")
INDEL_TYPE_VECTOR_PATH = file.path(BASE, "data", "variants", "type_vector_indels.rdata")
PANEL_CASE_PATH = file.path(BASE, "data", "metadata", "panelcases.txt")
MATCHED_PT_CHAR_MBTCGA = file.path(BASE, "data","metadata", "matched_patients_characteristics.txt")
IDMAP = file.path(BASE, "data","metadata", "tumouridmap_MB.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "rot.lab.R"))

# LOAD DATA ---------------------------------------------------------------
#survival_tcga = read.delim(TCGA_CLIN_PATH) 
matched_patients_char = read.delim(MATCHED_PT_CHAR_MBTCGA)
tumour_IDs = read.delim(IDMAP)
ck_included_cases_vault = read.xlsx(CASES_INCL_PATH, colNames = F)[,1]
matched_patients = c(matched_patients_char$PATIENT_ID, tumour_IDs$sample[match(matched_patients_char$PATIENT_ID, tumour_IDs$metabricId)])
matched_patients = matched_patients[!is.na(matched_patients)]
#survival_vault = read.delim(VAULT_CLIN_PATH)
#maf_out = qread(VAULT_VARIANTS_PATH)
signatures_out = qread(SIGNATURES_PATH)
signatures_out_pretcgafilt = qread(SIGNATURES_PREFILT_PATH)
signatures_out_indels = qread(SIGNATURES_INDEL_PATH)
mat_indel = qread(VAULT_RAWINDEL_PATH)
TCGA_indel = qread(TCGA_RAWINDEL_PATH)
MSK_indel = qread(MSK_RAWINDEL_PATH)
mat = qread(VAULT_RAWSNV_PATH)
type_vector = qread(SNV_TYPE_VECTOR_PATH)
type_vector_prefilt = qread(SNV_TYPE_VECTOR_PREFILT_PATH)
TCGA_data_out_anno = qread(TCGA_SNV_PATH)
MSK_data_out_anno = qread(MSK_SNV_PATH)
nikzainal_sigs = read.csv(NIKZAINAL_SIGS_PATH,check.names = F)
ffpe_sig = read.csv(FFPESIG_PATH)
time_in_formalin = read.delim(VAULT_ENDPOINT_PATH)
type_vector_indel = qread(INDEL_TYPE_VECTOR_PATH)
panel_cases = read.delim(PANEL_CASE_PATH)[,1]


# FUNCTIONS ---------------------------------------------------------------

# options -----------------------------------------------------------------
basesize = 20
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"


ffpe_samples = mat$SampleID%>%unique%>%grep("_B|_FT", ., value = T)%>%grep(paste(substr(ck_included_cases_vault, 1, 5), collapse = "|"), ., value = T)

included_cases = c(ffpe_samples, ck_included_cases_vault)[!c(ffpe_samples, ck_included_cases_vault)%in% panel_cases]
mat = mat%>%
  dplyr::filter(SampleID %in% included_cases)

mat_filtered_data = mat%>%
  dplyr::filter(
    (
      vaf_tumour>= 0.05
      &alt_count_tumour > 4
      &msec_artefact ==FALSE
      & (is.na(p_ADfit) | p_ADfit<0.05)
      & ROQ>20
      &!(grepl("clustered_events", FILTER) & vaf_tumour <0.1)
    )
    |
      (vaf_tumour> 0.01 &FILTER == "PASS"
       & msec_artefact ==FALSE
       &alt_count_tumour > 4
       # #     & (is.na(p_ADfit) | p_ADfit<0.05)
       & ROQ>20
       & variant_in_deepsomatic ==T
      )
  )

mut.counts = countMutTypes(mutTable = mat_filtered_data,
                           mutType_colName = "mutType",
                           sample_colName = "SampleID")

mat_filtered = mutSignatures::as.data.frame(mut.counts)

mat_prefiltered = mat_filtered
mat_filtered = mat_filtered[, colSums(mat_filtered) >=50]
#sum(is.na(mat$altmat_filtered#sum(is.na(mat$alt_count))



names(mat_filtered) = paste("filtered", names(mat_filtered), sep = "_")





# create distance cosine umap for snvs ------------------------------------


distance = proxy::dist(t(signatures_out), method = "cosine")

umap = uwot::umap(distance)

umap = as.data.frame(umap)

colnames(umap) = c("UMAP1", "UMAP2")


names(signatures_out) = paste(names(signatures_out), type_vector, sep = "#")

mutcount_matrix = as.data.frame(cbind(counts = as.numeric(colSums(signatures_out_pretcgafilt)), samples = names(signatures_out_pretcgafilt)))
variants = signatures_out
names(variants) = paste(type_vector, names(signatures_out), sep = "#")

colours = c("FFIX_raw" = "#E31A1Cb7" ,
            "FFPE_raw" = "#FF7F00ff", 
            "Breast_TCGA" = "#1F78B4db" ,
            "FFIX_filt" = "#6A3D9Adb",
            "FFPE_filt" = "#33A02Cb7",
            "Breast_MSK" ="#FB9A99db"
)

levels_titles = c("FFIX_raw",
                  "FFPE_raw", 
                  "FFIX_filt",
                  "FFPE_filt",
                  "Breast_TCGA",
                  "Breast_MSK"
)

names(levels_titles) = c("Raw", "Raw",rep("Filtered", 4))


# Output the plot ---------------------------------------------------------
row_order_mat = variants%>%  
  mutate(variants = substr(rownames(.), 3, 5))%>%
  pivot_longer(c(everything(), -variants))%>%
  mutate(type = sub("#.*", "", name))%>%
  group_by(type, name, variants)%>%
  reframe(count = sum(value))%>%
  group_by(name)%>%
  mutate(sum = sum(count),
         Prop = count/sum)%>%
  dplyr::select(-c(sum, count, type))%>%
  pivot_wider(names_from = name, values_from = Prop)%>%
  dplyr::select(-c(variants))%>%
  as.data.frame()%>%as.matrix()%>%t()

hc = hclust(dist(row_order_mat), method = "ward.D2")
row_order <- rownames(row_order_mat)[hc$order]


# generate snv plots ------------------------------------------------------
counts_snvs_plot = mutcount_matrix%>%
  mutate(study = substr(samples, 1, 1), 
         counts = as.numeric(counts),
         col = type_vector_prefilt,
         col = sub("CK_sort", "FFIX", col),
         col = factor(col, levels = levels_titles),
         filtering = names(levels_titles)[match(col, levels_titles)],
         filtering = factor(filtering, levels = c("Raw", "Filtered")))%>%
  group_by(col)%>%
  mutate(count = n(),
         count = if_else(duplicated(paste0(col, count)), NA, count))%>%
  ggplot(aes(col, counts))+
  geom_violin(width = 0.8, trim = F, scale = "width", aes(fill = col))+
  
  geom_jitter(width = 0.1, alpha = 0.2)+
  geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
  geom_label(aes(col, 5, label = count))+
  ylab("SNV counts")+
  scale_y_log10()+
  theme_classic()+
  theme(axis.title.x = element_blank())+
  rot.lab()+
  facet_grid(~filtering, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = colours)+
  stat_compare_means(
    data = \(d) d[d$filtering == "Filtered", ], label = "p.signif",
    comparisons =list(c("FFIX_filt", "Breast_TCGA"),
                      c("FFPE_filt", "Breast_TCGA")),
    label.y= 4.5
    )+
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )


umap_snv_all = umap%>%
  mutate(type = substr(rownames(umap), 7, 8))%>%
  mutate(colt = type_vector, 
         colt = sub("CK_sort", "FFIX", colt))%>%
  dplyr::filter(!grepl("filter", rownames(.)))%>%
  ggplot(aes(UMAP1, UMAP2))+
  
  #  geom_hdr(aes(fill = colt), alpha = 1, method = "kde", n = 200) +
  #stat_density_2d(geom = "polygon", aes(fill = colt, alpha = ..level..), color = NA, bins = 100) +
  scale_fill_brewer(palette = "Set2")+
  geom_point(aes(col = colt))+
  theme_classic(base_size = basesize)+
  scale_color_manual(values = colours)+
  guides(alpha = "none",
         col = guide_legend(title = NULL, nrow = 1))+
  theme(legend.position = "none")

umap_snv_shaded = umap%>%
  mutate(type = substr(rownames(umap), 7, 8),
         alpha = if_else(grepl("filter", rownames(.)), 1, .2))%>%
  mutate(col = type_vector,
         col = sub("CK_sort", "FFIX", col),
         col = factor(col, levels = levels_titles))%>%
  ggplot(aes(UMAP1, UMAP2))+
  geom_point(aes(col = col, alpha = alpha))+
  theme_classic(base_size = basesize)+
  guides(alpha = "none",
         col = guide_legend(title = NULL, nrow = 1))+
  scale_color_manual(values = colours)+
  
  theme(legend.position = "top")



snv_type_plot = variants%>%  
  mutate(variants = substr(rownames(.), 3, 5))%>%
  pivot_longer(c(everything(), -variants))%>%
  mutate(type = sub("#.*", "", name))%>%
  group_by(type, name, variants)%>%
  reframe(count = sum(value))%>%
  group_by(name)%>%
  mutate(sum = sum(count),
         Prop = count/sum,
         name = factor(name, levels = row_order),
         type = sub("CK_sort", "FFIX", type),
         type = factor(type, levels = levels_titles),
         filtering = names(levels_titles)[match(type, levels_titles)],
         filtering = factor(filtering, levels = c("Raw", "Filtered"))
  )%>%
  ggplot(aes(name, Prop, fill = variants))+
  geom_col(width = 1)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  #rot.lab()+
  scale_fill_brewer(palette = "Set2")+
  guides(fill = guide_legend(nrow =1, title = NULL))+
  #facet_grid(~type+filtering, scales = "free_x")
  facet_nested_wrap(~filtering+ type,  scales = "free_x",nrow = 1
                    #space = "free_x",
                    #labeller = labeller(col = filtering)
  )+
  ylab("Proportion of total SNVs")


snv_vaf_plot = mat_filtered_data%>%
  mutate(type = if_else(grepl("ET", SampleID), "FFIX_filt", "FFPE_filt"),
         mut = substr(mutType, 3, 5),
         VAF = vaf_tumour
  )%>%bind_rows(.,
                mat%>%
                  mutate(type = if_else(grepl("ET", SampleID), "FFIX_raw", "FFPE_raw"),
                         mut = substr(mutType, 3, 5),
                         VAF = vaf_tumour
                  )
  )%>%
  dplyr::select(type, mut, VAF)%>%
  bind_rows(.,
            TCGA_data_out_anno%>%
              mutate(VAF = t_alt_count/t_depth,
                     mut = substr(mutType, 3, 5),
                     type = "Breast_TCGA"
              ),
            MSK_data_out_anno%>%
              mutate(VAF = Mutant.allele.fraction,
                     mut = substr(mutType, 3, 5),
                     type = "Breast_MSK"
              )
  )%>%
  #dplyr::filter(type == "FFPE")%>%
  ggplot(aes(VAF, col = type, fill = type))+
  stat_density(geom = "area", position = "identity", alpha = 0.2, aes(y = ..scaled..))+
  theme_classic()+
  theme(legend.position = "bottom")+
  #scale_x_log10()+
  facet_wrap(~mut, ncol = 1)+
  #ggtitle("FFPE - Filtered")+
  guides(col = guide_legend(nrow =1))+
  theme(legend.position = "none",
        legend.justification = c("right", "top"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
  )+
  scale_color_manual(values = c(colours))+
  scale_fill_manual(values = c(colours))+
  guides(col = guide_legend(ncol = 1,title = NULL))



# Indel plots -------------------------------------------------------------


distance_indel = proxy::dist(t(signatures_out_indels), method = "cosine")

umap_indel = uwot::umap(distance_indel)

umap_indel = as.data.frame(umap_indel)

colnames(umap_indel) = c("UMAP1", "UMAP2")

# Get counts for the cohort -----------------------------------------------
mutcount_matrix_indel = as.data.frame(cbind(counts = as.numeric(colSums(signatures_out_indels)), samples = names(signatures_out_indels)))

mat_indel_filtered = mat_indel%>%
  #mutate(t_VAF = as.numeric(alt_count)/as.numeric(alt_count+ref_count))%>%
  dplyr::filter(
    
    (vaf_tumour>0.1
     &alt_count_tumour >4
     &!grepl("bp_insertion$", muttype) 
     &(is.na(p_ADfit)
       | p_ADfit<0.05
     ) 
     &msec_artefact == FALSE
     &!(grepl("clustered_events", FILTER) & vaf_tumour <0.3)
     #&(alt_count_tumour +alt_count_normal >30)
     #& FILTER == "PASS"
    )
    |
      (variant_in_deepsom == T 
       &vaf_tumour>0.1
       &alt_count_tumour >4
       &(is.na(p_ADfit)| p_ADfit<0.05) 
       &msec_artefact == FALSE
       #& FILTER == "PASS"
      )
  )


# arrange signatures plot -------------------------------------------------

order_facets = sub("_[^_]+$", "", rownames(signatures_out_indels))
order_facets = order_facets[!duplicated(order_facets)]

type_selected = 	"TCGA"

labels = substr(order_facets, 1, 1)
names(labels) = order_facets 

group = c("Del-C", "Del-T",
          "Ins-C", "Ins-T",
          "Del-Repeats", "Del-Repeats",
          "Del-Repeats", "Del-Repeats",
          "Ins-Repeats", "Ins-Repeats",
          "Ins-Repeats", "Ins-Repeats",                  
          "Del-MH", "Del-MH", 
          "Del-MH", "Del-MH")
xaxis_values = sub("^.*_", "", rownames(signatures_out))
names(xaxis_values) = rownames(signatures_out)


colours_indels = c("FFIX_raw" = "#E31A1Cb7" ,
            "FFPE_raw" = "#FF7F00ff", 
            "TCGA" = "#1F78B4db" ,
            "FFIX_filtered" = "#6A3D9Adb",
            "FFPE_filtered" = "#33A02Cb7",
            "MSK" = "#FB9A99db"
)

# signatures_out%>%
#   mutate(label = rownames(.))

# Output Plot -------------------------------------------------------------

levels_titles_indels = c("FFIX_raw",
                  "FFPE_raw", 
                  "FFIX_filtered",
                  "FFPE_filtered",
                  "TCGA",
                  "MSK"
)


names(levels_titles_indels) = c("Raw", "Raw",rep("Filtered", 4))
row_order_mat = t(signatures_out_indels)%>%
  as.data.frame()%>%
  mutate(type = type_vector_indel,
         samples = rownames(.)
  )%>%
  #dplyr::filter(type == type_selected)%>%
  pivot_longer(c(everything(), -c(type, samples)))%>%
  mutate(col = sub("_[^_]+$", "", name),
         col = factor(col, levels = order_facets),
         grouping = group[match(col, order_facets)],
         grouping = factor(grouping, levels = group[!duplicated(group)]),
         name = factor(name, levels = rownames(signatures_out)),
         samples = paste0(type, samples)
  )%>%dplyr::select(-c(col, name))%>%
  group_by(type, samples, grouping)%>%
  dplyr::mutate(value = if_else(grepl("raw", type), 0, value))%>%
  reframe(sum = sum(value),
          #prop = value/sum
  )%>%  group_by(samples)%>%
  mutate(sum_total = sum(sum),
         prop = sum/sum_total)%>%
  dplyr::select(-c(sum, sum_total, type))%>%
  pivot_wider(names_from = samples, values_from = prop)%>%
  dplyr::select(-c(grouping))%>%
  as.data.frame()%>%as.matrix()%>%t()

distance_mat = dist(row_order_mat[!grepl("raw", rownames(row_order_mat)),])

distance_mat[grep("FFPE_filt", rownames(distance_mat)), grep("FFPE_r", rownames(distance_mat))]

hc = hclust(distance_mat, method = "complete")


row_order <- (rownames(row_order_mat)[!grepl("raw", rownames(row_order_mat))])[hc$order]
row_order = c(row_order, rownames(row_order_mat)[grepl("raw", rownames(row_order_mat))])

indel_count_plot = mutcount_matrix_indel%>%
  mutate(study = substr(samples, 1, 1), 
         counts = as.numeric(counts),
         col = type_vector_indel,
         col = sub("RepSeq", "FFIX", col),
         col = factor(col, levels = levels_titles_indels),
         filtering = names(levels_titles)[match(col, levels_titles_indels)],
         filtering = factor(filtering, levels = c("Raw", "Filtered"))
  )%>%
  group_by(col)%>%
  mutate(count = n(),
         count = if_else(duplicated(paste0(col, count)), NA, count))%>%
  ggplot(aes(col, counts))+
  geom_violin(aes(fill = col), width = 0.8, trim = F, scale = "width")+
  
  geom_jitter(width = 0.1, alpha = 0.2)+
  geom_boxplot(width = 0.4, outlier.alpha = 0, varwidth = F)+
  geom_label(aes(col, 0.5, label = count))+
  scale_y_log10()+theme_classic()+
  scale_fill_manual(values = colours_indels)+
  stat_compare_means(data = \(d) d[d$filtering == "Filtered", ], label = "p.signif",
                     comparisons = list(c("FFIX_filtered", "TCGA"),
                                          c("FFPE_filtered", "TCGA")),
                     label.y = 4)+
  facet_grid(~filtering, scales = "free_x", space = "free_x")+
  theme(axis.title.x = element_blank())+
  rot.lab()+
  ylab("INDEL Counts")+
  theme(legend.position = "none")


umap_indel_plot = umap_indel%>%
  mutate(col = type_vector_indel)%>%
  dplyr::filter(!grepl("filter", col))%>%
  ggplot(aes(UMAP1, UMAP2, col = col))+
  geom_point()+
  scale_color_manual(values = colours)+
  theme_classic()+
  theme(legend.position = "none")

umap_indel_plot_shaded = umap_indel%>%
  mutate(col = type_vector_indel,
         alpha = if_else(grepl("filter", col), 1, 0.2))%>%
  ggplot(aes(UMAP1, UMAP2, col = col))+
  
  #  geom_hdr(aes(fill = colt), alpha = 1, method = "kde", n = 200) +
  #stat_density_2d(geom = "polygon", aes(fill = colt, alpha = ..level..), color = NA, bins = 100) +
  scale_color_manual(values = colours)+
  geom_point(aes(alpha = alpha))+
  theme_classic()+
  theme(legend.position = "bottom")+
  guides(alpha = "none",
         col = guide_legend(title = NULL, nrow = 1))+
  theme(legend.position = "top")

indel_plot_vaf = mat_indel_filtered%>%
  mutate(type = if_else(grepl("ET", SampleID), "FFIX_filtered", "FFPE_filtered"),
         VAF = vaf_tumour)%>%dplyr::select(type, VAF, muttype)%>%
  bind_rows(.,
            mat_indel%>%mutate(type = if_else(grepl("ET", SampleID), "FFIX_raw", "FFPE_raw"),
                               VAF = vaf_tumour),
            TCGA_indel%>%mutate(type = "TCGA",
                                VAF = t_alt_count/t_depth),
            MSK_indel%>%mutate(type = "MSK",
                               VAF = Mutant.allele.fraction)
  )%>%dplyr::select(VAF, type, muttype)%>%
  mutate(muttype = if_else(as.numeric(sub("bp.*", "", muttype))>=5 &!is.na(as.numeric(sub("bp.*", "", muttype))), sub(".*bp", "5+bp", muttype), muttype),
         mut_group = group[match(muttype, names(labels))],
         mut_group = factor(mut_group, levels = group[!duplicated(group)])
         
  )%>%
  ggplot(aes(VAF, col = type, fill = type))+
  stat_density(geom = "area", position = "identity", alpha = 0.2, aes(y = ..scaled..))+
  theme_classic()+
  theme(legend.position = "bottom")+
  facet_wrap(~mut_group, ncol = 1)+
  guides(col = guide_legend(nrow =1))+
  theme(legend.position = "none",
        legend.justification = c("right", "top"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_color_manual(values = colours_indels)+
  scale_fill_manual(values = colours_indels)+
  guides(col = guide_legend(ncol = 1,title = NULL))

indel_types_plot = t(signatures_out_indels)%>%
  as.data.frame()%>%
  mutate(type = type_vector_indel,
         samples = rownames(.)
  )%>%
  #dplyr::filter(type == type_selected)%>%
  pivot_longer(c(everything(), -c(type, samples)))%>%
  mutate(col = sub("_[^_]+$", "", name),
         col = factor(col, levels = order_facets),
         grouping = group[match(col, order_facets)],
         grouping = factor(grouping, levels = group[!duplicated(group)]),
         name = factor(name, levels = rownames(signatures_out_indels)),
         samples = paste0(type, samples)
  )%>%dplyr::select(-c(col, name))%>%
  group_by(type, samples, grouping)%>%
  reframe(sum = sum(value),
          #prop = value/sum
  )%>%  group_by(samples)%>%
  mutate(sum_total = sum(sum),
         prop = sum/sum_total,
         samples = factor(samples, levels = row_order),
         type = sub("RepSeq", "FFIX", type),
         )%>%
  mutate(
         type = factor(type, levels = levels_titles_indels),
         filtering = names(levels_titles)[match(type, levels_titles_indels)],
         filtering = factor(filtering, levels = c("Raw", "Filtered")))%>%
  ggplot(aes(samples, prop, fill = grouping))+
  geom_col(width = 1)+
  #facet_wrap(~type, scales = "free_x", nrow = 1)+
  facet_nested_wrap(~filtering+ type,  scales = "free_x",nrow = 1
                    #space = "free_x",
                    #labeller = labeller(col = filtering)
  )+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.line = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill = guide_legend(nrow =1, title = NULL))+
  ylab("Proportion of total INDELs")



# WRITE THE PLOTS ---------------------------------------------------------





plot_grid(counts_snvs_plot, indel_count_plot, ncol = 1, rel_heights = c(1, 1.4))

ggsave(file.path(OUT_DIR, "SupplementaryFigure_snv_indel_counts_by_filtering.png"))
umap_snv_all
ggsave(file.path(OUT_DIR, "SupplementaryFigure_umap_by_filtering.png"))
umap_snv_shaded
ggsave(file.path(OUT_DIR, "SupplementaryFigure_umap_by_filtering_shaded.png"))

plot_grid(snv_type_plot, indel_types_plot, ncol = 1)
ggsave(file.path(OUT_DIR, "SupplementaryFigure_variant_type_by_filtering.png"))

plot_grid(snv_vaf_plot, indel_plot_vaf)
ggsave(file.path(OUT_DIR, "SupplementaryFigure_by_filtering.png"))

time_in_formalin_plot
ggsave(file.path(OUT_DIR, "Supplementaryfigure_time_in_formalin.png"))



umap_indel_plot
ggsave(file.path(OUT_DIR, "SupplementaryFigure_indel_umap.png"))
umap_indel_plot_shaded
ggsave(file.path(OUT_DIR, "SupplementaryFigure_indel_umap_shaded.png"))

