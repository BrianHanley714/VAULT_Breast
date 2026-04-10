# This is the script for reproducing analysis and figures in in Figure 2 of the manuscript
rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------
library(qs)
library(openxlsx)
library(mutSignatures)
library(ggh4x)
library(ComplexHeatmap)
# PATHS -------------------------------------------------------------------

BASE = here::here()
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
get_exposure_counts = function(dataframe){
  mut.counts_FFPE = countMutTypes(mutTable = dataframe,
                                  mutType_colName = "mutType",
                                  sample_colName = "SampleID")
  
  mut.counts_FFPE = mutSignatures::as.data.frame(mut.counts_FFPE)
  
  normalised_ffpe = sweep((mut.counts_FFPE), 2, colSums(mut.counts_FFPE, na.rm = TRUE), FUN = "/")
  
  
  
  rownames(nikzainal_sigs) = nikzainal_sigs$Signature
  nikzainal_sigs = nikzainal_sigs[,3:length(colnames(nikzainal_sigs))]
  nikzainal_sigs = t(nikzainal_sigs)
  nikzainal_sigs = nikzainal_sigs[,colnames(nikzainal_sigs) %in% paste0("SBS", c(1, 5, 2, 18, 8, 13, 3, 127))] # these are the sbs signatures seen in >10% of breast cancers according to https://signal.mutationalsignatures.com/
  cosmic_signatures = mutSignatures::getCosmicSignatures()
  
  cosmic_signatures_df = mutSignatures::as.data.frame(cosmic_signatures)
  cosmic_signatures_df = cosmic_signatures_df[,!colnames(cosmic_signatures_df) %in% c("COSMIC.1", "COSMIC.6")]
  
  nikzainal_sigs = nikzainal_sigs[match(cosmic_signatures@mutTypes$mutTypes, rownames(nikzainal_sigs)),]
  rownames(nikzainal_sigs) ==rownames( cosmic_signatures_df )
  
  
  cosmic_signatures_df = mutSignatures::as.data.frame(cosmic_signatures)
  cosmic_signatures_df = cosmic_signatures_df[,!colnames(cosmic_signatures_df) %in% c("COSMIC.1", "COSMIC.6")]
  cosmic_signatures_df_ffpesig = cbind(cosmic_signatures_df, ffpe_sig_rep = ffpe_sig_match[,2])
  
  nikzainal_sigs_ffpe = as.data.frame(cbind(nikzainal_sigs, ffpe_sig_rep = ffpe_sig_match[,2]))
  out = matchSignatures(mutSignatures::as.mutation.signatures(normalised_ffpe), mutSignatures::as.mutation.signatures(nikzainal_sigs_ffpe))
  hm_matrix = out$distanceDataFrame%>%
    dplyr::select(-dist)%>%
    pivot_wider(values_from = simil, names_from = refSign)%>%as.data.frame()
  rownames(hm_matrix) = hm_matrix$newSign
  hm_matrix = hm_matrix[,-1]

  return(hm_matrix)
}

get_cosffpesig_by_alt_count = function(min, max, variable){
  
  if(variable == "alt_count"){
    ffpesig_filtered = mat%>%
      dplyr::filter(
        msec_artefact ==FALSE
        , 
        alt_count_tumour >= min &alt_count_tumour<=max)
  }
  
  if(variable == "VAF"){
    ffpesig_filtered = mat%>%
      dplyr::filter(
        msec_artefact ==FALSE
        , 
        vaf_tumour >= min &vaf_tumour<=max
      )
  }
  
  mut.counts_FFPE = countMutTypes(mutTable = ffpesig_filtered,
                                  mutType_colName = "mutType",
                                  sample_colName = "SampleID")
  
  mut.counts_FFPE = mutSignatures::as.data.frame(mut.counts_FFPE)
  
  normalised_ffpe = sweep((mut.counts_FFPE), 2, colSums(mut.counts_FFPE, na.rm = TRUE), FUN = "/")
  
  
  matched_ffpe = matchSignatures(mutSignatures::as.mutation.signatures(normalised_ffpe), as.mutation.signatures(ffpe_sig_match[,2:3]))
  
  out_df = matched_ffpe$distanceDataFrame%>%
    dplyr::filter(refSign == "FFPE_Repaired_sig")%>%
    mutate(min = min,
           max = max,
           mut_counts = colSums(mut.counts_FFPE))
  
  
  return(out_df)
}
col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#1c9099", "#f0f0f0", "#762a83"))
plot_hm = function(mat, column_title){Heatmap(mat, name = "COS", 
                                              heatmap_legend_param = gpar(legend_height = unit(10, "cm")),
                                              col = col_fun, show_row_names = F,
                                              show_row_dend = F, show_column_dend = F, 
                                              row_names_side = "left", 
                                              row_names_gp = gpar(fontsize = 5), 
                                              column_title = column_title)}


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




# Filter out all MSK and un-matched TCGA cases ----------------------------


signatures_out = signatures_out[, grep("_MSK", type_vector, invert = T)]
type_vector = type_vector[grep("_MSK", type_vector, invert = T)]
signatures_out_pretcgafilt = signatures_out_pretcgafilt[, grep("_MSK", type_vector_prefilt, invert = T)]
type_vector_prefilt = type_vector_prefilt[grep("_MSK", type_vector_prefilt, invert = T)]
signatures_out_indels = signatures_out_indels [,grep("MSK", invert = T, type_vector_indel)]
type_vector_indel = type_vector_indel[grep("MSK", invert = T, type_vector_indel)]

matched_pt_grep = paste(c(matched_patients, "_HF", "^HF"), collapse = "|")

type_vector = type_vector[grep(matched_pt_grep, colnames(signatures_out), invert = F)]
signatures_out = signatures_out[, grep(matched_pt_grep, colnames(signatures_out), invert = F)]

type_vector_prefilt = type_vector_prefilt[grep(matched_pt_grep, colnames(signatures_out_pretcgafilt), invert = F)]
signatures_out_pretcgafilt = signatures_out_pretcgafilt[,grep(matched_pt_grep, colnames(signatures_out_pretcgafilt), invert = F)]

type_vector_indel = type_vector_indel[grep(matched_pt_grep, colnames(signatures_out_indels), invert = F)]
signatures_out_indels = signatures_out_indels [,grep(matched_pt_grep, colnames(signatures_out_indels), invert = F)]


# create distance cosine umap for snvs ------------------------------------


distance = proxy::dist(t(signatures_out), method = "cosine")

umap = uwot::umap(distance)

umap = as.data.frame(umap)

colnames(umap) = c("UMAP1", "UMAP2")


names(signatures_out) = paste(names(signatures_out), type_vector, sep = "#")

mutcount_matrix = as.data.frame(cbind(counts = as.numeric(colSums(signatures_out_pretcgafilt)), samples = names(signatures_out_pretcgafilt)))
variants = signatures_out
names(variants) = paste(type_vector, names(signatures_out), sep = "#")

colours = c("RepSeq_raw" = "#E31A1Cb7" ,
            "FFPE_raw" = "#FF7F00ff", 
            "Breast_TCGA" = "#1F78B4db" ,
            "RepSeq_filt" = "#6A3D9Adb",
            "FFPE_filt" = "#33A02Cb7",
            "Breast_MSK" ="#FB9A99db"
)

levels_titles = c("RepSeq_raw",
                  "FFPE_raw", 
                  "RepSeq_filt",
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
         col = sub("CK_sort", "RepSeq", col),
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
    comparisons =list(c("RepSeq_filt", "Breast_TCGA"),
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
         colt = sub("CK_sort", "RepSeq", colt))%>%
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
         col = sub("CK_sort", "RepSeq", col),
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
         type = sub("CK_sort", "RepSeq", type),
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
  mutate(type = if_else(grepl("ET", SampleID), "RepSeq_filt", "FFPE_filt"),
         mut = substr(mutType, 3, 5),
         VAF = vaf_tumour
  )%>%bind_rows(.,
                mat%>%
                  mutate(type = if_else(grepl("ET", SampleID), "RepSeq_raw", "FFPE_raw"),
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
            # MSK_data_out_anno%>%
            #   mutate(VAF = Mutant.allele.fraction,
            #          mut = substr(mutType, 3, 5),
            #          type = "Breast_MSK"
            #   )
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

mean_depths = mat%>%
  mutate(depth = alt_count_tumour + ref_count_tumour)%>%
  group_by(SampleID)%>%
  reframe(mean_depth = median(depth))

time_in_formalin_plot = mutcount_matrix%>%
  mutate(type = type_vector_prefilt)%>%
  dplyr::filter(grepl("CK_sort_raw", type))%>%
  
  mutate(time_in_formalin = time_in_formalin$Time..in.formalin..days.[match(substr(samples, 1, 5), time_in_formalin$Trial.ID)],
         mean_depth = mean_depths$mean_depth[match(sub("#.*", "", samples), mean_depths$SampleID)],
         counts = as.numeric(counts),
         norm_counts = counts/(mean_depth))%>%
  ggplot(aes(time_in_formalin, norm_counts, label = substr(samples, 1, 5)))+
  geom_point()+
  #  geom_text()+
  geom_smooth(col = "black", method = "lm")+
  stat_cor()+
  theme_classic(base_size = basesize)+
  ylab("Normalised Counts (Unfiltered Count / mean depth)")+
  xlab("Time in Formalin")



mut_type = paste0(substr(ffpe_sig$MutationType, 5, 5), "[", substr(ffpe_sig$MutationType, 1, 3), "]", substr(ffpe_sig$MutationType, 7, 7))
signatures = signatures_out[,grepl("raw", type_vector)]
ffpe_sig_match = ffpe_sig[match(rownames(signatures), mut_type),]
rownames(ffpe_sig_match) = rownames(signatures)


mutsigs_ffpe = as.mutation.signatures(ffpe_sig_match[2:3])

#mutSignatures::mutsigs_ffpe = rownames(signatures)
normalised_signatures = sweep((signatures), 2, colSums(signatures, na.rm = TRUE), FUN = "/")


matched_ffpe = matchSignatures(mutSignatures::as.mutation.signatures(normalised_signatures), as.mutation.signatures(ffpe_sig_match[,2:3]))




mat_out = rbind(
  get_cosffpesig_by_alt_count(0,0.01, variable = "VAF"),
  get_cosffpesig_by_alt_count(0.01,0.02, variable = "VAF"),
  get_cosffpesig_by_alt_count(0.02,0.03, variable = "VAF"),
  get_cosffpesig_by_alt_count(0.03,0.04, variable = "VAF"),
  get_cosffpesig_by_alt_count(0.04,0.05, variable = "VAF"),
  get_cosffpesig_by_alt_count(0.05,1, variable = "VAF")
)

FFPE_sig_plot = mat_out%>%
  mutate(type = if_else(newSign %in% ck_included_cases_vault, "CK_sort", "FFPE"),
         max = factor(max, levels = c(0.01,0.02, 0.03,0.04,0.05,1)))%>%
  group_by(max)%>%
  mutate(mean_count = round(mean(mut_counts), 0),
         mean_count = if_else(duplicated(mean_count), NA, mean_count))%>%
  
  ggplot(aes(max, simil))+
  #geom_point()+
  #geom_line()+
  geom_violin(aes(fill = max))+
  geom_boxplot(outlier.alpha = 0, width = 0.2)+
  geom_text(aes(max, 1.05, label = mean_count))+
  theme_bw(base_size = basesize)+
  theme(legend.position = "none")+
  scale_color_brewer(palette = "Set2")+
  xlab("VAF")+
  scale_fill_brewer()+
  ylab("Cosine Similarity to FFPESig (repaired)")




unfiltered_exposures = get_exposure_counts(mat)
filtered_exposures = get_exposure_counts(mat_filtered_data)
colnames(filtered_exposures)[colnames(filtered_exposures) == "ffpe_sig_rep"] = "FFPEsig"
colnames(unfiltered_exposures)[colnames(unfiltered_exposures) == "ffpe_sig_rep"] = "FFPEsig"

heatmap_plot_signatures = ComplexHeatmap::draw(plot_hm(unfiltered_exposures, column_title =  "Pre Filtering")%v% plot_hm(filtered_exposures, column_title =  "Post Filtering"))


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


colours_indels = c("RepSeq_raw" = "#E31A1Cb7" ,
            "FFPE_raw" = "#FF7F00ff", 
            "TCGA" = "#1F78B4db" ,
            "RepSeq_filtered" = "#6A3D9Adb",
            "FFPE_filtered" = "#33A02Cb7",
            "MSK" = "#FB9A99db"
)

# signatures_out%>%
#   mutate(label = rownames(.))

# Output Plot -------------------------------------------------------------

levels_titles_indels = c("RepSeq_raw",
                  "FFPE_raw", 
                  "RepSeq_filtered",
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
         col = sub("RepSeq", "RepSeq", col),
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
                     comparisons = list(c("RepSeq_filtered", "TCGA"),
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
  mutate(type = if_else(grepl("ET", SampleID), "RepSeq_filtered", "FFPE_filtered"),
         VAF = vaf_tumour)%>%dplyr::select(type, VAF, muttype)%>%
  bind_rows(.,
            mat_indel%>%mutate(type = if_else(grepl("ET", SampleID), "RepSeq_raw", "FFPE_raw"),
                               VAF = vaf_tumour),
            TCGA_indel%>%mutate(type = "TCGA",
                                VAF = t_alt_count/t_depth),
            # MSK_indel%>%mutate(type = "MSK",
            #                    VAF = Mutant.allele.fraction)
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
         type = sub("RepSeq", "RepSeq", type),
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

ggsave(file.path(OUT_DIR, "2B_C_snv_indel_counts_by_filtering.png"))
umap_snv_all
ggsave(file.path(OUT_DIR, "2F_umap_by_filtering.png"))
umap_snv_shaded
ggsave(file.path(OUT_DIR, "2F_umap_by_filtering_shaded.png"))

plot_grid(snv_type_plot, indel_types_plot, ncol = 1)
ggsave(file.path(OUT_DIR, "2D_E_variant_type_by_filtering.png"))

plot_grid(snv_vaf_plot, indel_plot_vaf)
ggsave(file.path(OUT_DIR, "2I_VAF_by_filtering.png"))

time_in_formalin_plot
ggsave(file.path(OUT_DIR, "Supplementaryfigure_time_in_formalin.png"))

FFPE_sig_plot
ggsave(file.path(OUT_DIR, "2H_FFPEsig_by_VAF.png"))

heatmap_plot_signatures
ggsave(file.path(OUT_DIR, "2G_FFPEsig_by_VAF.png"))


umap_indel_plot
ggsave(file.path(OUT_DIR, "SupplementaryFigure_indel_umap.png"))
umap_indel_plot_shaded
ggsave(file.path(OUT_DIR, "SupplementaryFigure_indel_umap_shaded.png"))

