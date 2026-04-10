#devtools::install_github("im3sanger/dndscv")
library(dndscv)
library(tidyverse)
library(qs)
library(ggnewscale) 
RS_col = "#6A3D9Adb" 
SR_col = "#33A02Cb7"



# PATHS -------------------------------------------------------------------
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "figures")
VAULT_VARIANTS_PATH = file.path(BASE, "data", "variants", "variant_calls_VAULT.rdata")

maf_out = qread(VAULT_VARIANTS_PATH)

mut = maf_out%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault)%>%
  dplyr::select(sampleID = Tumor_Sample_Barcode,
                chr = Chromosome,
                pos = Start_Position,
                ref = Reference_Allele,
                mut = Tumor_Seq_Allele2
  )

res <- dndscv(mut)
head(res$sel_cv)
res$globaldnds
sel_cv = res$sel_cv
sign_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

maf_out%>%
  mutate(global_q_value = sign_genes$qglobal_cv[match(Hugo_Symbol, sign_genes$gene_name)],
         global_q_value = if_else(global_q_value == 0, 1e-09,global_q_value)
  )%>%
  dplyr::filter(Tumor_Sample_Barcode %in% ck_included_cases_vault,
                grepl("CLON", clonality),
                Hugo_Symbol %in% sign_genes$gene_name)%>%
  mutate(x_axis = paste0(Tumor_Sample_Barcode, Hugo_Symbol, Start_Position))%>%
  ggplot()+
  geom_point(aes(reorder(x_axis, 1-ccf_expected_copies), ccf_expected_copies, col = clonality))+
  geom_errorbar(aes(reorder(x_axis, 1-ccf_expected_copies), ccf_expected_copies, col = clonality, ymin = ccf_expected_copies_lower, ymax = ccf_expected_copies_upper),
                width = 0.2)+
  facet_wrap(~Hugo_Symbol, scales = "free_x", nrow = 1)+
  theme_classic(base_size = basesize)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom")+
  scale_color_manual(values = c("CLONAL" = RS_col, "SUBCLONAL" = SR_col))+
  new_scale_colour() +
  geom_hline((aes(reorder(x_axis, 1-ccf_expected_copies), yintercept = 1.1, size = 10,  col = global_q_value)), , shape = 22)+
  scale_color_viridis_c( trans  = "log10",
                         limits = c(1e-09, max(sign_genes$qglobal_cv)),
                         breaks = c(1e-9, 1e-6, 1e-3, 1e-1),
                         labels = scales::label_log()
  )+
  guides(size = "none")+
  ylab("Clonality")


ggsave(file.path(OUT_DIR, "3E_dndscv.png"))

