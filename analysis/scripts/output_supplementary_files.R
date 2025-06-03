# Create a supplmentary excel document for submission with paper

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
#library(tidyverse)
library(tools)
library(openxlsx)
# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "data", "metadata")
DATA = file.path(BASE, "data")


paths = list.files(DATA, recursive = T, pattern = ".txt$|.tsv$|.xlsx$", full.names = T)
paths = grep("simulations", paths, invert = T, value = T)
paths = grep("VAULT_manuscript_supplementary.xlsx", paths, invert = T, value = T)
text_files = lapply(paths, function(fn){
    print(fn)
    if(file_ext(fn) == "xlsx"){
    return(read.xlsx(fn, fillMergedCells = T))
                   
  }
  else{
  return(read.delim(fn))
    }
  })

names(text_files) = basename(paths)
names(text_files)%>%sort

out_list = list(
  "ST1_Excluded_Cases" = text_files$Cases_excluded.xlsx,
  "ST2_VAULTPatients" = text_files$clinical_data.txt,
  "ST3_Recurrences" = text_files$Recurrence_data.xlsx,
  "ST4_ImageAnalysis" = text_files$tumour_cell_counts_w_calcs.tsv,
  "ST5_VAULTvariants" = text_files$variant_calls_VAULT.txt,
  "ST6_TCGA_MBvariants" = text_files$variant_calls_TCGA_MB.txt,
  "ST7_MatchedPatients" = text_files$matched_patients_characteristics.txt,
  "ST8_Ki67expression" = text_files$Ki67_pathologist_scores.txt,
  "ST9_FixedFlow" = text_files$fixed_flow.txt,
  "ST10_CustomPanel" = text_files$custom_panel_coordinates.txt,
  "ST11_Clinical_Trial" = text_files$trials_by_gene.txt
)



hs <- createStyle(
  textDecoration = "BOLD",wrapText = T, fontColour = "black", fontSize = 12,border = "TopBottomLeftRight", borderColour = "black", borderStyle = rep("thick", 4),
  fontName = "Arial Narrow", fgFill = "white"
)





# set the border colour
openxlsx::write.xlsx(out_list, file.path(OUT_DIR, "VAULT_manuscript_supplementary.xlsx"), headerStyle = hs, borderColour = "black", borders = "all")

