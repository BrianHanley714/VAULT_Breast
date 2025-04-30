# Draw the phylogenies from the simulated tumour which was the best fit for breast cancer

rm(list = ls(all = TRUE))
# LIBRARIES ---------------------------------------------------------------
library(ape)
library(tidyverse)
# PATHS -------------------------------------------------------------------
run = "run001"

BASE = here::here()
BASE = "/Users/hanleyb/Dropbox (The Francis Crick)/HoLSTF_Breast/Github_Repo"
OUT_DIR = file.path(BASE, "analysis", "figures")


# LOAD DATA ---------------------------------------------------------------

input_dir = file.path(BASE, "data/simulations", run, "out/")
  if(file.exists(input_dir)){
    mat_out = read.delim(paste0(input_dir, "/matrix_ancestry.txt"), header = T)
    tree = read.tree(paste0(input_dir, "/clone_tree.new"))
    sample_matrix = read.delim(paste0(input_dir, "/sample_matrix.txt"), header = T)
    mat_test = as_tibble(read.delim(paste0(input_dir, "/sorted_matrix.txt"), header = T))
  }




# ORGANISE PLOT DATA ------------------------------------------------------
# set seed for reproducible sampling
set.seed(100)
names(mat_test) = c("V1", "V2", "V3")
rm(index_out)
for( i in sample(1:nrow(sample_matrix), 10)){
  
  indices = mat_test%>%
    mutate(index = 1:nrow(mat_test))%>%
    filter(V1>sample_matrix$Var1[i]&V1<sample_matrix$Var4[i])%>%
    filter(V2>sample_matrix$Var2[i]&V2<sample_matrix$Var5[i])%>%
    filter(V3>sample_matrix$Var3[i]&V3<sample_matrix$Var6[i])%>%
    pull(index)
  indices = mat_out$ID[indices]
  indices = paste(indices, collapse = ",")
  if(exists("index_out")==TRUE){index_out = rbind(index_out,indices)}
  if(exists("index_out")==FALSE){index_out = indices}
}  

colours <- c(
  "#FF0000",  # Red
  "#0000FF",  # Blue
  "#008000",  # Green
  "#800080",  # Purple
  "#FFA500",  # Orange
  "#A52A2A",  # Brown
  "#000000",  # Black
  "#FFD700",  # Gold
  "#DC143C",  # Crimson
  "#2E8B57"   # Sea Green
)
rm(mat_col)
for(i in 1:nrow(index_out)){
  vec = unique(c(str_split(index_out[i,], ",", simplify = T)))
  mat = cbind(vec, i)
  if(exists("mat_col")==TRUE){mat_col = as_tibble(rbind(mat_col,mat))}
  if(exists("mat_col")==FALSE){mat_col = mat}
}

mat_col = bind_cols(mat_col, colours[as.numeric(mat_col$i)])

ids = c(tree$tip.label, tree$node.label)
ids = ids[tree$edge[,2]]

edge_col = if_else(ids %in% mat_col$vec, mat_col$...3[match(ids, mat_col$vec)], "#F0F0F0")

tip_col = if_else(tree$tip.label %in% mat_col$vec, mat_col$...3[match(tree$tip.label, mat_col$vec)], "#F0F0F0")

no.repeats = as_tibble(cbind(as.numeric(table(mat_col$vec)), as.numeric(names(table(mat_col$vec)))))
tree$tip.label = if_else(tree$tip.label %in% no.repeats$V2, as.character(no.repeats$V1[match(tree$tip.label, no.repeats$V2)]), "")
tree$node.label = if_else(tree$node.label %in% no.repeats$V2, as.character(no.repeats$V1[match(tree$node.label, no.repeats$V2)]), "")
####### Get colours for the rep_sample
blocks_to_remove = sample_matrix[sample(c(1:nrow(sample_matrix)), size = nrow(sample_matrix)*0.1, replace = F),]
mat_left_over = mat_test%>%
  mutate(index = 1:nrow(mat_test))
for(x in 1:nrow(blocks_to_remove)){  
  mat_left_over = mat_left_over%>%
    filter(!(V1>blocks_to_remove$Var1[x]&V1<blocks_to_remove$Var4[x]&V2>blocks_to_remove$Var2[x]&V2<blocks_to_remove$Var5[x]&V3>blocks_to_remove$Var3[x]&V3<blocks_to_remove$Var6[x]))
  
}  
indices_rs = mat_left_over%>%
  sample_n(115000, replace = T)%>% # you can sample from the same population of cells twice
  pull(index)
tree_rs = tree
ids_rs = c(tree_rs$tip.label, tree_rs$node.label)
ids_rs = ids[tree_rs$edge[,2]]

rep_seq_col = "#54278F"

edge_col_rs = if_else(ids_rs %in% indices_rs, rep_seq_col, "#F0F0F0")

tip_col_rs = if_else(tree$tip.label %in% indices_rs, rep_seq_col, "#F0F0F0")




# WRITE PLOT --------------------------------------------------------------
png(file.path(OUT_DIR, "Figure2A_illustrative_phylogenies.png"), width = 1000, height = 600, res = 150)

layout(matrix(c(1, 2), nrow = 1), widths = c(3, 3))
par(mar = c(5, 4, 4, 2))
plot.phylo(tree, "fan", use.edge.length = T,
           show.tip.label = F,
           edge.color = edge_col_rs,
           tip.color = tip_col_rs, 
           edge.width = 1,show.node.label = F)

plot.phylo(tree, "fan", 
           use.edge.length = T,
           show.tip.label = F,
           edge.color = edge_col,
           tip.color = "light grey",
           node.color = "light grey",
           edge.width = 1,
           show.node.label = F)



dev.off()


# DRAW LEGEND -------------------------------------------------------------
plot.new()
legend("bottomleft",title = "Subclones captured in:",
       legend = c( "RepSamp",paste0("SingReg", c(1:10)),
                   "Subclone missed"),
       fill = c(rep_seq_col, colours, "lightgrey"),
       horiz = T,
       #border = "black",
       #cex = 1,
       text.width = 0.06,
       x.intersp = 0.1,
       bty = "n",
       #pch = 15,
       #lwd = 2,
       lty = FALSE,pt.lwd = FALSE,
       
       xpd = F
)

