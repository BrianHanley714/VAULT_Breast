# this is a custom CONIPHER plotting function forked from CONIPHER package (https://github.com/McGranahanLab/CONIPHER)
# edits were made by Brian Hanley

# Plot trees for CK -------------------------------------------------------
plot_pyclone_tree = function(sample_pyclone_tree, node_size_reduction, ccf_y, sampleID, driver_line_width = 0.04){
  pyclone_tree <- sample_pyclone_tree$graph_pyclone
  g <- graph.data.frame(pyclone_tree$default_tree, directed = FALSE)
  indx <- V(g)$name
  vcol <- setNames(color.tree(pyclone_tree$edgelength), 
                   names(pyclone_tree$edgelength))[indx]
  l <- layout_as_tree(g, root = pyclone_tree$trunk)
  pie.size <- ncol(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)
  node.shape <- setNames(rep("pie", length(vcol)), names(vcol))
  pie.slices <- lapply(1:length(vcol), function(x) rep(1, 
                                                       pie.size))
  empty.col = "light grey"
  node_size_factor <- log2(max(pyclone_tree$edgelength))/30
  node.size <- log2(pyclone_tree$edgelength)/node_size_factor
  node.size <- node.size[names(node.shape)]
  pie.colors <- sample_pyclone_tree$nested_pyclone$ccf_cluster_table[match(names(vcol), 
                                                                           rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)), 
                                                                     , drop = F]
  
  pie.colors <- ifelse(pie.colors >= 90, 99, pie.colors)
  pie.colors <- ifelse(pie.colors < 10 & pie.colors >= 
                         1, 10, pie.colors)
  
  ### BH edits
  names = colnames(pie.colors)
  unenriched_Rs = which(grepl("RT001", names))
  ck_enriched = which(grepl("ET001", names))
  ki67_enriched = which(grepl("ET002", names))
  ffpe_info = read.xlsx('/Users/hanleyb/The Francis Crick Dropbox/Brian Hanley/HoLSTF_Breast/FFPE_sample_info.xlsx')
  
  met_sr = which(names %in% ffpe_info$Sample_ID[ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  prim_sr = which(names %in% ffpe_info$Sample_ID[!ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  
  cols = vector(length = length(names))
  cols[unenriched_Rs] = "#E31A1Cdb"
  cols[ck_enriched] = "#1F78B4db"
  cols[ki67_enriched] = "#6A3D9Aff"
  cols[prim_sr] = "#33A02Cff"
  cols[met_sr] = "#FF7F00db"
  
  colour_frames = rep("grey", length(prim_sr))
  
  
  if(length(met_sr >=1)){
    if(length(met_sr) == 1){
      colour_frames = if_else(sample_pyclone_tree$nested_pyclone$ccf_cluster_table[,met_sr]>0, "black", "grey")  
    }
    
    if(length(met_sr) > 1){
      colour_frames = if_else(rowSums(sample_pyclone_tree$nested_pyclone$ccf_cluster_table[,met_sr])>0, "black", "grey")
    }
    
    colour_frames = colour_frames[match(names(vcol),
                                        rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table))
                                  , drop = F]
  }
  
  
  pie.colors <- lapply(1:nrow(pie.colors), function(x) {
    if (!all(is.na(pie.colors[x, ]))) {
      tmp <- pie.colors[x, ]
      tmp2 <- tmp
      
      tmp = sapply(1:length(tmp), function(fn){
        colfunc = colorRampPalette(c("white", cols[fn]))
        speccolours <- colfunc(100)
        if(tmp[fn]>0){ return(speccolours[tmp[fn]])}
        if(tmp[fn]==0){ return(empty.col)}
        
      })
      
      
    }
  })
  
  pie.colors
  
  
  g_dir <- graph.data.frame(pyclone_tree$default_tree, 
                            directed = TRUE)
  edges <- get.edgelist(g_dir)
  ecol <- setNames(rep("#bdbdbd", nrow(edges)), edges[, 
                                                      2])
  ewidth <- rep(1, length(ecol))
  ecol[paste(edges[, 1], edges[, 2], sep = ":") %in% pyclone_tree$consensus_relationships] <- "#000000"
  ewidth[paste(edges[, 1], edges[, 2], sep = ":") %in% 
           pyclone_tree$consensus_relationships] <- 150
  
  vcol[indx]
  
  
  
  clusters = unique(sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster)
  drivers = drivers%>%mutate(id = paste(Trial_ID, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":"))
  
  drivers_tumour = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    
    if(length(unique(drivers%>%
                     dplyr::filter(id %in% variants_in_cluster)%>%pull(Hugo_Symbol)))>100){drivers_in_cluster = ">10drivers"}
    else{drivers_in_cluster = paste(unique(drivers%>%
                                             dplyr::filter(id %in% variants_in_cluster)%>%
                                             mutate(out = paste(Hugo_Symbol, HGVSp_short, sep = "_"))%>%
                                             pull(out)), collapse = "\n")}
    
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  called_in_rs = lapply(clusters, function(cluster){
    factor_order = as.vector(str_split(drivers_tumour[[which(clusters == cluster)]], "\n", simplify = T))
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = drivers%>%
      dplyr::filter(id %in% variants_in_cluster)%>%
      mutate(factor = factor(paste(Hugo_Symbol, HGVSp_short, sep = "_"), levels = factor_order))%>%
      dplyr::filter(!duplicated(factor))%>%
      arrange(factor)%>%pull(called_in_rs)%>%paste(., collapse = "\n")
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  
  
  called_in_sr = lapply(clusters, function(cluster){
    factor_order = as.vector(str_split(drivers_tumour[[which(clusters == cluster)]], "\n", simplify = T))
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = drivers%>%
      dplyr::filter(id %in% variants_in_cluster)%>%
      mutate(factor = factor(paste(Hugo_Symbol, HGVSp_short, sep = "_"), levels = factor_order))%>%
      dplyr::filter(!duplicated(factor))%>%
      arrange(factor)%>%pull(called_in_ffpe)%>%paste(., collapse = "\n")
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  driver_colours = if_else(unlist(called_in_rs)=="TRUE"& unlist(called_in_sr) == "TRUE", "black", 
          if_else(unlist(called_in_rs)=="TRUE"& unlist(called_in_sr) == "FALSE", "blue",
                  if_else(unlist(called_in_rs)=="FALSE"& unlist(called_in_sr) == "TRUE", "red", NA
  )))
  
  driver_colours_rs = unlist(called_in_rs)[match(V(g)$name, as.character(clusters))]
  driver_colours_sr = unlist(called_in_sr)[match(V(g)$name, as.character(clusters))]
  driver_per_clust = do.call("c", drivers_tumour)[match(V(g)$name, as.character(clusters))]
  #driver_type_per_cluster = do.call("c", driver_type_cluster)[match(V(g)$name, as.character(clusters))]
  #driver_type_per_cluster = as.numeric(driver_type_per_cluster)
  
  
  # change orientation of trees   
  l[,2] = -l[,2]
  l[,1] = -l[,1]
  
  
   return(
  {
    par(mar = c(0, 0, 0, 0))
    coords <- norm_coords(l, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
    plot(g, layout = coords, main = "", vertex.color = vcol[indx], 
              vertex.frame.color = colour_frames, vertex.shape = node.shape, 
         vertex.label = NA, 
         #vertex.label.color= driver_colours,
              vertex.lwd = 5, vertex.pie.lwd = 3, vertex.pie = pie.slices, 
              vertex.pie.color = lapply(pie.colors, rev), vertex.size = node.size/node_size_reduction, 
              edge.color = ecol, edge.size = ewidth, vertex.label.cex = 0.75, 
              vertex.label.degree = pi,
         margin = 0,
              vertex.label.pos = 2, vertex.label.dist = 5, vertex.label.family = "Helvetica", 
              vertex.label.font = 1, vertex.label.color = "black")
    
    #coords <- l          # layout matrix (x, y)
    labels <- driver_per_clust
    
    for (i in seq_along(labels)) {
      
      lines <- strsplit(labels[i], "\n")[[1]]
      if(length(lines)>0){
      lines = str_split(lines, "_", simplify = T)[,1]
      }
      rs_calls = strsplit(driver_colours_rs[i], "\n")[[1]]
      sr_calls = strsplit(driver_colours_sr[i], "\n")[[1]]
      x_offset = (node.size/node_size_reduction)[i]
      # choose colours per line (example)
      
      #line_cols <- c("red", "blue", "darkgreen")[seq_along(lines)]
      line_cols = if_else(rs_calls == "TRUE" &sr_calls == "TRUE", "black",
                          if_else(rs_calls == "FALSE" &sr_calls == "TRUE", "#4575B4",
                                  if_else(rs_calls == "TRUE" &sr_calls == "FALSE", "#D73027",
                                          if_else(rs_calls == "FALSE" &sr_calls == "FALSE", "darkgreen", NA
                          ))))
      
      # vertical offsets so lines don't overlap
      y_offsets <- seq(0, by = -driver_line_width, length.out = length(lines))
      title(sampleID, line = -2)
      for (j in seq_along(lines)) {
       text(
          x = coords[i, 1]-x_offset/250,
          y = coords[i, 2] + y_offsets[j]+(0.02*length(lines) -0.02),
          labels = lines[j],
          col = line_cols[j],
          cex = 0.5,
          font = 1,
          family = "Helvetica",
          pos = 2
        )
      }
    }
}
   )
 }
sampleID = "HF039"
node_size_reduction = 1


# Plotting trees for Ki67 -------------------------------------------------

plot_pyclone_tree_Ki67 = function(sample_pyclone_tree, node_size_reduction, ccf_y, sampleID, driver_line_width = 0.04, fontsize =1){
  pyclone_tree <- sample_pyclone_tree$graph_pyclone
  g <- graph.data.frame(pyclone_tree$default_tree, directed = FALSE)
  indx <- V(g)$name
  vcol <- setNames(color.tree(pyclone_tree$edgelength), 
                   names(pyclone_tree$edgelength))[indx]
  l <- layout_as_tree(g, root = pyclone_tree$trunk)
  pie.size <- ncol(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)
  node.shape <- setNames(rep("pie", length(vcol)), names(vcol))
  pie.slices <- lapply(1:length(vcol), function(x) rep(1, 
                                                       pie.size))
  empty.col = "light grey"
  node_size_factor <- log2(max(pyclone_tree$edgelength))/30
  node.size <- log2(pyclone_tree$edgelength)/node_size_factor
  node.size <- node.size[names(node.shape)]
  pie.colors <- sample_pyclone_tree$nested_pyclone$ccf_cluster_table[match(names(vcol), 
                                                                           rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table)), 
                                                                     , drop = F]
  
  pie.colors <- ifelse(pie.colors >= 90, 99, pie.colors)
  pie.colors <- ifelse(pie.colors < 10 & pie.colors >= 
                         1, 10, pie.colors)
  
  
  names = colnames(pie.colors)
  unenriched_Rs = which(grepl("RT001", names))
  ck_enriched = which(grepl("ET001", names))
  ki67_enriched = which(grepl("ET002", names))
  #ffpe_info = read.delim(file.path(BASE,"data", "metadata","FFPE_sample_info.txt"))
  ffpe_info = read.xlsx('/Users/hanleyb/The Francis Crick Dropbox/Brian Hanley/HoLSTF_Breast/FFPE_sample_info.xlsx')
  met_sr = which(names %in% ffpe_info$Sample_ID[ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  prim_sr = which(names %in% ffpe_info$Sample_ID[!ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  
  cols = vector(length = length(names))
  cols[unenriched_Rs] = "#E31A1Cdb"
  cols[ck_enriched] = "#1F78B4db"
  cols[ki67_enriched] = "#6A3D9Aff"
  cols[prim_sr] = "#33A02Cff"
  cols[met_sr] = "#FF7F00db"
  
  colour_frames = rep("grey", length(prim_sr))
  
  
  if(length(met_sr >=1)){
    if(length(met_sr) == 1){
      colour_frames = if_else(sample_pyclone_tree$nested_pyclone$ccf_cluster_table[,met_sr]>0, "black", "grey")  
    }
    
    if(length(met_sr) > 1){
      colour_frames = if_else(rowSums(sample_pyclone_tree$nested_pyclone$ccf_cluster_table[,met_sr])>0, "black", "grey")
    }
    
    colour_frames = colour_frames[match(names(vcol),
                                        rownames(sample_pyclone_tree$nested_pyclone$ccf_cluster_table))
                                  , drop = F]
  }
  
  pie.colors <- lapply(1:nrow(pie.colors), function(x) {
    if (!all(is.na(pie.colors[x, ]))) {
      tmp <- pie.colors[x, ]
      tmp2 <- tmp
      
      tmp = sapply(1:length(tmp), function(fn){
        colfunc = colorRampPalette(c("white", cols[fn]))
        speccolours <- colfunc(100)
        if(tmp[fn]>0){ return(speccolours[tmp[fn]])}
        if(tmp[fn]==0){ return(empty.col)}
        
      })
      
    }
  })
  
  
  
  
  g_dir <- graph.data.frame(pyclone_tree$default_tree, 
                            directed = TRUE)
  edges <- get.edgelist(g_dir)
  ecol <- setNames(rep("#bdbdbd", nrow(edges)), edges[, 
                                                      2])
  ewidth <- rep(1, length(ecol))
  ecol[paste(edges[, 1], edges[, 2], sep = ":") %in% pyclone_tree$consensus_relationships] <- "#000000"
  ewidth[paste(edges[, 1], edges[, 2], sep = ":") %in% 
           pyclone_tree$consensus_relationships] <- 150
  
  vcol[indx]
  
  
  
  clusters = unique(sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster)
  drivers = drivers%>%mutate(id = paste(Trial_ID, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":"))
  
  drivers_tumour = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    
    if(length(unique(drivers%>%
                     dplyr::filter(id %in% variants_in_cluster)%>%pull(Hugo_Symbol)))>100){drivers_in_cluster = ">10drivers"}
    else{drivers_in_cluster = paste(unique(drivers%>%
                                             dplyr::filter(id %in% variants_in_cluster)%>%
                                             mutate(out = paste(Hugo_Symbol, HGVSp_short, sep = "_"))%>%
                                             pull(out)), collapse = "\n")}
    
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  # driver_type_cluster = lapply(clusters, function(cluster){
  #   variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
  #   drivers_in_cluster = unique(drivers%>%
  #                                 dplyr::filter(id %in% variants_in_cluster)%>%pull(tiered_clinical_importance))
  #   if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
  #   return(drivers_in_cluster)})
  # 
  
  called_in_rs = lapply(clusters, function(cluster){
    factor_order = as.vector(str_split(drivers_tumour[[which(clusters == cluster)]], "\n", simplify = T))
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = drivers%>%
      dplyr::filter(id %in% variants_in_cluster)%>%
      mutate(factor = factor(paste(Hugo_Symbol, HGVSp_short, sep = "_"), levels = factor_order))%>%
      dplyr::filter(!duplicated(factor))%>%
      arrange(factor)%>%pull(called_in_rs)%>%paste(., collapse = "\n")
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  
  called_in_rs_ki67 = lapply(clusters, function(cluster){
    factor_order = as.vector(str_split(drivers_tumour[[which(clusters == cluster)]], "\n", simplify = T))
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = drivers%>%
      dplyr::filter(id %in% variants_in_cluster)%>%
      mutate(factor = factor(paste(Hugo_Symbol, HGVSp_short, sep = "_"), levels = factor_order))%>%
      dplyr::filter(!duplicated(factor))%>%
      arrange(factor)%>%pull(called_in_rs_ki67)%>%paste(., collapse = "\n")
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  called_in_sr = lapply(clusters, function(cluster){
    factor_order = as.vector(str_split(drivers_tumour[[which(clusters == cluster)]], "\n", simplify = T))
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = drivers%>%
      dplyr::filter(id %in% variants_in_cluster)%>%
      mutate(factor = factor(paste(Hugo_Symbol, HGVSp_short, sep = "_"), levels = factor_order))%>%
      dplyr::filter(!duplicated(factor))%>%
      arrange(factor)%>%pull(called_in_ffpe)%>%paste(., collapse = "\n")
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
   driver_per_clust = do.call("c", drivers_tumour)[match(V(g)$name, as.character(clusters))]
  # driver_type_per_cluster = do.call("c", driver_type_cluster)[match(V(g)$name, as.character(clusters))]
  # driver_type_per_cluster = as.numeric(driver_type_per_cluster)
  driver_colours_rs = unlist(called_in_rs)[match(V(g)$name, as.character(clusters))]
  driver_colours_sr = unlist(called_in_sr)[match(V(g)$name, as.character(clusters))]
  driver_colours_ki67 = unlist(called_in_rs_ki67)[match(V(g)$name, as.character(clusters))]
  
  
  l[,2] = -l[,2]
  if(#sampleID %in% c("HF041", "HF057", "HF039", "HF299")
    sampleID %in% c("HF283", "HF388", "HF420", "HF105")){
    l[,1] = -l[,1]  
  }
  
  return(
    {
    par(mar = c(0, 0, 0, 0))
    coords <- norm_coords(l, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
    plot(g, layout = coords, main = "", vertex.color = vcol[indx], 
         vertex.frame.color = colour_frames, vertex.shape = node.shape, 
         vertex.label = NA, 
         #vertex.label.color= driver_colours,
         vertex.lwd = 5, vertex.pie.lwd = 3, vertex.pie = pie.slices, 
         vertex.pie.color = lapply(pie.colors, rev), vertex.size = node.size/node_size_reduction, 
         edge.color = ecol, edge.size = ewidth, vertex.label.cex = 0.75, 
         vertex.label.degree = pi,
         margin = 0,
         vertex.label.pos = 2, vertex.label.dist = 5, vertex.label.family = "Helvetica", 
         vertex.label.font = 1, vertex.label.color = "black")
    
    #coords <- l          # layout matrix (x, y)
    labels <- driver_per_clust
    
    for (i in seq_along(labels)) {
      
      lines <- strsplit(labels[i], "\n")[[1]]
      if(length(lines)>0){
        lines = str_split(lines, "_", simplify = T)[,1]
      }
      rs_calls = strsplit(driver_colours_rs[i], "\n")[[1]]
      sr_calls = strsplit(driver_colours_sr[i], "\n")[[1]]
      ki67_calls = strsplit(driver_colours_ki67[i], "\n")[[1]]
      
      x_offset = (node.size/node_size_reduction)[i]
      # choose colours per line (example)
      
      #line_cols <- c("red", "blue", "darkgreen")[seq_along(lines)]
      line_cols = if_else(rs_calls == "TRUE" &sr_calls == "TRUE" & ki67_calls=="TRUE", "black",
                          if_else(rs_calls == "FALSE" &sr_calls == "TRUE", "#4575B4",
                                  if_else(ki67_calls=="TRUE" & rs_calls == "FALSE" & sr_calls == "FALSE", "#6A3D9Aff",
                                  if_else(rs_calls == "TRUE" &sr_calls == "FALSE", "#D73027",
                                          if_else(rs_calls == "FALSE" &sr_calls == "FALSE" &ki67_calls == "FALSE", "darkgreen", NA
                                          )))))
      
      # vertical offsets so lines don't overlap
      y_offsets <- seq(0, by = -driver_line_width, length.out = length(lines))
      title(sampleID, line = -2)
      for (j in seq_along(lines)) {
        text(
          x = coords[i, 1]-x_offset/250,
          y = coords[i, 2] + y_offsets[j]+(0.02*length(lines) -0.02),
          labels = lines[j],
          col = line_cols[j],
          cex = fontsize,
          font = 1,
          family = "Helvetica",
          pos = 2
        )
      }
    }
    }
)
}
  


