# this is a custom CONIPHER plotting function forked from CONIPHER package (https://github.com/McGranahanLab/CONIPHER)
# edits were made by Brian Hanley


# Plot trees for CK -------------------------------------------------------
plot_pyclone_tree = function(sample_pyclone_tree, node_size_reduction, ccf_y, sampleID){
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
  ffpe_info = read.delim(file.path(BASE,"data", "metadata","FFPE_sample_info.txt"))
  
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
  
  drivers_tumour = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    
    if(length(unique(drivers%>%
                     filter(id %in% variants_in_cluster)%>%pull(Hugo_Symbol)))>10){drivers_in_cluster = ">10drivers"}
    else{drivers_in_cluster = paste(unique(drivers%>%
                                             filter(id %in% variants_in_cluster)%>%
                                             mutate(out = paste(tiered_clinical_importance, Hugo_Symbol, sep = "-"))%>%
                                             pull(out)), collapse = "\n")}
    
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  driver_type_cluster = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = unique(drivers%>%
                                  filter(id %in% variants_in_cluster)%>%pull(tiered_clinical_importance))
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  driver_per_clust = do.call("c", drivers_tumour)[match(V(g)$name, as.character(clusters))]
  driver_type_per_cluster = do.call("c", driver_type_cluster)[match(V(g)$name, as.character(clusters))]
  driver_type_per_cluster = as.numeric(driver_type_per_cluster)
  

# change orientation of trees   
  l[,2] = -l[,2]
  l[,1] = -l[,1]

  
  return(plot(g, layout = l, main = sampleID, vertex.color = vcol[indx], 
              vertex.frame.color = colour_frames, vertex.shape = node.shape, vertex.label = driver_per_clust,
              vertex.lwd = 5, vertex.pie.lwd = 3, vertex.pie = pie.slices, 
              vertex.pie.color = lapply(pie.colors, rev), vertex.size = node.size/node_size_reduction, 
              edge.color = ecol, edge.size = ewidth, vertex.label.cex = 0.75, 
              vertex.label.degree = pi,
              vertex.label.pos = 2, vertex.label.dist = 5, vertex.label.family = "Helvetica", 
              vertex.label.font = 1, vertex.label.color = "black")
  )
}




# Plotting trees for Ki67 -------------------------------------------------

plot_pyclone_tree_Ki67 = function(sample_pyclone_tree, node_size_reduction, ccf_y, sampleID){
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
  ffpe_info = read.delim(file.path(BASE,"data", "metadata","FFPE_sample_info.txt"))
  
  met_sr = which(names %in% ffpe_info$Sample_ID[ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  prim_sr = which(names %in% ffpe_info$Sample_ID[!ffpe_info$Sample_Type %in% c("MET", "RECURRENCE_NRS")])
  
  cols = vector(length = length(names))
  cols[unenriched_Rs] = "#E31A1Cdb"
  cols[ck_enriched] = "#1F78B4db"
  cols[ki67_enriched] = "#6A3D9Aff"
  cols[prim_sr] = "#33A02Cff"
  cols[met_sr] = "#FF7F00db"
  
  
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
  
  drivers_tumour = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    
    if(length(unique(drivers%>%
                     filter(id %in% variants_in_cluster)%>%pull(Hugo_Symbol)))>10){drivers_in_cluster = ">10drivers"}
    else{drivers_in_cluster = paste(unique(drivers%>%
                                             filter(id %in% variants_in_cluster)%>%
                                             mutate(out = paste(tiered_clinical_importance, Hugo_Symbol, sep = "-"))%>%
                                             pull(out)), collapse = "\n")}
    
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  driver_type_cluster = lapply(clusters, function(cluster){
    variants_in_cluster = rownames(sample_pyclone_tree$ccf_table_absolute_clean)[sample_pyclone_tree$ccf_table_absolute_clean$PycloneCluster == cluster]
    drivers_in_cluster = unique(drivers%>%
                                  filter(id %in% variants_in_cluster)%>%pull(tiered_clinical_importance))
    if(length(drivers_in_cluster) == 0){drivers_in_cluster = ""}
    return(drivers_in_cluster)})
  
  driver_per_clust = do.call("c", drivers_tumour)[match(V(g)$name, as.character(clusters))]
  driver_type_per_cluster = do.call("c", driver_type_cluster)[match(V(g)$name, as.character(clusters))]
  driver_type_per_cluster = as.numeric(driver_type_per_cluster)
  

  
  l[,2] = -l[,2]
  if(#sampleID %in% c("HF041", "HF057", "HF039", "HF299")
    sampleID %in% c("HF283", "HF388", "HF420", "HF105")){
    l[,1] = -l[,1]  
  }
  
  return(plot(g, layout = l, main = sampleID, vertex.color = vcol[indx], 
              vertex.frame.color = "black", vertex.shape = node.shape, vertex.label = driver_per_clust,
              vertex.lwd = 5, vertex.pie.lwd = 3, vertex.pie = pie.slices, 
              vertex.pie.color = lapply(pie.colors, rev), vertex.size = node.size/node_size_reduction, 
              edge.color = ecol, edge.size = ewidth, vertex.label.cex = 0.75, 
              vertex.label.degree = pi,
              vertex.label.pos = 2, vertex.label.dist = 5, vertex.label.family = "Helvetica", 
              vertex.label.font = 1, vertex.label.color = "black")
  )
}

