#### Andrea Del Cortona
#### 2023/05/25



#==============================================================================#
# 0 - Environment Setup                                                     ####

## install packages
#devtools::install_github("gaospecial/ggVennDiagram")
#BiocManager::install("PCAtools")
#install.packages("openxlsx")
#"tidytree", 

# import libraries
up_packages = c("ape", "circlize", "ComplexHeatmap", "corrplot", "DECIPHER", "dendextend", "doSNOW", "foreach",
                "ggdist", "ggplot2", "ggtree", "gridExtra", "phangorn", "phylogram", "phytools", "reshape",
                "RColorBrewer", "scales", "stringr", "TreeDist", "treeio")
lapply(up_packages, require, character.only = TRUE)

## Set working directory
mainDir = "/home/Andrea/Desktop/CP_MT/ANALYSES_AND_RESULTS/"
setwd(mainDir)



#==============================================================================#
# 1 - Functions                                                             ####


#==============================================================================#
# 2 - Tree comparisons                                                      ####

#------------------------------------------------------------------------------#
## 2.1 - import and sort trees                                              ####

# declare gene lists
cp_genes = c("petB", "psaA", "psaB", "psbB", "psbD", "rps2")
mt_genes = c("atp6", "cox1", "cox2", "rps3")

# import ref species tree
cp_mt_ML_dendro = ape::read.tree(file = "06_cp_mt_concat_ML/cp_mt_allgenes_concat.contree")

# create empty list of trees
combinatorial_trees = list("cp_mt_conc_2_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_3_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_4_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_5_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_6_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_7_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_8_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_9_genes" = list("cp_mt_ML" = cp_mt_ML_dendro),
                           "cp_mt_conc_10_genes" = list("cp_mt_ML" = cp_mt_ML_dendro))

# populate combinatorial trees
for(k in 2:10){
  
  # get file list
  tree_list_all = list.files(path = paste(mainDir, "/09_combinatorial_phylogeny/01_cp_mt_combinatorial/cp_mt_conc_", k, "_genes", sep = ""),
                             pattern = "\\.contree$")
  # keep only trees with both cp and mt markers
  tree_list_fltr = c()
  for(tree in tree_list_all){
    cp = FALSE
    mt = FALSE
    # check if cp is present
    for(gene in cp_genes){
      if(grepl(gene, tree, fixed = TRUE)){
        cp = TRUE
      }
    }
    # check if mt is present
    for(gene in mt_genes){
      if(grepl(gene, tree, fixed = TRUE)){
        mt = TRUE
      }
    }
    if(cp == TRUE & mt == TRUE){
      tree_list_fltr = c(tree_list_fltr, tree)
    }
  }
  
  # import trees
  for(tree in tree_list_fltr){
    combinatorial_trees[[paste("cp_mt_conc_", k, "_genes", sep = "")]][[length(combinatorial_trees[[paste("cp_mt_conc_", k, "_genes", sep = "")]]) + 1]] = ape::read.tree(file = paste(mainDir, "/09_combinatorial_phylogeny/01_cp_mt_combinatorial/cp_mt_conc_", k, "_genes/", tree, sep = ""))
    names(combinatorial_trees[[paste("cp_mt_conc_", k, "_genes", sep = "")]])[[length(combinatorial_trees[[paste("cp_mt_conc_", k, "_genes", sep = "")]])]] = stringr::str_remove(tree, ".aln.contree")
  }
  
  # clean
  rm(cp, gene, k, mt, tree, tree_list_all, tree_list_fltr)
}



#------------------------------------------------------------------------------#
## 2.2 - Get distance matrices                                              ####

# create distance lists
combinatorial_distances = list("cp_mt_conc_2_genes" = NULL,
                               "cp_mt_conc_3_genes" = NULL,
                               "cp_mt_conc_4_genes" = NULL,
                               "cp_mt_conc_5_genes" = NULL,
                               "cp_mt_conc_6_genes" = NULL,
                               "cp_mt_conc_7_genes" = NULL,
                               "cp_mt_conc_8_genes" = NULL,
                               "cp_mt_conc_9_genes" = NULL,
                               "cp_mt_conc_10_genes" = NULL)

# populate distance lists
combinatorial_list = c(25, 97, 195, 247, 210, 121, 46, 11, 2)
for(k in 1:length(combinatorial_distances)){
  
  # get size of the results matrix
  empty_matrix = matrix(nrow = combinatorial_list[[k]],
                        ncol = combinatorial_list[[k]])
  colnames(empty_matrix) = names(combinatorial_trees[[k]])
  rownames(empty_matrix) = names(combinatorial_trees[[k]])
  
  # populate
  combinatorial_distances[[k]] = list("common_nodes_corr" = empty_matrix,
                                      "RF" = empty_matrix,
                                      "KF" = empty_matrix,
                                      "MAST_%" = empty_matrix)
  
  # clean
  rm(empty_matrix)
}

# get pairwise distances
for(j in 1:length(combinatorial_trees)){
  for(i in 1:length(combinatorial_trees[[j]])){
    for(k in 1:length(combinatorial_trees[[j]])){
      
      # get tree
      tree1_tmp = combinatorial_trees[[j]][[i]]
      tree2_tmp = combinatorial_trees[[j]][[k]]

      # list of common species
      species_list = tree1_tmp$tip.label[which(tree1_tmp$tip.label %in% tree2_tmp$tip.label)]

      # get overlapping species
      tree1 = ape::keep.tip(tree1_tmp, species_list)
      tree2 = ape::keep.tip(tree2_tmp, species_list)
      
      # get distances
      distances = phangorn::treedist(tree1, tree2)
      combinatorial_distances[[j]][["RF"]][i, k] = distances[[1]]
      combinatorial_distances[[j]][["KF"]][i, k] = distances[[2]]
      combinatorial_distances[[j]][["MAST_%"]][i, k] = length(phangorn::mast(tree1, tree2, tree = FALSE)) / length(species_list)
      
      # clean
      rm(tree1, tree2, tree1_tmp, tree2_tmp, distances)
    }
  }
}


### get pairwise common_nodes_corr distances

# iterate all combinatorial trees
for(j in 1:length(combinatorial_trees)){
  
  # create CPU cluster
  cl = parallel::makeCluster(8, type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  
  # get chronos in parallel
  chronos_list = foreach(k = 1:length(combinatorial_trees[[j]]), .combine = "c") %dopar% {
    library("ape")
    # get tree
    tryCatch(chronos(combinatorial_trees[[j]][[k]]), error = function(e) { NA })
  }  
  parallel::stopCluster(cl)
  names(chronos_list) = names(combinatorial_trees[[j]])
  
  # get pairwise common_nodes_corr distances
  for(i in 1:length(combinatorial_trees[[j]])){
    for(k in 1:length(combinatorial_trees[[j]])){
      
      print(paste(j, i, k, sep = ","))
      
      # get trees
      tree1_tmp = chronos_list[[names(combinatorial_trees[[j]])[i]]]
      tree2_tmp = chronos_list[[names(combinatorial_trees[[j]])[k]]]
      
      # check trees and get distance
      if(all(!is.na(tree1_tmp)) & all(!is.na(tree2_tmp))){
        combinatorial_distances[[j]][["common_nodes_corr"]][i, k] = tryCatch(
          cor.dendlist(dendlist(as.dendrogram(ape::root(tree1_tmp, outgroup = c("Oviri", "Pakin"), resolve.root = TRUE)),
                                as.dendrogram(ape::root(tree2_tmp, outgroup = c("Oviri", "Pakin"), resolve.root = TRUE))),
                       method = "common_nodes")[2],
          error = function(e) { NA })
      } else {
        combinatorial_distances[[j]][["common_nodes_corr"]][i, k] = NA
      }
    }
  }
}


## export the distance matrices
for(i in 1:length(combinatorial_distances)){
  write.table(combinatorial_distances[[i]][["common_nodes_corr"]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".05_combinatorial_stats.cp_mt_conc_", i + 1, "_genes.common_nodes_corr.txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  write.table(combinatorial_distances[[i]][["RF"]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".05_combinatorial_stats.cp_mt_conc_", i + 1, "_genes.RF.txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  write.table(combinatorial_distances[[i]][["KF"]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".05_combinatorial_stats.cp_mt_conc_", i + 1, "_genes.KF.txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  write.table(combinatorial_distances[[i]][["MAST_%"]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".05_combinatorial_stats.cp_mt_conc_", i + 1, "_genes.MAST_%.txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
}



#------------------------------------------------------------------------------#
## 2.3 - Plot distance matrices                                             ####


# 
# for(k in 1:length(combinatorial_distances)){
#   
#   # remove possible rows/columns with only NAs
#   tmp_matrix = distance_sets[[i]][[k]]
#   tmp_matrix = tmp_matrix[, colSums(is.na(tmp_matrix)) < nrow(tmp_matrix)]
#   tmp_matrix = tmp_matrix[rowSums(is.na(tmp_matrix)) < ncol(tmp_matrix), ]
#   
#   # get names
#   plot_name = ifelse(names(distance_sets[[i]])[[k]] == "common_nodes_corr", "common\nnodes\ncorr", names(distance_sets[[i]])[[k]])
#   
#   # prep annotation
#   organell = as.data.frame(ifelse(rownames(tmp_matrix) %in% cp_genes_list, "cp",
#                                   ifelse(rownames(tmp_matrix) %in% mt_genes_list, "mt", "cp_mt")))
#   rownames(organell) = rownames(tmp_matrix)
#   colnames(organell) = "organell"
#   organell$color = ifelse(organell$organell == "cp", "#469d89", ifelse(organell$organell == "mt", "steelblue", "red"))
#   
#   # plot raw distance
#   distances_plots[[i]][[k]] = grid.grabExpr(draw(ComplexHeatmap::Heatmap(tmp_matrix,
#                                                                          col = colorRamp2(c(min(tmp_matrix[!is.na(tmp_matrix)]),
#                                                                                             max(tmp_matrix[!is.na(tmp_matrix)])),
#                                                                                           c("white", "#D55E00")),
#                                                                          row_names_gp = gpar(col = organell$color),
#                                                                          column_names_gp = gpar(col = organell$color),
#                                                                          top_annotation = HeatmapAnnotation(organell = as.matrix(organell$organell),
#                                                                                                             show_annotation_name = FALSE,
#                                                                                                             show_legend = FALSE,
#                                                                                                             col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
#                                                                          left_annotation = rowAnnotation(organell = as.matrix(organell$organell),
#                                                                                                          show_annotation_name = FALSE,
#                                                                                                          col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
#                                                                          name = plot_name)))
#   grid.newpage(recording = TRUE)
#   grid.draw(distances_plots[[i]][[k]])
# 
# 
# }


#==============================================================================#
# 99 - Save R image file                                                    ####

save.image(file = "02_markers_selection.RData")
#load("02_markers_selection.RData")




