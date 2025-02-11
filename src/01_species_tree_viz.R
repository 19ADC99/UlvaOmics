#### Andrea Del Cortona
#### 2022/12/20



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

### Create bootstrap vectors for colored nodes ---------------------------------
#' 
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note create a color vector for bootstrat nodes
#' @note 2022-12-22
#'
#' This function takes a dendrogram with node labels and, a string as ML or coalescent
#' and a vector with three colors.
#' It returns a vector of colors with:
#'   BS < 0.75          --> first color
#'   0.75 < BS < 0.95   --> second color
#'   BS > 0.95          --> third color
#' 
#' @param dendro        : a dendro to be recolored
#' @param type          : a string, "ML" or "coalescent"
#' @param colors        : a vector with three colors
color_nodes = function(dendro = dendro,
                       type = type,
                       colors = c("white", "grey75", "black")){
  
  # get nodes
  node_names = dendro %>% get_nodes_attr("label")
  # initialize output colors
  node_colors = c()
  
  # create breaks
  if(type == "ML"){
    my_breaks = c(75, 95)
  } else if(type == "coalescent"){
    my_breaks = c(0.75, 0.95)
  }    
  
  # populate node_colors vector
  for(i in 1:length(node_names)){
    if(node_names[[i]] == "Oviri"){
      node_colors = append(node_colors, colors[[3]])
    } else if(node_names[[i]] == "Pakin"){
      node_colors = append(node_colors, colors[[3]])
    } else if(startsWith(node_names[[i]], "U")){
      node_colors = append(node_colors, colors[[3]])
    } else {
      if(as.numeric(node_names[[i]]) < my_breaks[[1]]){
        node_colors = append(node_colors, colors[[1]])
      } else if(as.numeric(node_names[[i]]) > my_breaks[[2]]){
        node_colors = append(node_colors, colors[[3]])
      } else {
        node_colors = append(node_colors, colors[[2]])
      }
    }
  }
  
  # return colors
  return(node_colors)
}



### Read first N lines of a document -------------------------------------------
#' Read 'n' lines (ignoring comments and header) from a file.
#' 
#' Useful when you don't know the length/structure of a file
#' and want a useful sample to look at. Can skip ahead in the file too.
#' Copes well when there are less than 'n' lines in the file.
#'
#' @param fn name of the file(s) to get the length of
#' @param n number of valid lines to attempt to read
#'  looks at the top few lines (ignoring comments)
#' @param comment a comment symbol to ignore lines in files
#' @param skip number of lines to skip at top of file before processing
#' @param header whether to allow for, and skip, a header row
#' @return returns the first n lines of the file meeting the criteria,
#'  or if 'skip' implies lines beyond the length of the file, the 
#'  result,will be truncated - although in this case, the last 
#'  line will always be read.
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' orig.dir <- getwd(); setwd(tempdir()); # move to temporary dir
#' dat <- matrix(sample(100),nrow=10)
#' write.table(dat,"temp.txt",col.names=FALSE,row.names=FALSE)
#' n.readLines("temp.txt",n=2,skip=2,header=FALSE)
#' dat[3:4,]
#' unlink("temp.txt")
#' setwd(orig.dir) # reset working directory to original
n.readLines <- function(fn,n,comment="#",skip=0,header=TRUE)
{
  # read at least 'n' lines of a file, skipping lines and ignoring any starting with comment
  if(!file.exists(fn)) { warning("file doesn't exist"); return(NULL) }
  if(!is.character(comment)) { warning("illegal comment char, reverting to #"); comment <- "#" }
  rl <- 0; cc <- 0 + {if(is.numeric(skip)) skip else 0 }
  while(rl<n) { 
    test.bit <- readLines(fn,n+cc)
    if(skip>0 & length(test.bit>1)) { test.bit <- test.bit[-(1:(min((length(test.bit)-1),skip)))] }
    cmnt <- which(substr(test.bit,1,1)==comment)
    rl <- n+cc-length(cmnt)
    cc <- cc + length(cmnt)
  }
  if(length(cmnt)>0) { test.bit <- test.bit[-cmnt] } 
  if(length(test.bit)>1 & header) { test.bit <- test.bit[-1] }
  return(test.bit)
}



#==============================================================================#
# 2 - cp and mt trees comparison                                            ####

#------------------------------------------------------------------------------#
## 2.1 - import and sort trees                                              ####


# open printing file
pdf(file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".01_species_trees_comparisons.pdf", sep = ""),
    width = 12,
    height = 8)


### phylogenetic trees placeholder
phylogenetic_trees = list("cp_ML" = list("type" = "ML", "tree" = NA, "root" = NA, "dendro" = NA),
                          "cp_ASTRAL" = list("type" = "coalescent", "tree" = NA, "root" = NA, "dendro" = NA),
                          "mt_ML" = list("type" = "ML", "tree" = NA, "root" = NA, "dendro" = NA),
                          "mt_ASTRAL" = list("type" = "coalescent", "tree" = NA, "root" = NA, "dendro" = NA),
                          "cp_mt_ML" = list("type" = "ML", "tree" = NA, "root" = NA, "dendro" = NA),
                          "cp_mt_ASTRAL" = list("type" = "coalescent", "tree" = NA, "root" = NA, "dendro" = NA))


### import trees

# read ML concatenated trees
phylogenetic_trees[["cp_ML"]][["tree"]] = ape::read.tree(file = "04_cp_concat_ML/cp_allgenes_concat.contree")
phylogenetic_trees[["mt_ML"]][["tree"]] = ape::read.tree(file = "04_mt_concat_ML/mt_allgenes_concat.contree")
phylogenetic_trees[["cp_mt_ML"]][["tree"]] = ape::read.tree(file = "06_cp_mt_concat_ML/cp_mt_allgenes_concat.contree")

# read coalescence-based trees
phylogenetic_trees[["cp_ASTRAL"]][["tree"]] = ape::read.tree(file = "05_cp_ASTRAL/cp_genetrees.coalescence.tre")
phylogenetic_trees[["mt_ASTRAL"]][["tree"]] = ape::read.tree(file = "05_mt_ASTRAL/mt_genetrees.coalescence.tre")
phylogenetic_trees[["cp_mt_ASTRAL"]][["tree"]] = ape::read.tree(file = "07_cp_mt_ASTRAL/cp_mt_genetrees.coalescence.tre")


### reroot trees

for(i in 1:length(phylogenetic_trees)){
  # get root position
  phylogenetic_trees[[i]][["root"]] = tidytree::MRCA(phylogenetic_trees[[i]][["tree"]], "Pakin", "Oviri")
  # reroot set branch length NaN to 0 in coalescent trees and root length to 0.01 in ML trees
  if(phylogenetic_trees[[i]][["type"]] == "coalescent"){
    phylogenetic_trees[[i]][["tree"]] = treeio::root(phylogenetic_trees[[i]][["tree"]], node = phylogenetic_trees[[i]][["root"]], resolve.root = TRUE)
    phylogenetic_trees[[i]][["tree"]]$edge.length[is.na(phylogenetic_trees[[i]][["tree"]]$edge.length)] = 0
    phylogenetic_trees[[i]][["tree"]]$node.label = ifelse(phylogenetic_trees[[i]][["tree"]]$node.label == "Root", "1", phylogenetic_trees[[i]][["tree"]]$node.label)
    phylogenetic_trees[[i]][["tree"]]$node.label = ifelse(phylogenetic_trees[[i]][["tree"]]$node.label == "", "1", phylogenetic_trees[[i]][["tree"]]$node.label)
  } else {
    phylogenetic_trees[[i]][["tree"]] = treeio::root(phylogenetic_trees[[i]][["tree"]], node = phylogenetic_trees[[i]][["root"]] + 1, resolve.root = TRUE)
    phylogenetic_trees[[i]][["tree"]]$edge.length[which(phylogenetic_trees[[i]][["tree"]]$edge.length == 0)] = 0.1
    phylogenetic_trees[[i]][["tree"]]$node.label = ifelse(phylogenetic_trees[[i]][["tree"]]$node.label == "Root", "100", phylogenetic_trees[[i]][["tree"]]$node.label)
    phylogenetic_trees[[i]][["tree"]]$node.label = ifelse(phylogenetic_trees[[i]][["tree"]]$node.label == "", "100", phylogenetic_trees[[i]][["tree"]]$node.label)
  }
  # sort trees in descending order
  phylogenetic_trees[[i]][["tree"]] = ape::ladderize(phylogenetic_trees[[i]][["tree"]], right = TRUE)
  # set node label to numeric (BS and coalescence support)
  phylogenetic_trees[[i]][["tree"]]$node.label = as.numeric(phylogenetic_trees[[i]][["tree"]]$node.label)
  # plot tree
  plot(phylogenetic_trees[[i]][["tree"]])
  nodelabels(text = phylogenetic_trees[[i]][["tree"]]$node.label,
             node = 1:phylogenetic_trees[[i]][["tree"]]$Nnode + Ntip(phylogenetic_trees[[i]][["tree"]]),
             frame = "none", bg = "white")
}



#------------------------------------------------------------------------------#
## 2.2 - plot single trees                                                  ####

# ggtree(phylogenetic_trees[["cp_ML"]][["tree"]],
#        color = "black", size = 0.15, linetype = 1,  right = TRUE) + 
#   geom_tiplab(size = 4.5, align = TRUE) +
#   xlim(0, 5) + 
#   geom_point2(data = phylogenetic_trees[["cp_ML"]][["tree"]],
#               aes(subset = !isTip,
#                   fill = cut(node.label, c(0, 70, 90, 100))),
#               shape = 21, size = 4) +
#   scale_fill_manual(values = c("white"))
# 
#   scale_fill_manual(values = c("white", "grey", "black"), guide = "legend", 
#                     name = "Bootstrap Percentage (BP)", 
#                     breaks = c("(90,100]", "(70,90]", "(0,70]"))
# 
# 
#   scale_fill_manual(values = rep("white", 71))
#   
#   # geom_tiplab(size=4.5, hjust = -0.060, fontface = "bold") +  xlim(0, 0.09) + 
#   # geom_point2(aes(subset = !isTip & node != root,
#   #                 fill=cut(support, c(0, 70, 90, 100))),
#   #             shape=21, size=4) +
#   theme_tree(legend.position=c(0.2, 0.2)) + 
#   ggtitle("Example") +
#   scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
#                     name='Bootstrap Percentage(BP)', 
#                     breaks=c('(90,100]', '(70,90]', '(0,70]'), 
#                     labels=expression(BP>=90,70 <= BP * " < 90", BP < 70))


  
#------------------------------------------------------------------------------#
## 2.3 - cp_mt ML VS astral comparison                                      ####

### CP
# cp_ML dendro
cp_ML_dendro = as.dendrogram(chronos(phylogenetic_trees[["cp_ML"]][["tree"]]))
  
cp_ML_dendro = cp_ML_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(cp_ML_dendro, "ML")) %>%
  set("labels_colors", value = c(rep("black", 34), "red", "red")) %>% 
  set("leaves_cex", value = 0)

# cp_ASTRAL dendro  
cp_ASTRAL_dendro = as.dendrogram(chronos(phylogenetic_trees[["cp_ASTRAL"]][["tree"]]))

cp_ASTRAL_dendro = cp_ASTRAL_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(cp_ASTRAL_dendro, "coalescent")) %>%
  set("labels_colors", value = c(rep("black", 34), "red", "red")) %>% 
  set("leaves_cex", value = 0)

# save dendro
phylogenetic_trees[["cp_ML"]][["dendro"]] = cp_ML_dendro
phylogenetic_trees[["cp_ASTRAL"]][["dendro"]] = cp_ASTRAL_dendro

# plot
dendextend::tanglegram(cp_ML_dendro,
                       cp_ASTRAL_dendro,
                       color_lines = c(rep("grey", 4), rep("steelblue", 5),
                                       rep("grey", 18), rep("steelblue", 9)),
                       main_left = "ML",
                       main_right = "coalescence",
                       main = "cp, 70 genes",
                       common_subtrees_color_branches = FALSE,
                       fast = TRUE,
                       margin_bottom = 0,
                       margin_inner = 5,
                       margin_outer = 1.5,
                       axes = FALSE)

### MT
# mt_ML dendro
mt_ML_dendro = as.dendrogram(chronos(phylogenetic_trees[["mt_ML"]][["tree"]]))

mt_ML_dendro = mt_ML_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(mt_ML_dendro, "ML")) %>%
  set("labels_colors", value = c(rep("black", 37), "red", "red")) %>% 
  set("leaves_cex", value = 0)

mt_ASTRAL_dendro = as.dendrogram(chronos(phylogenetic_trees[["mt_ASTRAL"]][["tree"]]))

# mt_ASTRAL dendro
mt_ASTRAL_dendro = mt_ASTRAL_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(mt_ASTRAL_dendro, "coalescent")) %>%
  set("labels_colors", value = c(rep("black", 37), "red", "red")) %>% 
  set("leaves_cex", value = 0)

# save dendro
phylogenetic_trees[["mt_ML"]][["dendro"]] = mt_ML_dendro
phylogenetic_trees[["mt_ASTRAL"]][["dendro"]] = mt_ML_dendro

# plot
dendextend::tanglegram(mt_ML_dendro,
                       mt_ASTRAL_dendro,
                       color_lines = c(rep("grey", 11), rep("red", 3),
                                       rep("grey", 9), rep("red", 2),
                                       rep("grey", 3), "red"),
                       main_left = "ML",
                       main_right = "coalescence",
                       main = "mt, 29 genes",
                       common_subtrees_color_branches = TRUE,
                       fast = TRUE,
                       margin_bottom = 0,
                       margin_inner = 5,
                       margin_outer = 1.5,
                       axes = FALSE)

### CP + MT
# cp_mt_ML dendro
cp_mt_ML_dendro = as.dendrogram(chronos(phylogenetic_trees[["cp_mt_ML"]][["tree"]]))

cp_mt_ML_dendro = cp_mt_ML_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(cp_mt_ML_dendro, "ML")) %>%
  set("labels_colors", value = c(rep("black", 49), "red", "red")) %>% 
  set("leaves_cex", value = 0)

cp_mt_ASTRAL_dendro = as.dendrogram(chronos(phylogenetic_trees[["cp_mt_ASTRAL"]][["tree"]]))

# cp_mt_ASTRAL dendro
cp_mt_ASTRAL_dendro = cp_mt_ASTRAL_dendro %>%
  set("by_labels_branches_col", value = c("Oviri", "Pakin")) %>%
  set("by_labels_branches_lty", value = c("Oviri", "Pakin")) %>%
  set("nodes_cex", value = 1.75) %>% 
  set("nodes_pch", value = 21) %>% 
  set("nodes_bg", value = color_nodes(cp_mt_ASTRAL_dendro, "coalescent")) %>%
  set("labels_colors", value = c(rep("black", 49), "red", "red")) %>% 
  set("leaves_cex", value = 0)

# save dendro
phylogenetic_trees[["cp_mt_ML"]][["dendro"]] = cp_mt_ML_dendro
phylogenetic_trees[["cp_mt_ASTRAL"]][["dendro"]] = cp_mt_ASTRAL_dendro

# plot
dendextend::tanglegram(cp_mt_ML_dendro,
                       cp_mt_ASTRAL_dendro,
                       color_lines = c(rep("steelblue", 10), rep("grey", 4), rep("steelblue", 4),
                                       rep("grey", 4), "steelblue", rep("grey", 2),
                                       rep("steelblue", 4), rep("grey", 10), "steelblue",
                                       rep("grey", 5), rep("red", 2), rep("steelblue", 4)),
                       main_left = "ML",
                       main_right = "coalescence",
                       main = "cp + mt, 99 genes",
                       common_subtrees_color_branches = FALSE,
                       fast = TRUE,
                       margin_bottom = 0,
                       margin_inner = 5,
                       margin_outer = 1.5,
                       axes = FALSE)

### want to color the branches green???
#set("by_labels_branches_col", value = c(1,4))
#colours = colorRampPalette(brewer.pal(9, "Greens"))(11)



#------------------------------------------------------------------------------#
## 2.4 - calculate tree distances                                           ####

# get list of dendrograms to compare
dendro_list = dendlist("cp_ML" = cp_ML_dendro, "cp_ASTRAL" = cp_ASTRAL_dendro,
                       "mt_ML" = cp_ML_dendro, "mt_ASTRAL" = mt_ASTRAL_dendro,
                       "cp_mt_ML" = cp_ML_dendro, "cp_mt_ASTRAL" = cp_mt_ASTRAL_dendro)

# create empty lists for distances and correlations
empty_matrix = matrix(nrow = 6, ncol = 6)
colnames(empty_matrix) = c("cp_ML", "cp_ASTRAL", "mt_ML", "mt_ASTRAL", "cp_mt_ML", "cp_mt_ASTRAL")
rownames(empty_matrix) = c("cp_ML", "cp_ASTRAL", "mt_ML", "mt_ASTRAL", "cp_mt_ML", "cp_mt_ASTRAL")
dist_DB_pairwise = list("common_nodes_corr" = empty_matrix,
                        "RF" = empty_matrix,
                        "KF" = empty_matrix,
                        "MAST_%" = empty_matrix)
dist_DB_min_set_species = list("common_nodes_corr" = empty_matrix,
                               "RF" = empty_matrix,
                               "KF" = empty_matrix,
                               "MAST_%" = empty_matrix)

# common nodes correlation
common_nodes_corr = cor.dendlist(dendro_list, method = "common_nodes")
dist_DB_pairwise[["common_nodes_corr"]] = common_nodes_corr
dist_DB_min_set_species[["common_nodes_corr"]] = common_nodes_corr

# get pairwise distances
for(i in 1:length(phylogenetic_trees)){
  for(k in 1:length(phylogenetic_trees)){
    
    # get tree
    tree1_tmp = phylogenetic_trees[[i]][["tree"]]
    tree2_tmp = phylogenetic_trees[[k]][["tree"]]
    
    # list of common species
    species_list = tree1_tmp$tip.label[which(tree1_tmp$tip.label %in% tree2_tmp$tip.label)]

    # get overlapping species
    tree1 = ape::keep.tip(tree1_tmp, species_list)
    tree2 = ape::keep.tip(tree2_tmp, species_list)
    
    # get distances
    distances = phangorn::treedist(tree1, tree2)
    dist_DB_pairwise[["RF"]][i, k] = distances[[1]]
    dist_DB_pairwise[["KF"]][i, k] = distances[[2]]
    #dist_DB_pairwise[["path_diff"]][i, k] = distances[[3]]
    dist_DB_pairwise[["MAST_%"]][i, k] = length(phangorn::mast(tree1, tree2, tree = FALSE)) / length(species_list)
    
    # clean
    rm(tree1_tmp, tree2_tmp, tree1, tree2, species_list, distances)
  }
}

# make list of minimum set of species
min_species_list = phylogenetic_trees[["cp_mt_ML"]][["tree"]]$tip.label
for(i in 1:length(phylogenetic_trees)){
  min_species_list = min_species_list[which(min_species_list %in% phylogenetic_trees[[i]][["tree"]]$tip.label)]
}

# get distances based on minimum set of species
for(i in 1:length(phylogenetic_trees)){
  for(k in 1:length(phylogenetic_trees)){

    # get overlapping species
    tree1 = ape::keep.tip(phylogenetic_trees[[i]][["tree"]], min_species_list)
    tree2 = ape::keep.tip(phylogenetic_trees[[k]][["tree"]], min_species_list)
    
    # get distances
    distances = phangorn::treedist(tree1, tree2)
    dist_DB_min_set_species[["RF"]][i, k] = distances[[1]]
    dist_DB_min_set_species[["KF"]][i, k] = distances[[2]]
    #dist_DB_min_set_species[["path_diff"]][i, k] = distances[[3]]
    dist_DB_min_set_species[["MAST_%"]][i, k] = length(phangorn::mast(tree1, tree2, tree = FALSE)) / length(min_species_list)
    
    # clean
    rm(tree1, tree2, distances)
  }
}


### plots based on pairwise label sets between trees

# prepare plot layout
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))

for(i in 1:length(dist_DB_pairwise)){
  # parameters for correlation VS distances
  correlation = FALSE
  color_scale = list(rev(COL2("RdBu")), COL1("Blues"),
                     COL1("Blues"), COL1("YlOrRd"))
  tile_shape = "color"
  if(i == 1){
    correlation = TRUE
    tile_shape = "pie"
  }
  # plot
  corrplot(dist_DB_pairwise[[i]],
           title = names(dist_DB_pairwise)[[i]],
           tile_shape, 
           "upper", order = "original", 
           is.corr = correlation, col = color_scale[[i]],
           tl.srt = 45, tl.col = "grey15",
           mar = c(0, 0, 2, 0),
           cl.pos = "r", cl.align.text = "l")
}

# reset layout
par(mfrow = c(1, 1))


### plots based on minimum set of overlapping species between trees

# prepare plot layout
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))

for(i in 1:length(dist_DB_min_set_species)){
  # parameters for correlation VS distances
  correlation = FALSE
  color_scale = list(rev(COL2("RdBu")), COL1("Blues"),
                     COL1("Blues"), COL1("YlOrRd"))
  tile_shape = "color"
  if(i == 1){
    correlation = TRUE
    tile_shape = "pie"
  }
  # plot
  corrplot(dist_DB_min_set_species[[i]],
           title = names(dist_DB_min_set_species)[[i]],
           tile_shape, 
           "upper", order = "original", 
           is.corr = correlation, col = color_scale[[i]],
           tl.srt = 45, tl.col = "grey15",
           mar = c(0, 0, 2, 0),
           cl.pos = "r", cl.align.text = "l")
}

# reset layout
par(mfrow = c(1, 1))

# close the pdf file
dev.off()



#------------------------------------------------------------------------------#
## 2.5 - Figure 1: pairwise species                                         ####

# open printing file
pdf(file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".02_Figure1_pairwise.pdf", sep = ""),
    width = 8,
    height = 12)

# set layout
layout(matrix(c(1,1,2,3,3,
                1,1,2,3,3,
                4,5,6,7,8),
              nrow = 3, ncol = 5, byrow = TRUE))

# plot tanglegram
dendextend::tanglegram(cp_mt_ML_dendro,
                       cp_mt_ASTRAL_dendro,
                       color_lines = c(rep("steelblue", 10), rep("grey", 4), rep("steelblue", 4),
                                       rep("grey", 4), "steelblue", rep("grey", 2),
                                       rep("steelblue", 4), rep("grey", 10), "steelblue",
                                       rep("grey", 5), rep("red", 2), rep("steelblue", 4)),
                       main_left = "ML",
                       main_right = "coalescence",
                       main = "cp + mt, 99 genes",
                       common_subtrees_color_branches = FALSE,
                       fast = TRUE,
                       margin_bottom = 0,
                       margin_inner = 5,
                       margin_outer = 1.5,
                       axes = FALSE,
                       just_one = FALSE)

for(i in 1:length(dist_DB_pairwise)){
  # parameters for correlation VS distances
  correlation = FALSE
  color_scale = list(rev(COL2("RdBu")), COL1("Blues"),
                     COL1("Blues"), COL1("YlOrRd"))
  tile_shape = "color"
  if(i == 1){
    correlation = TRUE
    tile_shape = "pie"
  }
  # plot
  corrplot(dist_DB_pairwise[[i]],
           tile_shape, 
           "upper", order = "original", 
           is.corr = correlation, col = color_scale[[i]],
           tl.srt = 70, tl.cex = 0.5, tl.col = "grey15",
           cl.pos = "r", cl.align.text = "l")
  
}

# reset layout
par(mfrow = c(1, 1))

# close the pdf file
dev.off()



#------------------------------------------------------------------------------#
## 2.6 - Figure 1: common species only                                      ####

# open printing file
pdf(file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".03_Figure1_overlapping_species_only.pdf", sep = ""),
    width = 8,
    height = 12)

# set layout
layout(matrix(c(1,1,2,3,3,
                1,1,2,3,3,
                4,5,6,7,8),
              nrow = 3, ncol = 5, byrow = TRUE))

# plot tanglegram
dendextend::tanglegram(cp_mt_ML_dendro,
                       cp_mt_ASTRAL_dendro,
                       color_lines = c(rep("steelblue", 10), rep("grey", 4), rep("steelblue", 4),
                                       rep("grey", 4), "steelblue", rep("grey", 2),
                                       rep("steelblue", 4), rep("grey", 10), "steelblue",
                                       rep("grey", 5), rep("red", 2), rep("steelblue", 4)),
                       main_left = "ML",
                       main_right = "coalescence",
                       main = "cp + mt, 99 genes",
                       common_subtrees_color_branches = FALSE,
                       fast = TRUE,
                       margin_bottom = 0,
                       margin_inner = 5,
                       margin_outer = 1.5,
                       axes = FALSE,
                       just_one = FALSE)

for(i in 1:length(dist_DB_min_set_species)){
  # parameters for correlation VS distances
  correlation = FALSE
  color_scale = list(rev(COL2("RdBu")), COL1("Blues"),
                     COL1("Blues"), COL1("YlOrRd"))
  tile_shape = "color"
  if(i == 1){
    correlation = TRUE
    tile_shape = "pie"
  }
  # plot
  corrplot(dist_DB_min_set_species[[i]],
           tile_shape, 
           "upper", order = "original", 
           is.corr = correlation, col = color_scale[[i]],
           tl.srt = 70, tl.cex = 0.5, tl.col = "grey15",
           cl.pos = "r", cl.align.text = "l")
  
}

# reset layout
par(mfrow = c(1, 1))

# close the pdf file
dev.off()



#==============================================================================#
# 3 - single genes stats and selection                                      ####

#------------------------------------------------------------------------------#
## 3.1 - plot genes occupancies                                             ####

# open printing file
pdf(file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".04_Single_genes_stats.pdf", sep = ""),
    width = 12,
    height = 8)

# read tables
cp_mt_genes_occupancy = read.delim("cp_mt_genes_occupancy.tab")

# drop last line which is uninformative
cp_mt_genes_occupancy = head(cp_mt_genes_occupancy, -1)

# format headers
colnames(cp_mt_genes_occupancy) = str_remove_all(colnames(cp_mt_genes_occupancy), ".aln.fa")
colnames(cp_mt_genes_occupancy) = str_replace_all(colnames(cp_mt_genes_occupancy), "X..Sequences", "Species")

# melt and clean
cp_mt_genes_occupancy = melt(cp_mt_genes_occupancy)
cp_mt_genes_occupancy$organell = ifelse(str_detect(cp_mt_genes_occupancy$variable, "_cp_CDS_"), "Chloroplast", "Mitochondrion")
cp_mt_genes_occupancy$variable = str_remove_all(cp_mt_genes_occupancy$variable, "X02_cp_CDS_align.")
cp_mt_genes_occupancy$variable = str_remove_all(cp_mt_genes_occupancy$variable, "X02_mt_CDS_align.")
# set different number to assign different color to mitochondrial genes
cp_mt_genes_occupancy$value = ifelse(cp_mt_genes_occupancy$value == 1 & cp_mt_genes_occupancy$organell == "Mitochondrion",
                                        2, cp_mt_genes_occupancy$value)
cp_mt_genes_occupancy$value = as.character(cp_mt_genes_occupancy$value)

# plot
p1 = ggplot(cp_mt_genes_occupancy) +
  geom_tile(aes(x = variable, y = Species, fill = value), color = "grey75") +
  scale_fill_manual(values = c("white", "#469d89", "steelblue")) +
  #scale_fill_gradient(low = "white", high = "#469d89") +
  facet_grid(~organell, scales = "free", space = "free_x") +
  labs(#title = "Ulva organellar genes occupancy",
       x = "Genes",
       y = "Species") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        panel.background = element_blank(),
        #panel.spacing = unit(-0.2, "lines"),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_rect(colour = "NA", fill = "NA"),
        strip.placement = "outside",
        strip.text = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5))
plot(p1)


#------------------------------------------------------------------------------#
## 3.2 - plot gene alignments length distributions                          ####

# create gene distribution matrix
gene_length_distributions = matrix(ncol = 2)

# get cp length distributions
filenames = list.files("./02_cp_CDS_align/", pattern = "*.aln.fa", full.names = TRUE)
cp_genes_length = lapply(filenames, n.readLines, n = 1, skip = 1)
names(cp_genes_length) = str_remove_all(str_remove_all(filenames, "./02_cp_CDS_align//"), ".aln.fa")
for(i in 1:length(cp_genes_length)){
  gene_length_distributions = rbind(gene_length_distributions, c("cp", nchar(cp_genes_length[[i]][[1]])))
  rownames(gene_length_distributions)[nrow(gene_length_distributions)] = names(cp_genes_length)[i]
}

# get mt length distributions
filenames = list.files("./02_mt_CDS_align/", pattern = "*.aln.fa", full.names = TRUE)
mt_genes_length = lapply(filenames, n.readLines, n = 1, skip = 1)
names(mt_genes_length) = str_remove_all(str_remove_all(filenames, "./02_mt_CDS_align//"), ".aln.fa")
for(i in 1:length(mt_genes_length)){
  gene_length_distributions = rbind(gene_length_distributions, c("mt", nchar(mt_genes_length[[i]][[1]])))
  rownames(gene_length_distributions)[nrow(gene_length_distributions)] = names(mt_genes_length)[i]
}

# reformat table
gene_length_distributions = as.data.frame(gene_length_distributions)
gene_length_distributions = gene_length_distributions[-1, ]
names(gene_length_distributions) = c("Organell", "length (bp)")
gene_length_distributions$`length (bp)` = as.numeric(gene_length_distributions$`length (bp)`)

# plot
p2 = ggplot(gene_length_distributions, aes(x = Organell, y = `length (bp)`, fill = Organell)) +
  geom_point(aes(group = Organell, color = Organell), shape = "|", size = 5) +
  geom_boxplot(aes(fill = Organell), alpha = 0.75, width = 0.125,
               position = position_nudge(x = -0.1),
               outlier.colour = NA) +
  ggdist::stat_halfeye(aes(fill = Organell), alpha = 0.75, slab_color = "grey45", slab_size = 0.5, 
                       scale = 0.7, adjust = 1, justification = -0.05, .width = 0, point_colour = NA) +
  scale_color_manual(values = c("#469d89", "steelblue")) +
  scale_fill_manual(values = c("#469d89", "steelblue")) +
  scale_y_continuous(labels = label_number(suffix = " kbp", scale = 1e-3)) +
  coord_flip() +
  labs(#title = "Ulva organellar gene alignment length distribution",
       y = "Gene Length") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey85"),
        panel.grid.minor.x = element_line(colour = "grey85"),
        panel.border = element_rect(colour = "black", fill = NA))
plot(p2)

grid.arrange(p1, p2, heights = c(0.67, 0.33), ncol = 1)

# export gene length distributions table
write.table(gene_length_distributions, file = "./08_genes_distances/cp_mt_alignment_length_distributions.txt",
            quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)



#------------------------------------------------------------------------------#
## 3.3 - Genes distances                                                    ####

# get cp ML gene trees
cp_ML_genes = list()
filenames = list.files("./03_cp_singleGene_ML/", pattern = "*.contree", full.names = TRUE)
cp_ML_genes = lapply(filenames, ape::read.tree)
names(cp_ML_genes) = str_remove_all(str_remove_all(filenames, "./03_cp_singleGene_ML//"), ".aln.contree")
cp_ML_genes[["cp_mt_ML"]] = phylogenetic_trees[["cp_mt_ML"]][["tree"]]

# get mt ML gene trees
mt_ML_genes = list()
filenames = list.files("./03_mt_singleGene_ML/", pattern = "*.contree", full.names = TRUE)
mt_ML_genes = lapply(filenames, ape::read.tree)
names(mt_ML_genes) = str_remove_all(str_remove_all(filenames, "./03_mt_singleGene_ML//"), ".aln.contree")
mt_ML_genes[["cp_mt_ML"]] = phylogenetic_trees[["cp_mt_ML"]][["tree"]]

# create distance lists
empty_matrix = matrix(nrow = length(cp_ML_genes), ncol = length(cp_ML_genes))
colnames(empty_matrix) = names(cp_ML_genes)
rownames(empty_matrix) = names(cp_ML_genes)
cp_genes_distances = list("common_nodes_corr" = empty_matrix,
                          "RF" = empty_matrix,
                          "KF" = empty_matrix,
                          "MAST_%" = empty_matrix)

empty_matrix = matrix(nrow = length(mt_ML_genes), ncol = length(mt_ML_genes))
colnames(empty_matrix) = names(mt_ML_genes)
rownames(empty_matrix) = names(mt_ML_genes)
mt_genes_distances = list("common_nodes_corr" = empty_matrix,
                          "RF" = empty_matrix,
                          "KF" = empty_matrix,
                          "MAST_%" = empty_matrix)

empty_matrix = matrix(nrow = length(genes_sets[["cp_mt_ML"]]), ncol = length(genes_sets[["cp_mt_ML"]]))
colnames(empty_matrix) = names(genes_sets[["cp_mt_ML"]])
rownames(empty_matrix) = names(genes_sets[["cp_mt_ML"]])
cp_mt_genes_distances = list("common_nodes_corr" = empty_matrix,
                             "RF" = empty_matrix,
                             "KF" = empty_matrix,
                             "MAST_%" = empty_matrix)


# populate distances lists
genes_sets = list("cp_ML" = cp_ML_genes,
                  "mt_ML" = mt_ML_genes,
                  "cp_mt_ML" = c(cp_ML_genes, mt_ML_genes[which(names(mt_ML_genes) != "cp_mt_ML")]))
distance_sets = list("cp_genes_distances" = cp_genes_distances,
                     "mt_genes_distances" = mt_genes_distances,
                     "cp_mt_genes_distances" = cp_mt_genes_distances)

# get pairwise distances
for(j in 1:length(distance_sets)){
  for(i in 1:length(genes_sets[[j]])){
    for(k in 1:length(genes_sets[[j]])){
      
      # error handling
      # MAST call errors for very unresulved trees (e.g.: psaM)
      #   Error in root.phylo(x, bipart_x, resolve.root = TRUE) : 
      #   the specified outgroup is not monophyletic
      # catch error and move forward
      skip_to_next = FALSE
      
      # get tree
      tree1_tmp = genes_sets[[j]][[i]]
      tree2_tmp = genes_sets[[j]][[k]]
      
      # list of common species
      species_list = tree1_tmp$tip.label[which(tree1_tmp$tip.label %in% tree2_tmp$tip.label)]
      
      # get overlapping species
      tree1 = ape::keep.tip(tree1_tmp, species_list)
      tree2 = ape::keep.tip(tree2_tmp, species_list)
      
      # get distances
      distances = phangorn::treedist(tree1, tree2)
      distance_sets[[j]][["RF"]][i, k] = distances[[1]]
      distance_sets[[j]][["KF"]][i, k] = distances[[2]]
      #dist_DB_pairwise[["path_diff"]][i, k] = distances[[3]]
      distance_sets[[j]][["MAST_%"]][i, k] = tryCatch(length(phangorn::mast(tree1, tree2, tree = FALSE)) / length(species_list),
                                                      error = function(e) { skip_to_next <<- TRUE })
      
      # error handling
      if(skip_to_next){distance_sets[[j]][["MAST_%"]][i, k] = NA}
    }
  }
}


### get pairwise common_nodes_corr distances

# create CPU cluster
cl = parallel::makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

# get chronos in parallel
chronos_list = foreach(k = 1:length(genes_sets[["cp_mt_ML"]]), .combine = "c") %dopar% {
  library("ape")
  # get tree
  tryCatch(chronos(genes_sets[["cp_mt_ML"]][[k]]), error = function(e) { NA })
}  
parallel::stopCluster(cl)
names(chronos_list) = names(genes_sets[["cp_mt_ML"]])

# get pairwise common_nodes_corr distances
for(j in 1:length(distance_sets)){
  for(i in 1:length(genes_sets[[j]])){
    for(k in 1:length(genes_sets[[j]])){
      
      # get trees
      tree1_tmp = chronos_list[[names(genes_sets[[j]])[i]]]
      tree2_tmp = chronos_list[[names(genes_sets[[j]])[k]]]
      
      # check trees and get distance
      if(all(!is.na(tree1_tmp)) & all(!is.na(tree2_tmp))){
        distance_sets[[j]][["common_nodes_corr"]][i, k] = cor.dendlist(
          dendlist(as.dendrogram(tree1_tmp), as.dendrogram(tree2_tmp)),
          method = "common_nodes")[2]
      } else {
        distance_sets[[j]][["common_nodes_corr"]][i, k] = NA
      }
    }
  }
}


# rename matrices for cp/mt overlapping rpl/rps genes (for plotting reasons)
for(i in 1:length(distance_sets[["mt_genes_distances"]])){
  rownames(distance_sets[["mt_genes_distances"]][[i]]) = c("atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", "nad2",
                                        "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "rpl14 ", "rpl16 ", "rpl5 ", "rps10",
                                        "rps11 ", "rps12 ", "rps13", "rps14 ", "rps19 ", "rps2 ", "rps3 ", "rps4 ", "cp_mt_ML")
  colnames(distance_sets[["mt_genes_distances"]][[i]]) = c("atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", "nad2",
                                        "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "rpl14 ", "rpl16 ", "rpl5 ", "rps10",
                                        "rps11 ", "rps12 ", "rps13", "rps14 ", "rps19 ", "rps2 ", "rps3 ", "rps4 ", "cp_mt_ML")
}
for(i in 1:length(distance_sets[["cp_mt_genes_distances"]])){
  rownames(distance_sets[["cp_mt_genes_distances"]][[i]]) = c("accD", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI", "ccsA", "cemA", "chlI", "clpP",
                                           "infA", "petA", "petB", "petD", "petG", "petL", "psaA", "psaB", "psaC", "psaI", "psaJ",
                                           "psaM", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK",
                                           "psbL", "psbM", "psbN", "psbT", "psbZ", "rbcL", "rpl12", "rpl14", "rpl16", "rpl19", "rpl2",
                                           "rpl20", "rpl23", "rpl32", "rpl36", "rpl5", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps11",
                                           "rps12", "rps14", "rps18", "rps19", "rps2", "rps3",  "rps4", "rps7", "rps8", "rps9",
                                           "tufA", "ycf1", "ycf12", "ycf20", "ycf3", "ycf4", "cp_mt_ML", "atp1", "atp4", "atp6",
                                           "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4L",
                                           "nad5", "nad6", "nad7", "rpl14 ", "rpl16 ", "rpl5 ", "rps10", "rps11 ", "rps12 ", "rps13",
                                           "rps14 ", "rps19 ", "rps2 ", "rps3 ", "rps4 ")
  colnames(distance_sets[["cp_mt_genes_distances"]][[i]]) = c("accD", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI", "ccsA", "cemA", "chlI", "clpP",
                                           "infA", "petA", "petB", "petD", "petG", "petL", "psaA", "psaB", "psaC", "psaI", "psaJ",
                                           "psaM", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK",
                                           "psbL", "psbM", "psbN", "psbT", "psbZ", "rbcL", "rpl12", "rpl14", "rpl16", "rpl19", "rpl2",
                                           "rpl20", "rpl23", "rpl32", "rpl36", "rpl5", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps11",
                                           "rps12", "rps14", "rps18", "rps19", "rps2", "rps3",  "rps4", "rps7", "rps8", "rps9",
                                           "tufA", "ycf1", "ycf12", "ycf20", "ycf3", "ycf4", "cp_mt_ML", "atp1", "atp4", "atp6",
                                           "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", "nad2", "nad3", "nad4", "nad4L",
                                           "nad5", "nad6", "nad7", "rpl14 ", "rpl16 ", "rpl5 ", "rps10", "rps11 ", "rps12 ", "rps13",
                                           "rps14 ", "rps19 ", "rps2 ", "rps3 ", "rps4 ")
}

### plot heatmaps of distances
# prepare gene lists
cp_genes_list = c("accD", "atpA", "atpB", "atpE", "atpF", "atpH", "atpI", "ccsA", "cemA", "chlI", "clpP",
                  "infA", "petA", "petB", "petD", "petG", "petL", "psaA", "psaB", "psaC", "psaI", "psaJ",
                  "psaM", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK",
                  "psbL", "psbM", "psbN", "psbT", "psbZ", "rbcL", "rpl12", "rpl14", "rpl16", "rpl19", "rpl2",
                  "rpl20", "rpl23", "rpl32", "rpl36", "rpl5", "rpoA", "rpoB", "rpoC1", "rpoC2", "rps11",
                  "rps12", "rps14", "rps18", "rps19", "rps2", "rps3",  "rps4", "rps7", "rps8", "rps9",
                  "tufA", "ycf1", "ycf12", "ycf20", "ycf3", "ycf4")
mt_genes_list = c("atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3", "nad1", "nad2",
                  "nad3", "nad4", "nad4L", "nad5", "nad6", "nad7", "rpl14 ", "rpl16 ", "rpl5 ", "rps10",
                  "rps11 ", "rps12 ", "rps13", "rps14 ", "rps19 ", "rps2 ", "rps3 ", "rps4 ")

# prepare plot lists
distances_plots = list("cp" = list(), "mt" = list(), "cp_mt" = list())

# iterate plots
for(i in 1:length(distance_sets)){
  for(k in 1:length(distance_sets[[i]])){
    
    # remove possible rows/columns with only NAs
    tmp_matrix = distance_sets[[i]][[k]]
    tmp_matrix = tmp_matrix[, colSums(is.na(tmp_matrix)) < nrow(tmp_matrix)]
    tmp_matrix = tmp_matrix[rowSums(is.na(tmp_matrix)) < ncol(tmp_matrix), ]
    
    # get names
    plot_name = ifelse(names(distance_sets[[i]])[[k]] == "common_nodes_corr", "common\nnodes\ncorr", names(distance_sets[[i]])[[k]])
    
    # prep annotation
    organell = as.data.frame(ifelse(rownames(tmp_matrix) %in% cp_genes_list, "cp",
                                    ifelse(rownames(tmp_matrix) %in% mt_genes_list, "mt", "cp_mt")))
    rownames(organell) = rownames(tmp_matrix)
    colnames(organell) = "organell"
    organell$color = ifelse(organell$organell == "cp", "#469d89", ifelse(organell$organell == "mt", "steelblue", "red"))
    
    # plot raw distance
    distances_plots[[i]][[k]] = grid.grabExpr(draw(ComplexHeatmap::Heatmap(tmp_matrix,
                                                                           col = colorRamp2(c(min(tmp_matrix[!is.na(tmp_matrix)]),
                                                                                              max(tmp_matrix[!is.na(tmp_matrix)])),
                                                                                            c("white", "#D55E00")),
                                                                           row_names_gp = gpar(col = organell$color),
                                                                           column_names_gp = gpar(col = organell$color),
                                                                           top_annotation = HeatmapAnnotation(organell = as.matrix(organell$organell),
                                                                                                              show_annotation_name = FALSE,
                                                                                                              show_legend = FALSE,
                                                                                                              col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
                                                                           left_annotation = rowAnnotation(organell = as.matrix(organell$organell),
                                                                                                           show_annotation_name = FALSE,
                                                                                                           col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
                                                                           name = plot_name)))
    grid.newpage(recording = TRUE)
    grid.draw(distances_plots[[i]][[k]])

    # normalize
    normalized = t(scale(t(tmp_matrix)))
    # plot normalized distances
    distances_plots[[i]][[k + 4]] = grid.grabExpr(draw(ComplexHeatmap::Heatmap(normalized,
                                                                               col = colorRamp2(c(min(normalized[!is.na(normalized)]), 0, max(normalized[!is.na(normalized)])),
                                                                                                c("#56B4E9", "white", "#D55E00")),
                                                                               row_names_gp = gpar(col = organell$color),
                                                                               column_names_gp = gpar(col = organell$color),
                                                                               top_annotation = HeatmapAnnotation(organell = as.matrix(organell$organell),
                                                                                                                  show_annotation_name = FALSE,
                                                                                                                  show_legend = FALSE,
                                                                                                                  col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
                                                                               left_annotation = rowAnnotation(organell = as.matrix(organell$organell),
                                                                                                               show_annotation_name = FALSE,
                                                                                                               col = list(organell = c("cp_mt" = "red", "cp" = "#469d89", "mt" = "steelblue"))),
                                                                               name = paste(plot_name, "\nz-score", sep = ""))))
    grid.newpage(recording = TRUE)
    grid.draw(distances_plots[[i]][[k + 4]])
  }
}

# composite plots
for(i in 1:length(distances_plots)){
  grid.arrange(distances_plots[[i]][[1]], distances_plots[[i]][[2]],
               distances_plots[[i]][[3]], distances_plots[[i]][[4]],
               ncol = 2,
               top = textGrob(paste(names(distances_plots)[i], " raw distance", sep = "")))
  grid.arrange(distances_plots[[i]][[5]], distances_plots[[i]][[6]],
               distances_plots[[i]][[7]], distances_plots[[i]][[8]],
               ncol = 2,
               top = textGrob(paste(names(distances_plots)[i], " normalized distance", sep = "")))
  grid.arrange(distances_plots[[i]][[1]], distances_plots[[i]][[2]], distances_plots[[i]][[3]], distances_plots[[i]][[4]],
               distances_plots[[i]][[5]], distances_plots[[i]][[6]], distances_plots[[i]][[7]], distances_plots[[i]][[8]],
               ncol = 4,
               top = textGrob(paste(names(distances_plots)[i], " raw and normalized distance", sep = "")))
}

# export tables
for(i in 1:length(distance_sets)){
  for(k in 1:length(distance_sets[[i]])){
    write.table(distance_sets[[i]][[k]],
                file = paste("./08_genes_distances/", names(distance_sets)[i], ".distances.", names(distance_sets[[i]])[k], ".txt", sep = ""),
                quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  }
}

# clean up
cp_genes_distances = distance_sets[["cp_genes_distances"]]
mt_genes_distances = distance_sets[["mt_genes_distances"]]
cp_mt_genes_distances = distance_sets[["cp_mt_genes_distances"]]
rm(tree1_tmp, tree2_tmp, tree1, tree2, species_list, distances,
   distance_sets, chronos_list, plot_name, tmp_matrix, genes_sets, skip_to_next)

# close the pdf file
dev.off()

# export distance tables
for(i in 1:length(cp_genes_distances)){
  write.table(cp_genes_distances[[i]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".04_Single_genes_stats.cp.", names(cp_genes_distances)[i], ".txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
}
for(i in 1:length(mt_genes_distances)){
  write.table(mt_genes_distances[[i]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".04_Single_genes_stats.mt.", names(mt_genes_distances)[i], ".txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
}
for(i in 1:length(cp_mt_genes_distances)){
  write.table(cp_mt_genes_distances[[i]],
              file = paste("UlvaOmics.", format(Sys.Date(), format = "%Y%m%d"), ".04_Single_genes_stats.cp_mt.", names(cp_mt_genes_distances)[i], ".txt",  sep = ""),
              quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
}




#==============================================================================#
# 99 - Save R image file                                                    ####

save.image(file = "01_species_tree_viz.RData")
#load("01_species_tree_viz.RData")
