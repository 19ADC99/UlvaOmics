

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
