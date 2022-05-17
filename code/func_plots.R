### gene_filtering/code/func_plots.R
## R functions to facilitate nice plotting
# Caitlin Cherryh 2022

library(phytools) # Functions: nodeHeights

# Given a dataframe containing columns with treelikeness test statistic/statistical test values, this function creates a new column
#     for each test statistic/statistical test stating whether each value is treelike or non-treelike
classify.treelikeness.statistics <- function(df, tree_proportion_threshold){
  # Classify loci using p-values as cut off <- significant p-value = treelike
  df$X3SEQ_treelike <- as.numeric(df$X3SEQ_p_value)
  df$X3SEQ_treelike[df$X3SEQ_p_value <= 0.05] <- "TREELIKE"
  df$X3SEQ_treelike[df$X3SEQ_p_value > 0.05] <- "NON-TREELIKE"
  df$X3SEQ_treelike <- factor(df$X3SEQ_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  df$tree_proportion_p_value_treelike <- as.numeric(df$tree_proportion_p_value)
  df$tree_proportion_p_value_treelike[df$tree_proportion_p_value <= 0.05] <- "TREELIKE"
  df$tree_proportion_p_value_treelike[df$tree_proportion_p_value > 0.05] <- "NON-TREELIKE"
  df$tree_proportion_p_value_treelike <- factor(df$tree_proportion_p_value_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  
  # For Vanderpool data:
  # When cut-off is 0.7, 972/1730 are treelike
  # When cut-off is 0.75, 580/1730 are treelike
  # Median tree proportion is 0.7116612
  df$tree_proportion_treelike <- as.numeric(df$tree_proportion)
  df$tree_proportion_treelike[df$tree_proportion > tree_proportion_threshold] <- "TREELIKE"
  df$tree_proportion_treelike[df$tree_proportion <= tree_proportion_threshold] <- "NON-TREELIKE"
  df$tree_proportion_treelike <- factor(df$tree_proportion_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  
  # Sort loci by p-values for both tree proportion and 3seq
  df$sorted_p_value <- df$loci
  df$sorted_p_value[df$tree_proportion_p_value_treelike == "TREELIKE" & df$X3SEQ_treelike == "TREELIKE"] <- "Both"
  df$sorted_p_value[df$tree_proportion_p_value_treelike == "NON-TREELIKE" & df$X3SEQ_treelike == "TREELIKE"] <- "3seq only"
  df$sorted_p_value[df$tree_proportion_p_value_treelike == "TREELIKE" & df$X3SEQ_treelike == "NON-TREELIKE"] <- "Tree proportion only"
  df$sorted_p_value[df$tree_proportion_p_value_treelike == "NON-TREELIKE" & df$X3SEQ_treelike == "NON-TREELIKE"] <- "Neither"
  df$sorted_p_value <- factor(df$sorted_p_value, levels = c("Neither","3seq only","Tree proportion only","Both"), ordered = TRUE)
  
  # Sort loci by p-values for 3seq and test statistic value for tree proportion
  df$sorted <- df$loci
  df$sorted[df$tree_proportion_treelike == "TREELIKE" & df$X3SEQ_treelike == "TREELIKE"] <- "Both"
  df$sorted[df$tree_proportion_treelike == "NON-TREELIKE" & df$X3SEQ_treelike == "TREELIKE"] <- "3seq only"
  df$sorted[df$tree_proportion_treelike == "TREELIKE" & df$X3SEQ_treelike == "NON-TREELIKE"] <- "Tree proportion only"
  df$sorted[df$tree_proportion_treelike == "NON-TREELIKE" & df$X3SEQ_treelike == "NON-TREELIKE"] <- "Neither"
  df$sorted <- factor(df$sorted, levels = c("Neither","3seq only","Tree proportion only","Both"), ordered = TRUE)
  
  return(df)
}




# Open an alignment and return the number of phylogenetically informative sites and number of segregating sites in that alignment
get.DNA.site.info <- function(loci_name, alignment_dir){
  all_als <- list.files(alignment_dir)
  aln <- paste0(alignment_dir, grep(loci_name, all_als, value = TRUE))
  f <- read.dna(aln, format = "fasta")
  ss <- seg.sites(f)
  n_ss <- length(ss)
  n_pis <- ips::pis(f, what = "absolute")
  n_sites <- dim(f)[2]
  op <- c(n_ss, n_pis, n_sites)
  names(op) <- c("n_seg_sites", "n_pis", "n_sites")
  return(op)
}




# Rescale the total length (i.e. height) of a tree
rescale.tree.length <- function(tree, scaled_length){
  tree$edge.length <- tree$edge.length / max(nodeHeights(tree)[,2]) * scaled_length
  return(tree)
}


# Wrapper for feeding list of trees into rescale.tree.length using lapply
rescale.multiphylo <- function(trees, scaled_length){
  for (i in 1:length(trees)){
    trees[[i]] <- rescale.tree.length(trees[[i]], scaled_length)
  }
  return(trees)
}




