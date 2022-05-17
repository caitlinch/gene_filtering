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



# Quick function to take in vector of tomatoes species numbers and return 
rename.tomato.tips <- function(tip_vector){
  new_tip_names <- unlist(lapply(tip_vector, find.one.tomato.tip, include.genus = FALSE))
  return(new_tip_names)
}

# Quick function to take in a single species number from the tomatoes tree and return full scientific name for that species
find.one.tomato.tip <- function(single_tip, include.genus = FALSE){
  tip_list <- list("SL2.50" = c("Solanum.lyco. 'Heinz 1706'", "S. lyco. 'Heinz 1706'", "Solanum.lycopersicum 'Heinz 1706'", "S. lycopersicum 'Heinz 1706'"),
                   "LA3475" = c("Solanum lycopersicum 3475", "S. lycopersicum 3475"),
                   "LA3124" = c("Solanum cheesmaniae 3124", "S. cheesmaniae 3124"),
                   "LA0429" = c("Solanum galapagense 0429", "S. galapagense 0429"),
                   "LA0436" = c("Solanum cheesmaniae 0436", "S. cheesmaniae 0436"),
                   "LA3909" = c("Solanum galapagense 3909", "S. galapagense 3909"),
                   "LA2933" = c("Solanum lycopersicum 2933", "S. lycopersicum 2933", "Solanum lycopersicum var. cerasiforme 2933", "S. lycopersicum var. cerasiforme 2933"), 
                   "LA1269" = c("Solanum pimpinellifolium 1269", "S. pimpinellifolium 1269"),
                   "LA1589" = c("Solanum pimpinellifolium 1589", "S. pimpinellifolium 1589"),
                   "LA1028" = c("Solanum chmielewskii 1028", "S. chmielewskii 1028"),
                   "LA1316" = c("Solanum chmielewskii 1316", "S. chmielewskii 1316"),
                   "LA2172" = c("Solanum arcanum 2172", "S. arcanum 2172"),
                   "LA1322" = c("Solanum neorickii 1322", "S. neorickii 1322"), 
                   "LA2133" = c("Solanum neorickii 2133", "S. neorickii 2133"),
                   "LA4116" = c("Solanum sitiens 4116", "S. sitiens 4116"),
                   "LA2951" = c("Solanum lycopersicoides 2951", "S. lycopersicoides 2951"),
                   "LA4126" = c("Solanum lycopersicoides 4126", "S. lycopersicoides 4126"),
                   "LA0716" = c("Solanum pennellii 0716", "S. pennellii 0716"),
                   "LA3778" = c("Solanum pennellii 3778", "S. pennellii 3778"),
                   "LA0407" = c("Solanum habrochaites 0407", "S. habrochaites 0407"),
                   "LA1777" = c("Solanum habrochaites 1777", "S. habrochaites 1777"),
                   "LA1364" = c("Solanum huaylasense 1364", "S. huaylasense 1364"),
                   "LA1782" = c("Solanum chilense 1782", "S. chilense 1782"),
                   "LA4117" = c("Solanum chilense 4117A", "S. chilense 4117A"),
                   "LA2744" = c("Solanum peruvianum 2744", "S. peruvianum 2744"),
                   "LA2964" = c("Solanum peruvianum 2964", "S. peruvianum 2964"),
                   "LA0107" = c("Solanum corneliomuelleri 0107", "S. corneliomuelleri 0107"),
                   "LA0444" = c("Solanum corneliomuelleri 0444", "S. corneliomuelleri 0444"), 
                   "LA1358" = c("Solanum huaylasense 1358", "S. huaylasense 1358"),
                   "LA1360" = c("Solanum huaylasense 1360", "S. huaylasense 1360"))
  if (single_tip %in% names(tip_list)){
    # If this taxa is included in the tip list, use the taxa number to look up the taxa species name
    if (include.genus == FALSE){
      new_tip_name <- tip_list[[single_tip]][2]
    } else if (include.genus == TRUE){
      new_tip_name <- tip_list[[single_tip]][1]
    }
  } else {
    # If this taxa is not included in the tip list it may be because you have already renamed it. Return the name as is.
    new_tip_name <- single_tip
  }
  return(new_tip_name)
}



# Quick function to remove the sections of the tomato tree that remain identical
reformat.congruent.clades <- function(tom_tree){
  # Want to keep all Peruvianum species:
  Peruvianum = c("LA1364", "LA2744", "LA1358", "LA0107", "LA0444", "LA2964", "LA1782", "LA4117")
  # Want to keep one tip from each other clade:
  Esculentum = "SL2.50"
  Arcanum = "LA1322"
  Hirsutum = "LA0407"
  Clade_Outgroup = c("LA4116")
  root_outgroup = c("LA4116", "LA4126", "LA2951")
  # Collate the list of taxa to keep
  taxa_to_keep <- c(Peruvianum, Esculentum, Arcanum, Hirsutum, Clade_Outgroup)
  # Root at outgroup
  tom_tree <- root(tom_tree, root_outgroup)
  # Drop all other taxa
  tom_tree <- keep.tip(tom_tree, taxa_to_keep)
  # Rename the remaining taxa into the clade name
  tom_tree$tip.label[which(tom_tree$tip.label == "SL2.50")] <- "Esculentum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA1322")] <- "Arcanum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA0407")] <- "Hirsutum clade"
  tom_tree$tip.label[which(tom_tree$tip.label == "LA4116")] <- "Outgroup"
  # Return the reformatted tree
  return(tom_tree)
}





