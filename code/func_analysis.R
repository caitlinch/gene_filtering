### empirical_treelikeness/code/func_analysis.R
## R functions to analyse empirical treelikeness data
# Caitlin Cherryh 2021

library(ape)
library(phangorn)

# Given an alignment, this function will check the IQ-Tree .log file and .iqtree file for warnings and return all warnings and the corresponding alignment
check.for.IQTree.warnings <- function(alignment_path){
  # Open the .iqtree and .log files
  dotiqtree_path <- paste0(alignment_path, ".iqtree")
  dotiqtree_file <- readLines(dotiqtree_path)
  dotlog_path <- paste0(alignment_path, ".log")
  dotlog_file <- readLines(dotlog_path)
  
  # Check for any warnings
  dotiqtree_warnings <- dotiqtree_file[grep("WARNING", dotiqtree_file)]
  dotlog_warnings <- dotlog_file[grep("WARNING", dotlog_file)]
  
  # Extract loci name from alignment
  alignment_path_base <- basename(alignment_path)
  alignment_path_split <- strsplit(alignment_path_base, "\\.")[[1]]
  loci_name <- paste0(alignment_path_split[1:(length(alignment_path_split) - 1)], collapse = ".")
  
  # If one or more warnings was found, output a dataframe of the warnings
  if ((length(dotiqtree_warnings) + length(dotlog_warnings)) > 0){
    # Output warnings
    warning_df <- data.frame(loci_path = alignment_path, loci = loci_name, file = c(rep(".iqtree", length(dotiqtree_warnings)), rep(".log", length(dotlog_warnings))),
                             warnings = c(dotiqtree_warnings, dotlog_warnings))
    return(warning_df)
  }
}



# This function looks in a single IQ-Tree folder and checks which loci are present
get.loci.from.analysis <- function(folder, output_folder){
  # Create a filename for the output file
  op_filename <- paste0(output_folder, gsub("_IQ-Tree_partition","",basename(folder)), "_loci.csv")
  # Get all files in the folder
  folder_files <- list.files(folder)
  # Remove any partition files from the folder
  nonpartition_files <- grep("partition", folder_files, invert = TRUE, value = TRUE)
  # Get the location of each alignment
  loci_files <- paste0(folder, "/", nonpartition_files)
  # Get the loci name from each alignment
  loci_names <- gsub(".fa", "", nonpartition_files)
  # Combine the names and locations of each alignment
  loci_df <- data.frame(loci_name = loci_names, loci_file = loci_files)
  # Output this information about this analysis
  write.csv(loci_df, op_filename, row.names = FALSE)
}



# Test params for writing
species_tree_csv <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/Vanderpool2020/all_species_trees/p-value_categories_none_50loci_loci.csv"
sample_size = 100
loci_tree_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/Vanderpool2020/loci_trees/"
output_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_results/Vanderpool2020/" 

# This function takes a list of names of loci, randomly selects 100 pairs, compares those two trees and records the average difference between them
calculate.average.tree.distance <- function(species_tree_csv, loci_tree_folder, output_folder, sample_size = 100){
  # Get the record of which loci were included in this category
  species_tree_df <- read.csv(species_tree_csv, stringsAsFactors = FALSE)
  # Randomly sample loci to be compared
  tree1_loci <- sample(species_tree_df$loci_name, sample_size, replace = TRUE)
  tree2_loci <- sample(species_tree_df$loci_name, sample_size, replace = TRUE)
  # Find the path to each loci tree
  tree1_loci_paths <- unlist(lapply(tree1_loci, find.loci.tree.path, loci_tree_folder))
  tree2_loci_paths <- unlist(lapply(tree2_loci, find.loci.tree.path, loci_tree_folder))
  # Create a dataframe of all this information
  compare_df <- data.frame(tree1_loci, tree1_loci_paths, tree2_loci, tree2_loci_paths)
  # Feed into the wrapper program row by row
  distances_list <- lapply(1:nrow(compare_df), compare.trees.wrapper, compare_df)
  distances_df <- as.data.frame(do.call(rbind, distances_list))
  # Save distances df
  distances_op_name <- paste0(output_folder, gsub("_loci.csv", paste0("_",sample_size,"Samples_TreeDistances.csv"), basename(species_tree_csv)))
  write.csv(distances_df, distances_op_name, row.names = FALSE)
  # Calculate the mean for each of the distances
  avg_dists <- c(mean(as.numeric(distances_df$RF_distance)), mean(as.numeric(distances_df$normalized_RF_distance)), mean(as.numeric(distances_df$weighted_RF_distance)),
                 mean(as.numeric(distances_df$normalized_weighted_RF_distance)), mean(as.numeric(distances_df$path_difference)), 
                 mean(as.numeric(distances_df$weighted_path_difference)), mean(as.numeric(distances_df$KF_distance)))
  avg_dists <- signif(avg_dists, digits = 4)
  names(avg_dists) <- c("mean_RF_distance", "mean_normalized_RF_distance", "mean_weighted_RF_distance",
                        "mean_normalized_weighted_RF_distance", "mean_path_difference", 
                        "mean_weighted_path_difference","mean_KF_distance")
  # Return average distances
  return(avg_dists)
}



# Wrapper function to take one row containing two loci and the path to each loci and return a row of distances between the loci
compare.trees.wrapper <- function(i, df){
  row <- df[i,]
  # Extract information from row
  l1_name <- row$tree1_loci
  l2_name <- row$tree2_loci
  l1_path <- row$tree1_loci_paths
  l2_path <- row$tree2_loci_paths
  # Apply to distance function
  dists <- compare.two.trees(l1_path, l2_path)
  new_row <- c(l1_name, l1_path, l2_name, l2_path, dists)
  names(new_row) <- c("loci1_name", "loci1_path", "loci2_name", "loci2_path", "RF_distance", "normalized_RF_distance", "weighted_RF_distance",
                      "normalized_weighted_RF_distance", "path_difference", "weighted_path_difference","KF_distance")
  return(new_row)
}



# This function takes two trees, calculates the difference between them and returns it
compare.two.trees <- function(tree1_path, tree2_path){
  # Open the two trees
  tree1 <- read.tree(file = tree1_path)
  tree2 <- read.tree(file = tree2_path)
  # Calculate differences
  dist_RF <- signif(RF.dist(tree1, tree2, check.labels = TRUE), digits = 5)
  dist_nRF <- signif(RF.dist(tree1, tree2, check.labels = TRUE, normalize = TRUE), digits = 5)
  dist_wRF <- signif(wRF.dist(tree1, tree2, check.labels = TRUE), digits = 5)
  dist_nwRF <- signif(wRF.dist(tree1, tree2, check.labels = TRUE, normalize = TRUE), digits = 5)
  dist_pathDiff <- signif(path.dist(tree1, tree2, check.labels = TRUE), digits = 5)
  dist_wPathDiff <- signif(path.dist(tree1, tree2, check.labels = TRUE, use.weight = TRUE), digits = 5)
  dist_KF <- signif(KF.dist(tree1, tree2, check.labels = TRUE), digits = 5)
  # Return as a named vector
  dists <- c(dist_RF, dist_nRF, dist_wRF, dist_nwRF, dist_pathDiff, dist_wPathDiff, dist_KF)
  names(dists) <- c("RF_distance", "normalized_RF_distance", "weighted_RF_distance",
                    "normalized_weighted_RF_distance", "path_difference", "weighted_path_difference","KF_distance")
  return(dists)
}



# Function to take one loci name and select the loci tree from a folder of loci trees
find.loci.tree.path <- function(loci_name, loci_tree_folder){
  all_trees <- list.files(loci_tree_folder, full.names = TRUE)
  loci_tree_path <- grep(loci_name, all_trees, value = TRUE)
  loci_tree_path <- gsub("//","/", loci_tree_path) # Replace any double slashes (these occur if the folder name already has a slash at the end)
  if (file.exists(loci_tree_path)){
    return(loci_tree_path)
  } else {
    return("COULD NOT FIND PATH")
  }
}

