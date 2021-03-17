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



# This function takes a list of names of loci, randomly selects 100 pairs, compares those two trees and records the average difference between them
calculate.average.tree.distance <- function(species_tree_csv, loci_tree_folder, output_folder, sample_size = 100){
  # Extract the category name from the file name
  if (basename(species_tree_csv) == "all_loci_loci.csv"){
    category_name = "all_loci"
  } else {
    category_name <- gsub("p-value_categories_","",basename(species_tree_csv))
    category_name <- gsub("_loci.csv","", category_name)
    if (category_name == "both"){
      category_name = "both_significant"
    }
    if (category_name == "none"){
      category_name = "neither_significant"
    } 
  }
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
  avg_dists <- c(mean(as.numeric(distances_df$RF_distance)), mean(as.numeric(distances_df$normalized_RF_distance)), 
                 mean(as.numeric(distances_df$weighted_RF_distance)),mean(as.numeric(distances_df$normalized_weighted_RF_distance)), 
                 mean(as.numeric(distances_df$path_difference)),mean(as.numeric(distances_df$weighted_path_difference)), 
                 mean(as.numeric(distances_df$KF_distance)))
  avg_dists <- signif(avg_dists, digits = 4)
  avg_dists <- c(category_name, avg_dists)
  names(avg_dists) <- c("category","mean_RF_distance", "mean_normalized_RF_distance", 
                        "mean_weighted_RF_distance", "mean_normalized_weighted_RF_distance", 
                        "mean_path_difference", "mean_weighted_path_difference",
                        "mean_KF_distance")
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



# This function takes an alignment and calculates the AU test using IQ-Tree
perform.AU.test <- function(loci_name, data_folder, output_folder, csv_folder, three_trees_path, iqtree_path){
  copy_path <- paste0(output_folder, loci_name, ".fasta")
  alignment_path <- copy_path
  # If the log file doesn't exist, run IQ-Tree
  if (file.exists(paste0(alignment_path,".log")) == FALSE){
    # Find the loci file in the data folder
    all_data_files <- list.files(data_folder)
    loci_file <- paste0(data_folder, grep(loci_name, all_data_files, value = TRUE))
    # Copy the file into the output folder
    file.copy(from = loci_file, to = copy_path, overwrite = FALSE, copy.mode = TRUE)
    # # List all files and identify tree path
    # all_files <- list.files(loci_folder)
    # tree_path <- grep(".treefile", all_files, value = TRUE)
    # tree_path <- paste0(loci_folder, grep(".treefile.", tree_path, invert = TRUE, value = TRUE))
    # Conduct the analysis using the copied file
    # Construct IQ-tree command
    iqtree_command <- paste0(iqtree_path, " -s ", alignment_path, " -z ", three_trees_path, " -au")
    # run command
    system(iqtree_command)
  }
  # Collect the results of the AU test
  # Open the log file
  log_file <- paste0(alignment_path, ".log")
  log_lines <- readLines(log_file)
  ind <- grep("Reading trees in",log_lines) + 2
  # extract log likelihood values
  tree1_logl <- log_lines[ind]
  tree1_logl <- as.numeric(strsplit(tree1_logl, ":")[[1]][2])
  tree2_logl <- log_lines[ind + 1]
  tree2_logl <- as.numeric(strsplit(tree2_logl, ":")[[1]][2])
  tree3_logl <- log_lines[ind + 2]
  tree3_logl <- as.numeric(strsplit(tree3_logl, ":")[[1]][2])
  logl_sum <- tree1_logl + tree2_logl + tree3_logl
  tree1_prop_logl <- tree1_logl/logl_sum
  tree2_prop_logl <- tree2_logl/logl_sum
  tree3_prop_logl <- tree3_logl/logl_sum
  find_best_tree <- which(c(tree1_logl, tree2_logl, tree3_logl) == max(c(tree1_logl, tree2_logl, tree3_logl)))
  best_tree_number <- paste(c(find_best_tree), collapse = " + ")
  if (length(find_best_tree) == 1){
    best_tree_exists = "YES"
  } else if (length(find_best_tree) > 1){
    best_tree_exists = "NO"
  }
  if ((round(tree1_logl, 4) ==  round(tree2_logl, 4)) & (round(tree2_logl, 4) ==  round(tree3_logl, 4)) & (round(tree1_logl, 4) ==  round(tree3_logl, 4))){
    logl_equal = TRUE
  } else {
    logl_equal = FALSE
  }
  output_df <- data.frame(locus = loci_name, tree1_log_likelihood =  tree1_logl, tree2_log_likelihood = tree2_logl, 
                          tree3_log_likelihood = tree3_logl, sum_log_likelihood = logl_sum, best_tree = best_tree_number, 
                          one_tree_best = best_tree_exists, all_likelihoods_equal = logl_equal, tree1_likelihood_proportion = tree1_prop_logl, 
                          tree2_likelihood_proportion = tree2_prop_logl, tree3_likelihood_proportion = tree3_prop_logl)
  write.csv(output_df, file = paste0(csv_folder, loci_name, "_AU_test_results.csv"))
}




change.to.ternary.points <- function(ind, df){
  row <- df[ind,]
  point <- c(row$tree1_likelihood_proportion, row$tree2_likelihood_proportion, row$tree3_likelihood_proportion)
  return(point)
}

