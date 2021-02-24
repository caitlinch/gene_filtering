### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from loci with varying treelikeness
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics



##### Step 2: Set file paths and run variables #####
print("set filepaths")
# input_dir         <- the folder(s) containing the estimated trees from each loci
# alignment_dir     <- the folder(s) containing the alignments for each loci
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as input_dir and alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first in all three and same dataset last in all three)
# loci_to_remove    <- If there are any loci to exclude from the analysis, specify them here. 
#                   <- This parameter is a list containing a vector for each dataset from input_names. The vector names each loci to exclude from that dataset
# number_of_taxa    <- If there is a certain number of taxa required for your analysis, specify that here.
#                   <- This parameter is a list containing a vector for each dataset from input_names. The vector contains the allowable number of taxa for that dataset
# treelikeness_df   <- csv file containing collated treelikeness test statistics for each loci (created by running the code/1_TestStatistics_EmpiricalData.R file) 
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# cores_to_use      <- the number of cores to use for parametric bootstrap. 1 for a single core (wholly sequential), or higher if using parallelisation.
# exec_folder       <- the folder containing the software executables needed for analysis (ASTRAL and IQTREE)
# exec_paths        <- location to each executable within the folder
# datasets_to_copy_loci   <- Out of the input names, select which datasets to copy loci trees for tree estimation
# datasets_to_estimate_trees <- Out of the input names, select which datasets to estimate species trees based on treelikeness results
# partition.by.codon.position <- Whether to run analysis partitioing by codon position
#                             <- set TRUE if you want to estimate species trees partitioning by codon position, and FALSE if you don't
# loci_windows     <- select how many loci to include in each species tree estimate when ranking loci by treelikeness

# The SplitsTree executable path can be tricky to find: 
#       - in MacOS, the path is "SplitsTree.app/Contents/MacOS/JavaApplicationStub" (assuming you are in the same directory as the application)
#       - in Linux, after installing and navigating into the folder it's simply "SplitsTree"

# # UNCOMMENT THE FOLLOWING LINES AND ENTER YOUR FILE PATHS/VARIABLES
# input_dir <- ""
# input_names <- ""
# treelikeness_df <- ""
# output_dir <- ""
# treedir <- ""
# maindir <- ""
# cores_to_use <- 1
# # Create a vector with all of the executable file paths using the following lines as a template:
# # exec_folder <- "/path/to/executables/folder/"
# # exec_paths <- c("ASTRAL_executable","IQ-Tree_executable")
# # exec_paths <- paste0(exec_folder,exec_paths)
# # names(exec_paths) <- c("ASTRAL","IQTree")
# # To access a path: exec_paths[["name"]]
# exec_folder <- ""
# exec_paths <- c()
# exec_paths <- paste0(exec_folder, exec_paths)
# datasets_to_copy_loci <-  c()
# datasets_to_estimate_trees <- c()
# partition.by.codon.position = TRUE
# loci_windows     <- c(10, 50, 100, 250, 500)

### Caitlin's paths ###
# run_location = "local"
run_location = "server"

if (run_location == "local"){
  # Datasets/dataset information
  input_names <- "Vanderpool2020"
  loci_to_remove <- list("Vanderpool2020" = "ORTHOMCL14552")
  number_of_taxa <- list("Vanderpool2020" = 29)
  
  # File and directory locations
  input_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/Vanderpool2020_trees")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/")
  treelikeness_df_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/empiricalTreelikeness_Vanderpool2020_collated_results_20210120.csv"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths  in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("Astral/astral.5.7.5.jar","iqtree-2.0-rc1-MacOSX/bin/iqtree")
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # set number of cores for parallel processing
  cores_to_use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("Vanderpool2020")
  datasets_to_estimate_trees <- c("Vanderpool2020")
  partition.by.codon.position = TRUE
  
  # Select how many loci to include in each species tree estimate
  loci_windows     <- c(10, 50, 100, 250, 500)
  
} else if (run_location=="server"){
  # Datasets/dataset information
  input_names <- "Vanderpool2020"
  loci_to_remove <- list("Vanderpool2020" = "ORTHOMCL14552")
  number_of_taxa <- list("Vanderpool2020" = 29)
  
  # File and directory locations
  input_dir <- "/data/caitlin/empirical_treelikeness/Output/Vanderpool2020_trees/"
  alignment_dir <- "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/"
  treelikeness_df_file <- "/data/caitlin/empirical_treelikeness/Output/empiricalTreelikeness_Vanderpool2020_collated_results_20210120.csv"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness repository/folder is
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness repository/folder is 
  
  # Create a vector with all of the executable file paths in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/executables/ASTRAL/astral.5.7.5.jar","/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # set number of cores  for parallel processing
  cores_to_use = 15
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("Vanderpool2020")
  datasets_to_estimate_trees <- c("Vanderpool2020")
  partition.by.codon.position = TRUE
  
  # Select how many loci to include in each species tree estimate
  loci_windows     <- c(10, 50, 100, 250, 500)
}



##### Step 3: Source files for functions #####
# Source the functions using the filepaths from Step 2
source(paste0(maindir,"code/func_empirical.R"))
source(paste0(maindir,"code/func_analysis.R"))



##### Step 4: Set up dataframes for analyses #####
# Add dataset names the input_dir and alignment_dir variables so you can index by dataset name
names(input_dir) <- input_names
names(alignment_dir) <- input_names
# Create output folders for each data set if they don't exist
output_dirs <- paste0(output_dir,input_names,"/")
names(output_dirs) <- input_names
for (d in output_dirs){
  if (file.exists(d) == FALSE){
    dir.create(d)
  }
}

# Identify any warnings from the IQ-Tree loci tree estimation
# Use these warnings to select which loci to exclude
for (dataset in datasets){
  # Open this dataset's raw output file from the treelikeness analysis 
  all_csv_files <- grep(".csv",list.files(csv_data_dir), value = TRUE)
  all_untrimmed_csv_files <- grep("trimmed",all_csv_files, value = TRUE, invert = TRUE)
  dataset_csv_file <- grep(dataset, all_untrimmed_csv_files, value = TRUE)
  dataset_df <- read.csv(paste0(csv_data_dir,dataset_csv_file), stringsAsFactors = FALSE)
  
  # Take list of alignments from the raw output file
  all_alignments <- dataset_df$alignment_file
  
  # Collect warnings and write out as a csv file
  warning_df <- as.data.frame(do.call(rbind, (lapply(all_alignments, check.for.IQTree.warnings))))
  warning_df_file <- paste0(csv_data_dir, dataset, "_collated_IQ-Tree_warnings.csv")
  write.csv(warning_df, file = warning_df_file)
}

# Open the treelikeness results dataframe
treelikeness_df <- read.csv(treelikeness_df_file, stringsAsFactors = FALSE)
# Remove loci to remove
rm_inds <- c()
for (dataset in input_names){
  rm_loci_ds <- loci_to_remove[[dataset]]
  rm_inds_ds <- which(treelikeness_df$dataset == dataset & treelikeness_df$loci %in% rm_loci_ds)
  rm_inds <- c(rm_inds, rm_inds_ds)
}
keep_rows <- setdiff(c(1:nrow(treelikeness_df)),rm_inds)
treelikeness_df <- treelikeness_df[keep_rows,]
# Remove loci with the wrong number of taxa
keep_inds <- c()
for (dataset in input_names){
  keep_taxa_num_ds <- number_of_taxa[[dataset]]
  keep_inds_ds <- which(treelikeness_df$dataset == dataset & treelikeness_df$n_taxa %in% keep_taxa_num_ds)
  keep_inds <- c(keep_inds, keep_inds_ds)
}
treelikeness_df <- treelikeness_df[keep_inds,]
# Save trimmed treelikeness_df
trimmed_treelikeness_df_file <- gsub(".csv", "_trimmedLoci_trimmedNTaxa.csv", treelikeness_df_file)
write.csv(treelikeness_df, file = trimmed_treelikeness_df_file)



##### Step 5: Categorise loci by treelikeness test p-values (3seq and tree proportion) #####
# Iterate through each of the datasets
# Sort loci into four categories based on treelikeness p-values
#     * both 3seq and tree proportion significant
#     * neither 3seq nor tree proportion significant
#     * only 3seq significant
#     * only tree proportion significant
# Save the trees from a category into a separate folder and collect all trees from a category into a collated text file
for (dataset in datasets_to_copy_loci){
  # filter treelikeness_df by dataset
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  # split loci into four groups (neither, 3seq, tp or both), then copy all loci alignments from each group into a new folder and all trees from each group into a new collated text file
  # 3seq p-value and tree proportion p-value both >0.05 (not significant)
  cat_none <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
  copy.loci.trees(cat_none, dataset_df[dataset_df$loci %in% cat_none,]$tree, output_dirs[dataset], "p-value_categories_none_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_none, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_none_IQ-Tree/"), mc.cores = cores_to_use)
  # 3seq p-value and tree proportion p-value both <=0.05 (significant)
  cat_both <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_both, dataset_df[dataset_df$loci %in% cat_both,]$tree, output_dirs[dataset], "p-value_categories_both_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_both, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_both_IQ-Tree/"), mc.cores = cores_to_use)
  # Only 3seq p-value <=0.05 and significant, tree proportion p-value not significant
  cat_3seq <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
  copy.loci.trees(cat_3seq, dataset_df[dataset_df$loci %in% cat_3seq,]$tree, output_dirs[dataset], "p-value_categories_3seq_only_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_3seq, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_3seq_only_IQ-Tree/"), mc.cores = cores_to_use)
  # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant
  cat_tp   <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_tp, dataset_df[dataset_df$loci %in% cat_tp,]$tree, output_dirs[dataset], "p-value_categories_tree_proportion_only_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_tp, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_tree_proportion_only_IQ-Tree/"), mc.cores = cores_to_use)

  # Repeat for a random sample of 50 trees (because an uneven amount of loci is in each group)
  # split loci into four groups (neither, 3seq, tp or both), then copy all loci alignments from each group into a new folder and all trees from each group into a new collated text file
  # 3seq p-value and tree proportion p-value both >0.05 (not significant)
  cat_none_50 <- sample(cat_none, 50)
  copy.loci.trees(cat_none_50, dataset_df[dataset_df$loci %in% cat_none_50,]$tree, output_dirs[dataset], "p-value_categories_none_50loci_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_none_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_none_50loci_IQ-Tree/"), mc.cores = cores_to_use)
  # 3seq p-value and tree proportion p-value both <=0.05 (significant)
  cat_both_50 <- sample(cat_both, 50)
  copy.loci.trees(cat_both_50, dataset_df[dataset_df$loci %in% cat_both_50,]$tree, output_dirs[dataset], "p-value_categories_both_50loci_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_both_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_both_50loci_IQ-Tree/"), mc.cores = cores_to_use)
  # Only 3seq p-value <=0.05 and significant, tree proportion p-value not significant
  cat_3seq_50 <- sample(cat_3seq, 50)
  copy.loci.trees(cat_3seq_50, dataset_df[dataset_df$loci %in% cat_3seq_50,]$tree, output_dirs[dataset], "p-value_categories_3seq_only_50loci_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_3seq_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_3seq_only_50loci_IQ-Tree/"), mc.cores = cores_to_use)
  # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant
  cat_tp_50 <- sample(cat_tp, 50)
  copy.loci.trees(cat_tp_50, dataset_df[dataset_df$loci %in% cat_tp_50,]$tree, output_dirs[dataset], "p-value_categories_tree_proportion_only_50loci_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(cat_tp_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_tree_proportion_only_50loci_IQ-Tree/"), mc.cores = cores_to_use)

  # If partitioning by codon position, repeat the IQ-Tree set-up
  if (partition.by.codon.position == TRUE){
    # filter treelikeness_df by dataset
    dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
    # split loci into four groups (neither, 3seq, tp or both), then copy all loci alignments from each group into a new folder and all trees from each group into a new collated text file
    # 3seq p-value and tree proportion p-value both >0.05 (not significant)
    cat_none <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
    mclapply(cat_none, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_none_IQ-Tree_partition/"), mc.cores = cores_to_use)
    cat_none_50 <- sample(cat_none, 50)
    mclapply(cat_none_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_none_50loci_IQ-Tree_partition/"), mc.cores = cores_to_use)
    # 3seq p-value and tree proportion p-value both <=0.05 (significant)
    cat_both <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
    mclapply(cat_both, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_both_IQ-Tree_partition/"), mc.cores = cores_to_use)
    cat_both_50 <- sample(cat_both, 50)
    mclapply(cat_both_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_both_50loci_IQ-Tree_partition/"), mc.cores = cores_to_use)
    # Only 3seq p-value <=0.05 and significant, tree proportion p-value not significant
    cat_3seq <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
    mclapply(cat_3seq, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_3seq_only_IQ-Tree_partition/"), mc.cores = cores_to_use)
    cat_3seq_50 <- sample(cat_3seq, 50)
    mclapply(cat_3seq_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_3seq_only_50loci_IQ-Tree_partition/"), mc.cores = cores_to_use)
    # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant
    cat_tp   <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
    mclapply(cat_tp, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_tree_proportion_only_IQ-Tree_partition/"), mc.cores = cores_to_use)
    cat_tp_50 <- sample(cat_tp, 50)
    mclapply(cat_tp_50, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_tree_proportion_only_50loci_IQ-Tree_partition/"), mc.cores = cores_to_use)
  }
}

# Estimate a species tree for each of the four categories
for (dataset in datasets_to_estimate_trees){
  p_value_cat_files <- c("p-value_categories_none","p-value_categories_both","p-value_categories_3seq_only","p-value_categories_tree_proportion_only",
                         "p-value_categories_none_50loci","p-value_categories_both_50loci","p-value_categories_3seq_only_50loci","p-value_categories_tree_proportion_only_50loci")
  astral_inputs <- paste0(output_dirs[dataset], p_value_cat_files, "_ASTRAL.txt")
  iqtree_inputs <- paste0(output_dirs[dataset], p_value_cat_files,"_IQ-Tree")
  # Calculate the species tree using ASTRAL for each of the four categories
  # estimate.ASTRAL.species.tree(astral_inputs[1], gsub(".txt","_ASTRAL_species.tre",astral_inputs[1]), gsub(".txt","_ASTRAL_species.log",astral_inputs[1]), exec_paths["ASTRAL"]
  lapply(astral_inputs, ASTRAL.wrapper, exec_paths["ASTRAL"])
  # Calculate the species tree using IQ-Tree for each of the four categories
  lapply(iqtree_inputs, estimate.IQTREE.species.tree, exec_paths["IQTree"])

  # If partitioning by codon position, create a partition file for each folder and then estimate the tree using IQ-Tree
  if (partition.by.codon.position == TRUE){
    # Make list of folders for partition analysis
    partition_inputs <- paste0(iqtree_inputs,"_partition/")
    # Create the partition files
    mclapply(partition_inputs, make.partition.file, mc.cores = cores_to_use)
    # Run the analysis
    partition_files <- paste0(partition_inputs, "partitions.nex")
    lapply(partition_files, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"])
  }
}



##### Step 6: Categorise loci by treelikeness (tree proportion test statistic value) #####
# Save names of folders/files that were copied
astral_trees_to_save <- c()
iqtrees_to_save <- c()
partition_trees_to_save <- c()

# Set up alignments/trees to run species trees analysis
for (dataset in datasets_to_copy_loci){
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  # Rank loci by treelikeness using tree proportion test statistic values from low to high
  dataset_df <- dataset_df[order(dataset_df$tree_proportion),]
  # For each of the numbers in loci_windows, save that number of loci of the highest and of the lowest tree proportion values
  for (n in loci_windows){
    # Create name for folder/text file
    tl_name <- paste0("windows_treelike_",sprintf("%04d",n))
    not_tl_name <- paste0("windows_non-treelike_",sprintf("%04d",n))
    
    # Select best 'n' loci (last 'n' in list) and worst 'n' loci (first 'n' in list)
    tl_loci <- tail(treelikeness_df$loci, n)
    not_tl_loci <- treelikeness_df$loci[1:n]
    
    # Save best 'n' loci using the following functions:
    # copy.loci.trees(loci_names,loci_trees, output_folder, output_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
    # copy.loci.alignment(loci_name, dataset_loci_folder, new_alignment_location)
    copy.loci.trees(tl_loci, dataset_df[dataset_df$loci %in% tl_loci,]$tree, output_dirs[dataset], paste0(tl_name,"_ASTRAL"), copy.all.individually = FALSE, copy.and.collate = TRUE)
    mclapply(tl_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], tl_name, "_IQ-Tree/"), mc.cores = cores_to_use)
    
    # Save worst 'n' loci
    copy.loci.trees(not_tl_loci, dataset_df[dataset_df$loci %in% not_tl_loci,]$tree, output_dirs[dataset], paste0(not_tl_name,"_ASTRAL"), copy.all.individually = FALSE, copy.and.collate = TRUE)
    mclapply(not_tl_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], not_tl_name, "_IQ-Tree/"), mc.cores = cores_to_use)

    # Save names 
    astral_trees_to_save <- c(astral_trees_to_save, paste0(output_dirs[dataset], tl_name, "_ASTRAL.txt"), paste0(output_dirs[dataset], not_tl_name, "_ASTRAL.txt"))
    iqtrees_to_save <- c(iqtrees_to_save, paste0(output_dirs[dataset], tl_name, "_IQ-Tree/"), paste0(output_dirs[dataset], not_tl_name, "_IQ-Tree/"))
    
    if (partition.by.codon.position == TRUE){
      # Save best and worst 'n' trees
      mclapply(tl_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], tl_name, "_IQ-Tree_partition/"), mc.cores = cores_to_use)
      mclapply(not_tl_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], not_tl_name, "_IQ-Tree_partition/"), mc.cores = cores_to_use)
      # Save names
      partition_trees_to_save <- c(partition_trees_to_save , paste0(output_dirs[dataset], tl_name, "_IQ-Tree_partition/"), paste0(output_dirs[dataset], not_tl_name, "_IQ-Tree_partition/"))
    }
  }
}
# Save list of astral and iqtrees that were copied
write(astral_trees_to_save, file = paste0(output_dir,"ASTRAL_trees_saved.txt"))
write(iqtrees_to_save, file = paste0(output_dir,"iqtrees_saved.txt"))
if (partition.by.codon.position == TRUE){
  write(partition_trees_to_save, file = paste0(output_dir,"partition_trees_saved.txt"))
}

# Construct file/folder names for the datasets you want to estimate trees for 
astral_trees_to_estimate <- c()
iqtrees_to_estimate <- c()
partition_trees_to_estimate <- c()

for (dataset in datasets_to_estimate_trees){
  dataset_files <- list.files(output_dirs[[dataset]])
  window_files <- grep("windows_", dataset_files, value = TRUE)
  astral_trees_to_estimate <- c(astral_trees_to_estimate, paste0(output_dirs[[dataset]], grep("ASTRAL", window_files, value = TRUE)))
  window_iqtrees <- grep("IQ-Tree", window_files, value = TRUE)
  iqtrees_to_estimate <- c(iqtrees_to_estimate, paste0(output_dirs[[dataset]], grep("partition", window_iqtrees, invert = TRUE, value = TRUE),"/"))
  if (partition.by.codon.position == TRUE){
    partition_trees_to_estimate <- c(partition_trees_to_estimate, paste0(output_dirs[[dataset]], grep("partition", window_iqtrees, value = TRUE),"/")) 
  }
}
# Remove ASTRAL files from list
astral_trees_to_estimate <- grep("ASTRAL", astral_trees_to_estimate, value = TRUE)
astral_trees_to_estimate <- grep("_species.tre", astral_trees_to_estimate, invert = TRUE, value = TRUE)
astral_trees_to_estimate <- grep("_species.log", astral_trees_to_estimate, invert = TRUE, value = TRUE)
# Remove IQ-Tree files from list
iqtrees_to_estimate <- grep("IQ-Tree", iqtrees_to_estimate, value = TRUE) # Make sure these all contain "IQ-Tree"
iqtrees_to_estimate <- grep("_IQ-Tree\\.", iqtrees_to_estimate, invert = TRUE, value = TRUE) # Remove anything that's not a directory (only want the drectories)
# Save list of astral and iqtrees to estimate
write(astral_trees_to_estimate, file = paste0(output_dir,"ASTRAL_trees_to_estimate.txt"))
write(iqtrees_to_estimate, file = paste0(output_dir,"iqtrees_to_estimate.txt"))
if (partition.by.codon.position == TRUE){
  # Remove IQ-Tree files from list
  partition_trees_to_estimate <- grep("IQ-Tree", partition_trees_to_estimate, value = TRUE) # Make sure these all contain "IQ-Tree"
  partition_trees_to_estimate <- grep("_IQ-Tree_partition\\.", partition_trees_to_estimate, invert = TRUE, value = TRUE)
  write(partition_trees_to_estimate, file = paste0(output_dir,"partition_trees_to_estimate.txt")) 
}

# Run species trees analyses
# Run ASTRAL
lapply(astral_trees_to_estimate, ASTRAL.wrapper, exec_paths["ASTRAL"])
# Run IQ-Tree without partition model
lapply(iqtrees_to_estimate, estimate.IQTREE.species.tree, exec_paths["IQTree"])
if (partition.by.codon.position == TRUE){
  # Create partition files for IQ-Tree runs with partition model
  mclapply(partition_trees_to_estimate, make.partition.file, mc.cores = cores_to_use)
  partition_files <- paste0(partition_trees_to_estimate, "partitions.nex")
  lapply(partition_files, estimate.partitioned.IQTREE.species.tree, exec_paths[["IQTree"]]) 
}



##### Step 7: Species tree with all loci #####
for (dataset in datasets_to_copy_loci){
  # Subset dataframe to only this dataset
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  # Create name for ASTRAL and IQ-Tree runs
  astral_name <- "all_loci_ASTRAL"
  iqtree_name <- "all_loci_IQTree/"
  # Copy all loci and all trees
  all_loci <- dataset_df$loci
  all_trees <- dataset_df$tree
  copy.loci.trees(all_loci, all_trees, output_dirs[dataset], astral_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
  mclapply(all_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], iqtree_name), mc.cores = cores_to_use)
  
  if (partition.by.codon.position == TRUE){
    partition_name <- "all_loci_IQTree_partition/"
    mclapply(all_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], partition_name), mc.cores = cores_to_use)
  }
}

for (dataset in datasets_to_estimate_trees){
  # Calculate the species trees for all loci
  # estimate.ASTRAL.species.tree(astral_inputs[1], gsub(".txt","_ASTRAL_species.tre",astral_inputs[1]), gsub(".txt","_ASTRAL_species.log",astral_inputs[1]), exec_paths["ASTRAL"]
  ASTRAL.wrapper(paste0(output_dirs[[dataset]], astral_name,".txt"), exec_paths[["ASTRAL"]])
  # Calculate the species tree using IQ-Tree for each of the four categories
  estimate.IQTREE.species.tree(paste0(output_dirs[dataset],iqtree_name), exec_paths[["IQTree"]])
  
  if (partition.by.codon.position == TRUE){
    partition_name <- "all_loci_IQTree_partition/"
    make.partition.file(paste0(output_dirs[dataset], partition_name))
    estimate.partitioned.IQTREE.species.tree(paste0(output_dirs[[dataset]], partition_name, "partitions.nex"), exec_paths[["IQTree"]])
  }
}


