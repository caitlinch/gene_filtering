### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from loci with varying treelikeness
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
print("opening packages")
#library(ape) # analyses of phylogenetics and evolution
#library(parallel) # support for parallel computation
#library(phangorn) # phylogenetic reconstruction and analysis
#library(phytools) # tools for comparative biology and phylogenetics
#library(seqinr) # data analysis and visualisation for biological sequence data




##### Step 2: Set file paths and run variables #####
print("set filepaths")
# input_dir         <- the folder(s) containing the estimated trees from each loci
# alignment_dir     <- the folder(s) containing the alignments for each loci
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as input_dir and alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first in all three and same dataset last in all three)
# treelikeness_df   <- csv file containing collated treelikeness test statistics for each loci (created by running the code/1_TestStatistics_EmpiricalData.R file) 
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# cores_to_use      <- the number of cores to use for parametric bootstrap. 1 for a single core (wholly sequential), or higher if using parallelisation.
# exec_folder       <- the folder containing the software executables needed for analysis (ASTRAL and IQTREE)
# exec_paths        <- location to each executable within the folder
# datasets_to_copy_loci   <- Out of the input names, select which datasets to copy loci trees for tree estimation
# datasets_to_estimate_trees <- Out of the input names, select which datasets to estimate species trees based on treelikeness results

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
#   datasets_to_copy_loci <-  c()
# datasets_to_estimate_trees <- c()

### Caitlin's paths ###
run_location = "local"
# run_location = "server"

if (run_location == "local"){
  input_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/Vanderpool2020_trees")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/")
  input_names <- "Vanderpool2020"
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
  
  # set number of cores
  cores_to_use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("Vanderpool2020")
  datasets_to_estimate_trees <- c("Vanderpool2020")
  
} else if (run_location=="server"){
  input_dir <- c(NULL)
  input_names <- "Vanderpool2020"
  treelikeness_df_file <- NULL
  output_dir <- NULL
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness repository/folder is
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness repository/folder is 
  
  # Create a vector with all of the executable file paths in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c(NULL,"/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("Astral","IQTree")
  
  # set number of cores and reps for bootstraps
  cores_to_use = 25
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("Vanderpool2020")
  datasets_to_estimate_trees <- c("Vanderpool2020")
}



##### Step 3: Source files for functions #####
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

# Open the treelikeness results dataframe
treelikeness_df <- read.csv(treelikeness_df_file)

# Source the functions using the filepaths from Step 2
#source(paste0(treedir,"code/func_test_statistic.R"))
#source(paste0(treedir,"code/func_process_data.R"))
#source(paste0(treedir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_empirical.R"))

##### Step 4: Partition loci by treelikeness test p-values (3seq and tree proportion) #####
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
  copy.loci.trees(cat_none, dataset_df[dataset_df$loci %in% cat_none,]$tree, output_dirs[dataset], "p-value_categories_none", copy.all.individually = FALSE, copy.and.collate = TRUE)
  lapply(cat_none, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_none/"))
  # 3seq p-value and tree proportion p-value both <=0.05 (significant)
  cat_both <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_both, dataset_df[dataset_df$loci %in% cat_both,]$tree, output_dirs[dataset], "p-value_categories_both", copy.all.individually = FALSE, copy.and.collate = TRUE)
  lapply(cat_both, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_both/"))
  # Only 3seq p-value <=0.05 and significant, tree proportion p-value not significant
  cat_3seq <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
  copy.loci.trees(cat_3seq, dataset_df[dataset_df$loci %in% cat_3seq,]$tree, output_dirs[dataset], "p-value_categories_3seq_only", copy.all.individually = FALSE, copy.and.collate = TRUE)
  lapply(cat_3seq, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_3seq_only/"))
  # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant
  cat_tp   <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_tp, dataset_df[dataset_df$loci %in% cat_tp,]$tree, output_dirs[dataset], "p-value_categories_tree_proportion_only", copy.all.individually = FALSE, copy.and.collate = TRUE)
  lapply(cat_tp, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset],"p-value_categories_tree_proportion_only/"))
}

# Estimate a species tree for each of the four categories
for (dataset in datasets_to_estimate_trees){
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  p_value_cat_files <- grep("p-value_categories",list.files(output_dirs[dataset]), value = TRUE)
  astral_inputs <- paste0(output_dirs[dataset], grep(".txt", p_value_cat_files, value = TRUE))
  iqtree_inputs <- paste0(output_dirs[dataset], grep(".txt", p_value_cat_files, value = TRUE, invert = TRUE))
  # Calculate the species tree using ASTRAL for each of the four categories
  estimate.ASTRAL.species.tree(astral_inputs[1], gsub(".txt","_ASTRAL_species.tre",astral_inputs[1]), gsub(".txt","_ASTRAL_species.log",astral_inputs[1]))
  estimate.ASTRAL.species.tree(astral_inputs[2], gsub(".txt","_ASTRAL_species.tre",astral_inputs[2]), gsub(".txt","_ASTRAL_species.log",astral_inputs[2]))
  estimate.ASTRAL.species.tree(astral_inputs[3], gsub(".txt","_ASTRAL_species.tre",astral_inputs[3]), gsub(".txt","_ASTRAL_species.log",astral_inputs[3]))
  estimate.ASTRAL.species.tree(astral_inputs[4], gsub(".txt","_ASTRAL_species.tre",astral_inputs[4]), gsub(".txt","_ASTRAL_species.log",astral_inputs[4]))
}

##### Step 5: Partition loci by treelikeness (tree proportion test statistic value) #####

