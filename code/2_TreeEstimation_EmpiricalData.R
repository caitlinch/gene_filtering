### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from loci with varying treelikeness
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
##     - 3SEQ (Lam et al 2018) (http://mol.ax/software/3seq/)
##     - Splitstree (Huson and Bryant 2006) (http://www.splitstree.org/) (need SplitsTree 4)
# Caitlin Cherryh 2019



##### Step 1: Open packages #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
library(seqinr) # data analysis and visualisation for biological sequence data
library(stringr) # wrappers for string operations



##### Step 2: Set file paths and run variables #####
print("set filepaths")
# input_dir         <- the folder(s) containing the estimated trees from each loci
# input_names       <- set name(s) for the dataset(s)
# treelikeness_df   <- csv file containing collated treelikeness test statistics for each loci (created by running the code/1_TestStatistics_EmpiricalData.R file) 
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# cores_to_use      <- the number of cores to use for parametric bootstrap. 1 for a single core (wholly sequential), or higher if using parallelisation.
# exec_folder       <- the folder containing the software executables needed for analysis (3SEQ, IQ-Tree and SplitsTree4)
# exec_paths        <- location to each executable within the folder. Attach the names of the executables so the paths can be accessed by name
# datasets_to_run   <- Out of the input names, select which datasets will have the treelikeness analysis run 

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
# # Create a vector with all of the executable file paths
# # To access a path: exec_paths[["name"]]
# exec_folder <- "/path/to/executables/folder/"
# exec_paths <- c("3seq_executable","IQ-Tree_executable","SplitsTree_executable")
# exec_paths <- paste0(exec_folder,exec_paths)
# names(exec_paths) <- c("3seq","IQTree","SplitsTree")
# datasets_to_run <- ""

### Caitlin's paths ###
# run_location = "local"
run_location = "server"

if (run_location == "local"){
  input_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/Vanderpool2020_trees")
  input_names <- "Vanderpool2020"
  treelikeness_df <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/empiricalTreelikeness_Vanderpool2020_collated_results_20210120.csv"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/SimulationsCodeAndResults/Executables/"
  # Create a vector with all of the executable file paths  in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree-2.0-rc1-MacOSX/bin/iqtree",
                  "SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
  
  # set number of cores and reps for bootstraps
  cores_to_use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_run <- c("Vanderpool2020")
  
} else if (run_location=="server"){
  input_dir <- c("/data/caitlin/empirical_treelikeness/Data_1KP/",
                 "/data/caitlin/empirical_treelikeness/Data_Misof2014/",
                 "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/")
  best_model_paths <- c("/data/caitlin/empirical_treelikeness/Data_inputFiles/OKP_loci_bestmodel.txt",
                        "/data/caitlin/empirical_treelikeness/Data_inputFiles/Misof2014_orthologousGenes_bestmodel.txt",
                        NA)
  input_names <- c("1KP", "Misof2014","Vanderpool2020")
  output_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness repository/folder is
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness repository/folder is 
  
  # Create a vector with all of the executable file paths in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
  
  # set number of cores and reps for bootstraps
  cores_to_use = 25
  cores_for_iqtree = 1
  reps_to_do= 99
  sCF_replicates = 1000
  
  # Select datasets to run analysis and collect results
  datasets_to_run <- c()
  datasets_to_collate <- c()
  datasets_to_collect_trees <- c("Vanderpool2020")
}



##### Step 3: Source files for functions #####
# Attach the names used throughout the code and functions to the exec_paths object
# Do not edit these names - it will cause internal functions and thus calculations to fail
names(exec_paths) <- c("3seq","IQTree","SplitsTree")
# Attach the input_names to the input_files 
names(input_dir) <- input_names
# Attach the names to the best model files
names(best_model_paths) <- input_names
# Create a set of output folders
output_dirs <- paste0(output_dir,input_names,"/")
names(output_dirs) <- input_names

# Source the functions using the filepaths from Step 2
source(paste0(treedir,"code/func_test_statistic.R"))
source(paste0(treedir,"code/func_process_data.R"))
source(paste0(treedir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_empirical.R"))