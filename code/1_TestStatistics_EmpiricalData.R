### empirical_treelikeness/code/1_TestStatistics_EmpiricalData.R
## R program to apply treelikeness statistics to transcriptomes from empirical data
## Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
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
print("initialising namespace")
# input_dir         <- the folder(s) containing the empirical data
# input_names       <- set name(s) for the dataset(s)
# output_dir        <- for collated output and results. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# cores_to_use      <- the number of cores to use for parametric bootstrap. 1 for a single core (wholly sequential), or higher if using parallelisation.
# cores_for_iqtree  <- specify the number of threads for IQ-Tree to use. If using parallelisation for parametric bootstrap, use 1. Otherwise, can use 1 or set to "AUTO" and let IQ-Tree select the best number
# reps_to_do        <- the number of of bootstrap replicates to perform
# exec_folder       <- the folder containing the software executables needed for analysis (3SEQ, IQ-Tree and SplitsTree4)
# exec_paths        <- location to each executable within the folder. Attach the names of the executables so the paths can be accessed by name

# The SplitsTree executable path can be tricky to find: 
#       - in MacOS, the path is "SplitsTree.app/Contents/MacOS/JavaApplicationStub" (assuming you are in the same directory as the application)
#       - in Linux, after installing and navigating into the folder it's simply "SplitsTree"

# # UNCOMMENT THE FOLLOWING LINES AND ENTER YOUR FILE PATHS/VARIABLES
# input_dir <- ""
# input_names <- ""
# output_dir <- ""
# treedir <- ""
# maindir <- ""
# cores_to_use <- 1
# cores_for_iqtree <- 1
# reps_to_do <- 199
# sCF_replicates <- 1000
# # Create a vector with all of the executable file paths
# # To access a path: exec_paths[["name"]]
# exec_folder <- "/path/to/executables/folder/"
# exec_paths <- c("3seq_executable","IQ-Tree_executable","SplitsTree_executable")
# exec_paths <- paste0(exec_folder,exec_paths)
# names(exec_paths) <- c("3seq","IQTree","SplitsTree")

### Caitlin's paths ###
# run_location = "local"
run_location = "server"

if (run_location == "local"){
  input_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                 "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Misof2014/loci/",
                 "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/")
  best_model_paths <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/OKP_loci_bestmodel.txt",
                        "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Misof2014/Misof2014_orthologousGenes_bestmodel.txt",
                        NA)
  input_names <- c("1KP", "Misof2014","Vanderpool2020")
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/")
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
  cores_for_iqtree = 1
  reps_to_do = 6
  sCF_replicates = 1000
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



##### Step 4: Extract names and locations of loci #####
# Obtaining the list of file paths from 1KP is the messiest as each alignment in it a separate folder, where the folder's name is the gene number
# Then, extract the best model for each loci (to feed into IQ-Tree - because we want to use as many of the original paramaters as we can!)
OKP_paths <- paste0(input_dir[["1KP"]], list.files(input_dir[["1KP"]], recursive = TRUE, full.names = FALSE))
OKP_names <- list.files(input_dir[["1KP"]], full.names = FALSE)
if (is.na(best_model_paths[["1KP"]])){
  OKP_model <- rep("MFP", length(OKP_paths))
} else {
  OKP_model <- model.from.partition.scheme(OKP_names,best_model_paths[["1KP"]],"1KP")
}
OKP_allowed_missing_sites <- 0.5 # Remove any sequence that has less than half the sites present
# Obtaining the list of loci file paths from Misof 2014 is easy -- all the loci are in the same folder
Misof2014_paths <- paste0(input_dir[["Misof2014"]], list.files(input_dir[["Misof2014"]], full.names = FALSE))
Misof2014_names <- gsub(".nex","",grep(".nex",unlist(strsplit((Misof2014_paths), "/")), value = TRUE))
if (is.na(best_model_paths[["Misof2014"]])){
  Misof2014_model <- rep("MFP", length(Misof2014_paths))
} else {
  Misof2014_model <- model.from.partition.scheme(Misof2014_names,best_model_paths[["Misof2014"]],"Misof2014")
}
Misof2014_allowed_missing_sites <- 0.5 # Remove any sequence that has less than half the sites present
# Obtaining the list of loci file paths from Vanderpool 2020 is easy -- all the loci are in the same folder
Vanderpool2020_paths <- paste0(input_dir[["Vanderpool2020"]], list.files(input_dir[["Vanderpool2020"]], full.names = FALSE))
Vanderpool2020_names <- gsub("_NoNcol.Noambig.fa","",grep(".fa",unlist(strsplit((Vanderpool2020_paths), "/")), value = TRUE))
if (is.na(best_model_paths[["Vanderpool2020"]])){
  # no set of models for these loci. 
  # Use "-m MFP" in IQ-Tree to automatically set best model - see Vanderpool et al (2020)
  Vanderpool2020_model <- rep("MFP", length(Vanderpool2020_paths))
}
Vanderpool2020_allowed_missing_sites <- NA # Allow the Vanderpool dataset to run as it it
# Create a dataframe of loci information for all three datasets: loci name, alphabet type, model, dataset, path, output path
loci_df <- data.frame(loci_name = c(Vanderpool2020_names, Misof2014_names, OKP_names),
                      alphabet = c(rep("dna", length(Vanderpool2020_paths)), rep("protein", length(Misof2014_paths)), rep("protein",length(OKP_paths))),
                      best_model = c(Vanderpool2020_model, Misof2014_model, OKP_model),
                      dataset = c(rep("Vanderpool2020", length(Vanderpool2020_paths)), rep("Misof2014",length(Misof2014_paths)), rep("1KP",length(OKP_paths))),
                      loci_path = c(Vanderpool2020_paths, Misof2014_paths, OKP_paths),
                      output_folder = c(rep(output_dirs[["Vanderpool2020"]], length(Vanderpool2020_paths)), rep(output_dirs[["Misof2014"]], length(Misof2014_paths)), rep(output_dirs[["1KP"]],length(OKP_paths))),
                      allowable_proportion_missing_sites = c(rep(Vanderpool2020_allowed_missing_sites, length(Vanderpool2020_paths)),
                                                             rep(Misof2014_allowed_missing_sites, length(Misof2014_paths)),
                                                             rep(OKP_allowed_missing_sites,length(OKP_paths))),
                      stringsAsFactors = FALSE)

# Remove rows with best_model == NA <- these are the Misof2014 clans, protein domains and voids (whereas above, we extracted the models only for the orthologous genes)
# "!is.na()" means "is not NA" <- we want to keep only the rows where best_model is not NA
loci_df <- loci_df[!is.na(loci_df$best_model),]
# output loci_df <- save a record of the input parameters you used!
loci_df_name <- paste0(output_dir,"input_loci_parameters.csv")
write.csv(loci_df, file = loci_df_name)

# # Remove 1KP data for now (until I work out which species to include, or how to deal with the huge amount of species)
# loci_df <- data.frame(loci_name = c(Vanderpool2020_names, Misof2014_names), 
#                       alphabet = c(rep("dna", length(Vanderpool2020_paths)), rep("protein", length(Misof2014_paths))),
#                       best_model = c(Vanderpool2020_model, Misof2014_model),
#                       dataset = c(rep("Vanderpool2020", length(Vanderpool2020_paths)), rep("Misof2014",length(Misof2014_paths))),
#                       loci_path = c(Vanderpool2020_paths, Misof2014_paths), 
#                       output_folder = c(rep(output_dirs[["Vanderpool2020"]], length(Vanderpool2020_paths)), rep(output_dirs[["Misof2014"]], length(Misof2014_paths))),
#                       stringsAsFactors = FALSE)



##### Step 5: Calculate the test statistics and run the parametric bootstraps  #####
print("starting analysis")
print("apply treelikeness test statistics")
# To run locally for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths,
#                                                        number_of_replicates = 9, iqtree.num_threads = AUTO, iqtree.num_quartets = 100,
#                                                        num_of_cores = 1)
# Parameter choice explanation:
#       ~ iqtree.num_thread = 1: allows parallelisation higher up in the workflow 
#                              - i.e. program will run parametric bootstrap for test statistics with multiple threads
#                              - (number of simulatenous bootstraps set by choice of cores_to_use value)
#       ~ iqtree.num_quartets = 1000: greater than 100 quartets necessary for stable sCF values
# V_ids <- which(loci_df$dataset == "Vanderpool2020")
# lapply(V_ids, empirical.bootstraps.wrapper, loci_df, program_paths = exec_paths, number_of_replicates = reps_to_do, iqtree.num_threads = cores_for_iqtree,
#        iqtree.num_quartets = sCF_replicates, num_of_cores = cores_to_use)

M_ids <- which(loci_df$dataset == "Misof2014")
lapply(M_ids, empirical.bootstraps.wrapper, loci_df, program_paths = exec_paths, number_of_replicates = reps_to_do, iqtree.num_threads = cores_for_iqtree,
       iqtree.num_quartets = sCF_replicates, num_of_cores = cores_to_use)

O_ids <- which(loci_df$dataset == "1KP")
lapply(O_ids, empirical.bootstraps.wrapper, loci_df, program_paths = exec_paths, number_of_replicates = reps_to_do, iqtree.num_threads = cores_for_iqtree,
       iqtree.num_quartets = sCF_replicates, num_of_cores = cores_to_use)



##### Step 6: Collate test statistic results #####
print("collate results")
# Collate all the dataframes together
results_file <- paste0(output_dir,"empiricalAlignments_collated_testStatisticResults.csv")
results_df <- collate.bootstraps(directory = output_dir, file.name = "pValues", id = "", output.file.name = results_file)
 
 

# ##### Step 7: Extract additional information about the alignments and trees #####
# # Extracting more information for statistical tests and plots
# # Extract total tree depth from each alignment's iqtree file and add to dataframe
# results_df$total_tree_depth <- unlist(lapply(results_df$alignment_file, extract.total.tree.length))
# # Extract the number of parsimony informative sites from each alignment's iqtree file and add to dataframe
# results_df$num_parsimony_informative_sites <- unlist(lapply(results_df$alignment_file, extract.num.informative.sites))
# # Extract the GC content for each species in the alignment and add summary statistics to dataframe
# GC_vector <- unlist(lapply(results_df$alignment_file,calculate.GC.content))
# results_df$GC_content_mean <- GC_vector[c(TRUE,FALSE,FALSE)]
# results_df$GC_content_variance <- GC_vector[c(FALSE,TRUE,FALSE)]
# results_df$GC_content_sd <-  GC_vector[c(FALSE,FALSE,TRUE)]
# # Extract the gene trees (.treefile for each locus) from each alignment's IQ-Tree run and add them to the dataframe
# results_df$newick_tree <- unlist(lapply(results_df$alignment_file, extract.treefile.tree))
# 
# # Save the newly extended dataframe
# results_file <- paste0(output_dir,basename(input_dir),"_completeResults.csv")
# write.csv(results_df, file = results_file, row.names = FALSE)


