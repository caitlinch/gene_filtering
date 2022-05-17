### gene_filtering/code/1_TestStatistics_EmpiricalData.R
## R program to apply recombination detection methods to empirical data
# Caitlin Cherryh 2022

## Additional software packages are required:
##     - IQTREE2 (Minh et. al. 2020) (http://www.iqtree.org/)
##     - PHIPack (Bruen et. al. 2006) (https://www.maths.otago.ac.nz/~dbryant/software.html)
##     - GeneConv (Sawyer 1989) (https://www.math.wustl.edu/~sawyer/geneconv/)

## This script:
# 1. Generates a data frame with one row per locus, including all parameters necessary to estimate a gene tree for that locus
# 2. Applies the recombination detection methods to each locus: PHI and MaxChi (implemented in PHIPack; GeneConv (implemented in GeneConv)
# 3. Estimates a gene tree for each locus
# 4. Saves all gene trees to a separate folder
# 5. Identifies loci to exclude from further analysis (based on warnings from IQ-Tree .log and .iqtree files)



##### Step 1: Set file paths and run variables #####
# input_names       <- set name(s) for the dataset(s)
# input_dir         <- the folder(s) containing the empirical data
# best_model_paths  <- set path to file containing the best model of substitution for each loci. Set to NA to allow ModelFinder in IQ-Tree to choose best model.
# output_dir        <- for collated output and results. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# maindir           <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# cores_to_use      <- the number of cores to use for parallelisation. 1 for a single core (wholly sequential), or higher if using parallelisation.
# iqtree_num_threads<- specify number of cores for IQ-Tree to use during tree estimation. May be number, or set as "AUTO" for IQ-Tree to choose best number of threads
# exec_folder       <- the folder containing the software executables needed for analysis (PhiPack, GeneConv)
# exec_paths        <- location to each software executables needed for analysis (PhiPack, GeneConv). To attach names: names(exec_paths) <- c("PHIPack","GeneConv","IQTree"). To access a path: exec_paths[["name"]]

# Set which datasets you want to run through which analyses
# If do not want to run that part of the analysis, assign empty vector i.e. datasets_to_run <- c()
# If want to run specific datasets through that part of the analysis, assign only those. E.g. if you have datasets "Trees", "Animals" and "Fungi" and
#    want to run only "Trees" and "Fungi": datasets_to_run <- c("Trees", "Fungi")
# If want to run all of the datasets, assign all names i.e. datasets_to_run <- input_names

#   create_information_dataframe <- whether to gather information about each dataset required to run further analysis
#                              <- this information includes loci name, best model for each loci, location of alignment for each loci, etc
#                              <- if TRUE, program will collect all these variables in a dataframe and output a .csv containing this dataframe
#                              <- if FALSE, this step will be skipped
#   datasets_to_run   <- Out of the input names, select which datasets will have the treelikeness analysis run and the results collated. If running all, set datasets_to_run <- input_names
#   datasets_to_collect_trees <- Out of the input names, select which datasets will have the maximum likelihood trees from IQ-Tree collected and saved in a separate folder for easy downloading. 
#                             If saving all, set datasets_to_run <- input_names
#   datasets_to_check <- Out of the input names, select whuich datasets to collect the warnings from the IQ-Tree gene tree estimation log files

### Caitlin's paths ###
run_location = "local"

if (run_location == "local"){
  input_names <- c("1KP", "Whelan2017","Vanderpool2020", "Pease2016")
  input_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                 "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2017/genes/",
                 "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                 "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  best_model_paths <- c(NA, NA, NA, NA)
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/" # where the empirical treelikeness code is
  
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  # Create a vector with all of the executable file paths  in this order: PHIPack, GeneConv, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("Phi", "GENECONV_v1.81_unix.source/geneconv","iqtree-2.0-rc1-MacOSX/bin/iqtree")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("PHIPack","GeneConv","IQTree")
  
  # set number of cores for parallelisation
  cores_to_use = 1
  iqtree_num_threads = "AUTO"
} else if (run_location=="server"){
  input_names <- c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  input_dir <- c("/data/caitlin/empirical_treelikeness/Data_1KP/",
                 "/data/caitlin/empirical_treelikeness/Data_Whelan2017/",
                 "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/",
                 "/data/caitlin/empirical_treelikeness/Data_Pease2016/")
  best_model_paths <- c(NA, NA, NA, NA)
  output_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness repository/folder is 
  
  # Create a vector with all of the executable file paths in this order: PHIPack, GeneConv, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/executables/GENECONV_v1.81_unix.source/geneconv", 
                  "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("PHIPack","GeneConv","IQTree")
  # Extend the executable paths
  exec_paths <- paste0(exec_folder, exec_paths)
  
  # set number of cores for parallelisation
  cores_to_use = 30
  iqtree_num_threads = "AUTO"
}

# Select which analyses to apply to each dataset
create_information_dataframe <- TRUE
datasets_to_run <- input_names
datasets_to_collect_trees <- input_names
datasets_to_check <- input_names
### End Caitlin's paths ###



##### Step 2: Source packages and files for functions #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(seqinr) # data analysis and visualisation for biological sequence data
print("sourcing functions")
source(paste0(maindir,"code/func_recombination_detection.R"))
source(paste0(maindir,"code/func_empirical.R"))
source(paste0(maindir,"code/func_analysis.R"))



##### Step 3: Extract names and locations of loci #####
# Attach the input_names to the input_files and model paths
names(input_dir) <- input_names
names(best_model_paths) <- input_names
# Create a set of output folders
output_dirs <- paste0(output_dir,input_names,"/")
names(output_dirs) <- input_names
for (d in output_dirs){
  if (file.exists(d) == FALSE){
    dir.create(d)
  }
}

if (create_information_dataframe == TRUE){
  ### One Thousand Plants dataset
  # Obtaining the list of file paths from 1KP is the messiest as each alignment in it a separate folder, where the folder's name is the gene number
  # Then, extract the best model for each loci (to feed into IQ-Tree - because we want to use as many of the original paramaters as we can!)
  OKP_paths <- paste0(input_dir[["1KP"]], list.files(input_dir[["1KP"]], recursive = TRUE, full.names = FALSE))
  OKP_names <- list.files(input_dir[["1KP"]], full.names = FALSE)
  if (is.na(best_model_paths[["1KP"]])){
    # Select the same models as used for the original 1000 Plants runs (in the original paper)
    OKP_model <- rep("MFP", length(OKP_paths))
  } else {
    # Allow ModelFinder in IQ-Tree to select the best model
    OKP_model <- model.from.partition.scheme(OKP_names,best_model_paths[["1KP"]],"1KP")
  }
  OKP_allowed_missing_sites <- NA # Don't remove any sequences based on the number of gaps/missing sites 
  
  ### Whelan 2015 dataset
  Whelan2017_paths <- paste0(input_dir[["Whelan2017"]], list.files(input_dir[["Whelan2017"]], full.names = FALSE))
  Whelan2017_names <- gsub("\\.fa","",basename(Whelan2017_paths))
  if (is.na(best_model_paths[["Whelan2017"]])){
    # no set of models for these loci.
    # Use "-m MFP" in IQ-Tree to automatically set best model
    Whelan2017_model <- rep("MFP", length(Whelan2017_paths))
  }
  Whelan2017_allowed_missing_sites <- NA # Don't remove any sequences based on the number of gaps/missing sites 
  
  ### Vanderpool 2020 dataset
  # Obtaining the list of loci file paths from Vanderpool 2020 is easy -- all the loci are in the same folder
  Vanderpool2020_paths <- paste0(input_dir[["Vanderpool2020"]], list.files(input_dir[["Vanderpool2020"]], full.names = FALSE))
  Vanderpool2020_names <- gsub("_NoNcol.Noambig.fa","",grep(".fa",unlist(strsplit((Vanderpool2020_paths), "/")), value = TRUE))
  if (is.na(best_model_paths[["Vanderpool2020"]])){
    # no set of models for these loci.
    # Use "-m MFP" in IQ-Tree to automatically set best model - see Vanderpool et al (2020)
    Vanderpool2020_model <- rep("MFP", length(Vanderpool2020_paths))
  }
  Vanderpool2020_allowed_missing_sites <- NA # Don't remove any sequences based on the number of gaps/missing sites 
  
  ### Pease 2016 dataset
  Pease2016_paths <- paste0(input_dir[["Pease2016"]], list.files(input_dir[["Pease2016"]], full.names = FALSE))
  Pease2016_names <- gsub("_", "-", gsub("cSL2.50", "", gsub("_100kb_windows.fasta", "", basename(Pease2016_paths))))
  if (is.na(best_model_paths[["Pease2016"]])){
    # Allow ModelFinder in IQ-Tree to select the best model
    Pease2016_model <- rep("MFP", length(Pease2016_paths))
  } else {
    # Select the same model used for the original Pease 2016 runs.
    # Open the info_paths[["Pease2016"]] file that contains extra information about the Pease 2016 loci
    Pease2016_info <- read.csv(file = info_paths[["Pease2016"]])
    Pease2016_model <- Pease2016_info$RAxML_model_input
    # The models taken from the RAxML info file are "GTRGAMMA", which is not a recognised input for model specification in IQ-Tree
    # The models are renamed "GTR+G", which is the terminology for the same model in IQ-Tree 
    Pease2016_model <- gsub("GTRGAMMA","GTR+G", Pease2016_model)
  }
  Pease2016_allowed_missing_sites <- NA # Don't remove any sequences based on the number of gaps/missing sites 
  
  ### Compile datasets into one dataframe
  # Create a dataframe of loci information for all three datasets: loci name, alphabet type, model, dataset, path, output path
  loci_df <- data.frame(loci_name = c(Vanderpool2020_names, Whelan2017_names, OKP_names, Pease2016_names),
                        alphabet = c(rep("dna", length(Vanderpool2020_paths)), rep("protein", length(Whelan2017_paths)), 
                                     rep("protein",length(OKP_paths)), rep("dna", length(Pease2016_paths))),
                        best_model = c(Vanderpool2020_model, Whelan2017_model, OKP_model, Pease2016_model),
                        dataset = c(rep("Vanderpool2020", length(Vanderpool2020_paths)), rep("Whelan2017",length(Whelan2017_paths)), 
                                    rep("1KP",length(OKP_paths)), rep("Pease2016", length(Pease2016_paths))),
                        loci_path = c(Vanderpool2020_paths, Whelan2017_paths, OKP_paths, Pease2016_paths),
                        output_folder = c(rep(output_dirs[["Vanderpool2020"]], length(Vanderpool2020_paths)), rep(output_dirs[["Whelan2017"]], length(Whelan2017_paths)), 
                                          rep(output_dirs[["1KP"]],length(OKP_paths)), rep(output_dirs[["Pease2016"]], length(Pease2016_paths))),
                        allowable_proportion_missing_sites = c(rep(Vanderpool2020_allowed_missing_sites, length(Vanderpool2020_paths)),
                                                               rep(Whelan2017_allowed_missing_sites, length(Whelan2017_paths)),
                                                               rep(OKP_allowed_missing_sites,length(OKP_paths)),
                                                               rep(Pease2016_allowed_missing_sites, length(Pease2016_paths))),
                        stringsAsFactors = FALSE)
  
  # output loci_df <- save a record of the input parameters you used!
  if (setequal(datasets_to_run, input_names) == TRUE){
    loci_df_name <- paste0(output_dir,"01_AllDatasets_RecombinationDetection_input_loci_parameters.csv")
    write.csv(loci_df, file = loci_df_name, row.names = FALSE)
  } else {
    loci_df_name <- paste0(output_dir,"01_",paste(sort(datasets_to_run), collapse = "_"),"_RecombinationDetection_input_loci_parameters.csv")
    write.csv(loci_df, file = loci_df_name, row.names = FALSE)
  }
}



##### Step 4: Apply the recombination detection methods  #####
if (length(datasets_to_run) > 0){
  # Rerun recombination.detection.wrapper to iterate through loci, run recombination detection methods,
  # extract all RecombinationDetection_results files, and save the output
  print("run recombination detection methods and collect output csv files")
  dataset_ids <- which(loci_df$dataset %in% datasets_to_run)
  run_list <- mclapply(dataset_ids, recombination.detection.wrapper, df = loci_df, executable_paths = exec_paths, iqtree_num_threads, mc.cores = cores_to_use)
  run_df <- as.data.frame(do.call(rbind, run_list))
  print("save results file as a .csv")
  if (setequal(datasets_to_run, input_names) == TRUE){
    results_file <- paste0(output_dir,"01_AllDatasets_RecombinationDetection_complete_collated_results.csv")
    write.csv(run_df, file = results_file, row.names = FALSE)
  } else {
    results_file <- paste0(output_dir,"01_",paste(sort(datasets_to_run), collapse = "_"),"_RecombinationDetection_complete_collated_results.csv")
    write.csv(run_df, file = results_file, row.names = FALSE)
  }
}



##### Step 5: Collate trees #####
print("collate trees")
# Check if there are any trees to collate
if (length(datasets_to_collect_trees) > 0){
  # Iterate through each dataset
  for (dataset in datasets_to_collect_trees){
    # Create output folder path for trees and create folder if it doesn't exist
    op_tree_folder <- paste0(output_dir,dataset,"_trees/")
    if (!dir.exists(op_tree_folder)){
      dir.create(op_tree_folder)
    }
    # Start by getting each loci folder
    all_ds_folder <- paste0(output_dirs[[dataset]], list.dirs(output_dirs[[dataset]], recursive = FALSE, full.names = FALSE))
    # Want to go through each loci folder and save tree into op_tree_folder
    # Any missing trees will be ignored
    lapply(all_ds_folder, save.tree, trees_folder = op_tree_folder)
  }
}



##### Step 6: Identify warnings from IQ-Tree runs #####
# Iterate through the warning files from estimating each gene tree and record all IQ-Tree warnings
if (length(datasets_to_check)>0){
  # Identify any warnings from the IQ-Tree loci tree estimation
  # Use these warnings to select which loci to exclude
  for (dataset in datasets_to_check){
    # Collect alignment folders
    all_folders <- list.dirs(output_dirs[[dataset]])
    all_al_folders <- setdiff(all_folders, output_dirs[[dataset]])
    # Remove double slash in folder names
    all_al_folders <- gsub("//", "/", all_al_folders)
    # Check each folder for warnings
    warning_df <- as.data.frame(do.call(rbind, (lapply(all_al_folders, check.folder.for.IQTree.warnings))))
    # Save the warnings as a dataframe
    warning_df_file <- paste0(output_dir,"01_", dataset, "_collated_IQ-Tree_warnings.csv")
    write.csv(warning_df, file = warning_df_file, row.names = FALSE)
  }
}



##### Step 7: Identify loci to exclude from species tree estimation #####
# Select the csv files from datasets run through datasets_to_check
all_csvs <- list.files(output_dir)
dataset_warning_csvs <- paste0(output_dir, grep("IQ-Tree_warnings.csv", all_csvs, value = TRUE))
# Initialise vectors to store loci for exclusion in
exclusion_loci_name <- c()
exclusion_loci_dataset <- c()
exclusion_warning <- c()
# Iterate through the datasets one at a time
for (dataset in datasets_to_check){
  # Open IQ-Tree warnings record for this dataset
  w_df_file <- grep(dataset, dataset_warning_csvs, value = TRUE)
  w_df <- read.csv(w_df_file, stringsAsFactors = FALSE)
  # Remove duplicate rows from the dataframe
  w_df <- w_df[!duplicated(w_df),]
  # Check for the different kinds of common warnings
  # Check for WARNING: Number of threads seems too high/low for the alignment given the length
  # Recommends using -nt AUTO to select the best number of threads
  nt_auto_ids <- grep("Number of threads seems too",w_df$warnings)
  # Check for WARNING: x near-zero internal branches - appears when some internal branches are very small
  near_zero_ids <- grep("near-zero internal branches", w_df$warnings)
  # Check for WARNING: bootstrap analysis did not converge (needed higher number of iterations using -nm option)
  bs_ids <- grep("bootstrap analysis did not converge", w_df$warnings)
  # Check for WARNING: Log-likelihood for model x worse than for model y
  # This warning appears during model selection when a model has a better AIC but worse likelihood than a previous model
  logl_val_ids <- grep("WARNING: Log-likelihood ", w_df$warnings) 
  # Check for WARNING: x sequences contain more than 50% gaps/ambiguity
  gaps_ids <- grep("sequences contain more than 50% gaps/ambiguity", w_df$warnings)
  # Check for long strings of "****************************" - these indicate start and end of warnings 
  star_ids <- grep("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", w_df$warnings)
  # Check for difficult datasets (by looking for situations where the number of NNIs does not converge)
  # This indicates that it's a difficult dataset with low phylogenetic information, but doesn't necessarily indicate low quality gene trees
  nni_inds <- grep("NNI search needs unusual large number of steps", w_df$warnings)
  # Check for different sequence names
  seq_name_inds <- grep("Some sequence names are changed as follows", w_df$warnings)
  # Check for the number of sites that contain only gaps or ambiguous characters
  only_gap_inds <- grep("contain only gaps or ambiguous characters", w_df$warnings)
  # If any other type of warning exists, exclude that loci from further analysis
  # The kinds of warnings that remain are:
  #     - Estimated model parameters are boundary that can cause numerical instability
  #     - NNI search need unusually large number of steps to converge
  #     - Some pairwise ML distances are too long (saturated)
  # Exceptions:
  #     - Ignore "some pairwise ML distances too long (saturated)" when the vast majority of loci alignments have this warning
  #         - Both the deep datasets (1KP and Whelan2017)
  if ((dataset == "1KP") | (dataset == "Whelan2017")){
    sat_inds <- grep("WARNING: Some pairwise ML distances are too long", w_df$warnings)
  } else {
    sat_inds <- c()
  }
  remaining_ids <- setdiff(1:nrow(w_df), sort(c(nt_auto_ids, near_zero_ids, bs_ids, logl_val_ids, 
                                                gaps_ids, star_ids, nni_inds, seq_name_inds, only_gap_inds,
                                                sat_inds)) )
  remaining_warnings_df <- w_df[remaining_ids,]
  # Add to vectors for outputting as a csv
  exclusion_loci_name <- c(exclusion_loci_name, remaining_warnings_df$loci)
  exclusion_loci_dataset <- c(exclusion_loci_dataset, remaining_warnings_df$dataset)
  exclusion_warning <- c(exclusion_warning, remaining_warnings_df$warnings)
}
exclusion_df <- data.frame(dataset = exclusion_loci_dataset,
                           loci = exclusion_loci_name,
                           warning = exclusion_warning)
if (setequal(datasets_to_check, input_names) == TRUE){
  exclusion_op_name <- paste0(output_dir, "01_AllDatasets_IQ-Tree_warnings_LociToExclude.csv")
  write.csv(exclusion_df, exclusion_op_name, row.names = FALSE)
} else {
  exclusion_op_name <- paste0(output_dir, "01_", paste(datasets_to_check, collapse = "_"), "_IQ-Tree_warnings_LociToExclude.csv")
  write.csv(exclusion_df, exclusion_op_name, row.names = FALSE) 
}


