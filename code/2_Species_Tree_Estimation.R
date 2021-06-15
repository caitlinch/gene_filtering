### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from loci with varying treelikeness
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
# Caitlin Cherryh 2021

##### Step 1: Set file paths and run variables #####
# alignment_dir     <- the folder(s) containing the alignments for each loci
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first in all three and same dataset last in all three)
# loci_to_remove    <- If there are any loci to exclude from the analysis, specify them here. 
#                   <- This parameter is a list containing a vector for each dataset from input_names. The vector names each loci to exclude from that dataset
# number_of_taxa    <- If there is a certain number of taxa required for your analysis, specify that here.
#                   <- This parameter is a list containing a vector for each dataset from input_names. The vector contains the allowable number of taxa for that dataset
# treelikeness_df   <- csv file containing collated treelikeness test statistics for each loci (created by running the code/1_TestStatistics_EmpiricalData.R file) 
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
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
# input_names <- ""
# treelikeness_df <- ""
# output_dir <- ""
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
run_location = "local"

if (run_location == "local"){
  # Datasets/dataset information
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  loci_to_remove <- list("Vanderpool2020" = "ORTHOMCL14552")
  number_of_taxa <- list("Vanderpool2020" = NA)
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  treelikeness_df_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/RecombinationDetection_Vanderpool2020_collated_results_complete.csv"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
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
  check.for.warnings = FALSE # check IQ-Tree .log file and .iqtree file output for each gene tree for warnings
  estimate.species.trees.in.IQTREE = TRUE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd - based on position in alignment file) 
  
} else if (run_location=="server"){
  # Datasets/dataset information
  input_names <- c("1KP", "Strassert2021","Vanderpool2020")
  loci_to_remove <- list("Vanderpool2020" = c("ORTHOMCL14552"))
  number_of_taxa <- list("Vanderpool2020" = NA)
  
  # File and directory locations
  alignment_dir <- "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/"
  csv_data_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  
  treelikeness_df_file <- "/data/caitlin/empirical_treelikeness/Output/empiricalTreelikeness_Vanderpool2020_collated_results_20210120.csv"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness repository/folder is
  
  # Create a vector with all of the executable file paths in this order: 3SEQ, IQ-Tree, SplitsTree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/executables/ASTRAL/astral.5.7.5.jar","/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # set number of cores  for parallel processing
  cores_to_use = 15
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("1KP", "Strassert2021","Vanderpool2020")
  datasets_to_estimate_trees <- c("1KP", "Strassert2021","Vanderpool2020")
  check.for.warnings = FALSE # check IQ-Tree .log file and .iqtree file output for each gene tree for warnings
  estimate.species.trees.in.IQTREE = FALSE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
}



##### Step 2: Open packages and source files for functions #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
# Source the functions using the filepaths from Step 2
source(paste0(maindir,"code/func_empirical.R"))
source(paste0(maindir,"code/func_analysis.R"))



##### Step 3: Set up dataframes for analyses #####
# Add dataset names to the alignment_dir variable so you can index by dataset name
names(alignment_dir) <- input_names
# Create output folders for each data set if they don't exist
output_dirs <- paste0(output_dir,input_names,"/")
names(output_dirs) <- input_names
for (d in output_dirs){
  if (file.exists(d) == FALSE){
    dir.create(d)
  }
}

# Iterate through the warning files from estimating each gene tree and record all IQ-Tree warnings
if (check.for.warnings == TRUE){
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
}

# Open the treelikeness results dataframe
treelikeness_df <- read.csv(treelikeness_df_file, stringsAsFactors = FALSE)
# Trim treelikeness df to remove loci with IQ-Tree warnings and loci with too few taxa
# Remove loci to remove
rm_inds <- c()
for (dataset in input_names){
  if (dataset %in% names(loci_to_remove)){ 
    if (!is.na(loci_to_remove[[dataset]])){
      rm_loci_ds <- loci_to_remove[[dataset]]
      rm_inds_ds <- which(treelikeness_df$dataset == dataset & treelikeness_df$loci %in% rm_loci_ds)
      rm_inds <- c(rm_inds, rm_inds_ds)
    }
  }
}
keep_rows <- setdiff(c(1:nrow(treelikeness_df)),rm_inds)
treelikeness_df <- treelikeness_df[keep_rows,]
# Remove loci with the wrong number of taxa
keep_inds <- c()
for (dataset in input_names){
  if (dataset %in% names(number_of_taxa)){
    if (is.na(number_of_taxa[[dataset]])){
      # If you want all the loci regardless of how many taxa they have, keep all of the loci in that dataset
      # e.g. if you want all the loci for Vanderpool 2020: number_of_taxa <- list("Vanderpool2020" = NA)
      keep_taxa_num_ds <- number_of_taxa[[dataset]]
      keep_inds_ds <- which(treelikeness_df$dataset == dataset)
      keep_inds <- c(keep_inds, keep_inds_ds)
    } else if (class(number_of_taxa[[dataset]]) == "numeric") {
      # If you only want to keep loci with certain numbers of loci, extract all the loci with that number of taxa
      # e.g. if you want all the loci with all 29 taxa for Vanderpool 2020: number_of_taxa <- list("Vanderpool2020" = 29)
      # e.g if only allowing up to 2 missing taxa for Vanderpool 2020: number_of_taxa <- list("Vanderpool2020" = c(27, 28, 29))
      keep_taxa_num_ds <- number_of_taxa[[dataset]]
      keep_inds_ds <- which(treelikeness_df$dataset == dataset & treelikeness_df$n_taxa %in% keep_taxa_num_ds)
      keep_inds <- c(keep_inds, keep_inds_ds)
    }
  } else {
    # If this dataset is not included in the number_of_taxa list, the default assumption is to keep all loci regardless of how many taxa are in that alignment
    # Keep all the loci
    keep_taxa_num_ds <- number_of_taxa[[dataset]]
    keep_inds_ds <- which(treelikeness_df$dataset == dataset)
    keep_inds <- c(keep_inds, keep_inds_ds)
  }
}
treelikeness_df <- treelikeness_df[keep_inds,]
# Create new columns with pass/fail for each test
# 3SEQ p-value
treelikeness_df$pass_3seq <- treelikeness_df$X3SEQ_p_value
treelikeness_df$pass_3seq[treelikeness_df$X3SEQ_p_value > 0.05] <- "TRUE"
treelikeness_df$pass_3seq[treelikeness_df$X3SEQ_p_value <= 0.05] <- "FALSE"
treelikeness_df$pass_3seq <- as.logical(treelikeness_df$pass_3seq)
# PHI test p-value
treelikeness_df$pass_phi<- treelikeness_df$PHI_normal_p_value
treelikeness_df$pass_phi[treelikeness_df$PHI_normal_p_value > 0.05] <- "TRUE"
treelikeness_df$pass_phi[treelikeness_df$PHI_normal_p_value <= 0.05] <- "FALSE"
treelikeness_df$pass_phi <- as.logical(treelikeness_df$pass_phi)
# MaxChi test p-value
treelikeness_df$pass_maxchi <- treelikeness_df$max_chi_squared_p_value
treelikeness_df$pass_maxchi[treelikeness_df$max_chi_squared_p_value > 0.05] <- "TRUE"
treelikeness_df$pass_maxchi[treelikeness_df$max_chi_squared_p_value <= 0.05] <- "FALSE"
treelikeness_df$pass_maxchi <- as.logical(treelikeness_df$pass_maxchi)
# GeneConv inner fragments simulated p-value
treelikeness_df$pass_geneconv_inner <- treelikeness_df$geneconv_inner_fragment_simulated_p_value
treelikeness_df$pass_geneconv_inner[treelikeness_df$geneconv_inner_fragment_simulated_p_value > 0.05] <- "TRUE"
treelikeness_df$pass_geneconv_inner[treelikeness_df$geneconv_inner_fragment_simulated_p_value <= 0.05] <- "FALSE"
treelikeness_df$pass_geneconv_inner <- as.logical(treelikeness_df$pass_geneconv_inner)
# GeneConv outer fragments simulated p-value
treelikeness_df$pass_geneconv_outer <- treelikeness_df$geneconv_outer_fragment_simulated_p_value
treelikeness_df$pass_geneconv_outer[treelikeness_df$geneconv_outer_fragment_simulated_p_value > 0.05] <- "TRUE"
treelikeness_df$pass_geneconv_outer[treelikeness_df$geneconv_outer_fragment_simulated_p_value <= 0.05] <- "FALSE"
treelikeness_df$pass_geneconv_outer <- as.logical(treelikeness_df$pass_geneconv_outer)
# GeneConv both p-value (TRUE if both inner and outer > 0.05, FALSE if one or both <= 0.05)
treelikeness_df$pass_geneconv <- "FALSE"
treelikeness_df$pass_geneconv[((treelikeness_df$geneconv_outer_fragment_simulated_p_value > 0.05) & 
                                 (treelikeness_df$geneconv_inner_fragment_simulated_p_value > 0.05))] <- "TRUE"
treelikeness_df$pass_geneconv <- as.logical(treelikeness_df$pass_geneconv)


# Save the trimmed treelikeness_df
trimmed_treelikeness_df_file <- gsub(".csv", paste0("_trimmedLoci_trimmedTaxa.csv"), treelikeness_df_file)
write.csv(treelikeness_df, file = trimmed_treelikeness_df_file)
# Save a df of just the pass/fail info
pass_df <- treelikeness_df[,c("dataset", "loci_name", "alphabet", "n_taxa", "n_bp", "pass_3seq", "pass_phi", "pass_maxchi", 
                              "pass_geneconv_inner", "pass_geneconv_outer", "pass_geneconv")]
pass_df_file <- gsub(".csv", paste0("_PassFail_record_.csv"), treelikeness_df_file)
write.csv(pass_df, file = pass_df_file)



##### Step 4: Categorize loci by test results and estimate species trees #####
# Iterate through each of the datasets
# For each recombination detection test, record which loci pass the test
# For each recombination detection test, estimate a species tree from the loci that pass the test
# Estimate a species tree from all loci
# Estimate a species tree from loci that pass every recombination detection test

### Save the loci trees (for ASTRAL) and the loci alignment (for IQ-Tree) 
for (dataset in datasets_to_copy_loci){
  # Create a row to store information about all the different variables
  summary_row <- c(dataset)
  # Create new folders to put these tree files/loci files and records in
  category_output_folder <- paste0(output_dirs[dataset], "species_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  text_records_dir <- paste0(output_dirs[dataset], "species_tree_records/")
  if (dir.exists(text_records_dir) == FALSE){
    dir.create(text_records_dir)
  }
  
  ## filter treelikeness_df by dataset
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  
  ## Iterate through each var and save the loci/trees for a tree made from only loci that pass the test
  # make a list of the variables on which to filter the loci - should be columns from the treelikeness_df
  vars <- c("X3SEQ_p_value", "PHI_normal_p_value", "max_chi_squared_p_value", "geneconv_inner_fragment_simulated_p_value")
  # Assign output names for each of the variables in vars
  vars_names <- c("3SEQ", "PHI", "maxchi", "geneconv")
  names(vars_names) <- vars
  # Iterate through each var: 
  for (v in vars){
    # Get short version of name for output files
    v_name <- vars_names[v]
    # Make names for output files
    v_text_name <- paste0(text_records_dir, dataset, "_", v_name, "_loci_record.txt")
    v_ASTRAL_name <- paste0(dataset,"_",v_name,"_ASTRAL")
    if (estimate.species.trees.in.IQTREE == TRUE){
      if (partition.by.codon.position == TRUE){
        v_IQTree_name <- paste0(dataset,"_",v_name,"_IQTREE_partitioned")
      } else if (partition.by.codon.position == FALSE){
        v_IQTree_name <- paste0(dataset,"_",v_name,"_IQTREE") 
      }
    }
    # Break up dataframe into only loci that pass the test
    if (v == "geneconv_inner_fragment_simulated_p_value"){
      # For GeneConv results, want to get rows that have a non-significant p-value for both the inner and outer fragment simulated p-values
      v_inds <- which((dataset_df$geneconv_inner_fragment_simulated_p_value > 0.05) & (dataset_df$geneconv_outer_fragment_simulated_p_value > 0.05))
      # For geneconv, also record number of NA values (number of loci that did not run in GeneConv)
      geneconv_NA_inds <- length(which((is.na(dataset_df$geneconv_inner_fragment_simulated_p_value) | 
                                          is.na(dataset_df$geneconv_outer_fragment_simulated_p_value > 0.05))))
      summary_row <- c(summary_row, geneconv_NA_inds)
    } else {
      # For any other variable, get all rows with a non-significant p-value
      v_inds <- which(dataset_df[,c(v)] > 0.05)
    }
    # Use the indexes to subset the dataframe to just loci that pass the test 
    # (i.e. have a non significant p-value, meaning the null hypothesis of treelikeness cannot be rejected)
    v_df <- dataset_df[v_inds,]
    # Copy trees of all loci that pass the test into one file that can be fed into ASTRAL
    copy.loci.trees(v_df$loci_name, v_df$tree, category_output_folder, v_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
    # If running IQ-Tree analysis, copy all loci into a separate folder that can be fed into IQ-Tree
    if (estimate.species.trees.in.IQTREE == TRUE){
      # create the partition file required to run this IQ-Tree analysis
      # partition.file.from.loci.list(loci_list, directory, original_alignment_folder, add.charpartition = FALSE)
      partition.file.from.loci.list(loci_list = v_df$loci_name, directory = paste0(category_output_folder, v_IQTree_name, "/"),
                                    original_alignment_folder = alignment_dir[[dataset]], add.charpartition = FALSE)
    }
    # Create a record of which loci went into which analysis
    output_text <- v_df$loci_name
    write(output_text, file = v_text_name)
    # Add the number of loci in this category to the summary row
    summary_row <- c(summary_row, length(v_df$loci_name))
  }
  
  # Apply all four tests
  # Create output names
  all_text_name <- paste0(text_records_dir, dataset, "_all_loci_record.txt")
  all_ASTRAL_name <- paste0(dataset,"_all_ASTRAL")
  if (estimate.species.trees.in.IQTREE == TRUE){
    if (partition.by.codon.position == TRUE){
      all_IQTree_name <- paste0(dataset,"_all_IQTREE_partitioned")
    } else if (partition.by.codon.position == FALSE){
      all_IQTree_name <- paste0(dataset,"_all_IQTREE") 
    }
  }
  
  ## Subset the dataframe_df to loci that pass all four tests
  all_df <- dataset_df[((dataset_df$X3SEQ_p_value > 0.05) & (dataset_df$PHI_normal_p_value > 0.05) & (dataset_df$max_chi_squared_p_value > 0.05) & (dataset_df$geneconv_inner_fragment_simulated_p_value > 0.05)) ,]
  # Copy loci trees for ASTRAL
  copy.loci.trees(all_df$loci_name, all_df$tree, category_output_folder, all_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
  # Copy loci alignments for IQ-Tree
  if (estimate.species.trees.in.IQTREE == TRUE){
    # create the partition file required to run this IQ-Tree analysis
    # partition.file.from.loci.list(loci_list, directory, original_alignment_folder, add.charpartition = FALSE)
    partition.file.from.loci.list(loci_list = all_df$loci_name, directory = paste0(category_output_folder, all_IQTree_name, "/"),
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition = FALSE)
  }
  # Create a record of which loci went into which analysis
  output_text <- all_df$loci_name
  write(output_text, file = all_text_name)
  # Add the number of loci in this category to the summary row
  summary_row <- c(summary_row, length(all_df$loci_name))
  
  
  ## Apply no tests
  NoTest_text_name <- paste0(text_records_dir, dataset, "_NoTest_loci_record.txt")
  NoTest_ASTRAL_name <- paste0(dataset,"_NoTest_ASTRAL")
  if (estimate.species.trees.in.IQTREE == TRUE){
    if (partition.by.codon.position == TRUE){
      NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE_partitioned")
    } else if (partition.by.codon.position == FALSE){
      NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE") 
    }
  }
  # Copy loci trees for ASTRAL
  copy.loci.trees(dataset_df$loci_name, dataset_df$tree, category_output_folder, NoTest_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
  # Copy loci alignments for IQ-Tree
  if (estimate.species.trees.in.IQTREE == TRUE){
    # create the partition file required to run this IQ-Tree analysis
    # partition.file.from.loci.list(loci_list, directory, original_alignment_folder, add.charpartition = FALSE)
    partition.file.from.loci.list(loci_list = dataset_df$loci_name, directory = paste0(category_output_folder, NoTest_IQTree_name, "/"),
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition = FALSE)
  }
  # Create a record of which loci went into which analysis
  output_text <- dataset_df$loci_name
  write(output_text, file = NoTest_text_name)
  # Add the number of loci in this category to the summary row
  summary_row <- c(summary_row, length(dataset_df$loci_name))
  
  
  ## Write out the summary row as a dataframe
  names(summary_row) <- c("dataset", "n_pass_3SEQ", "n_pass_PHI", "n_pass_maxchi", "n_NA_geneconv", "n_pass_geneconv", "n_pass_all", "n_NoTest")
  summary_df <- data.frame(as.list(summary_row))
  summary_op_file <- paste0(output_dirs[dataset], dataset, "_species_tree_summary.csv")
  write.csv(summary_df, file = summary_op_file, row.names = FALSE)
}


### Estimate a species tree for each of the five categories
for (dataset in datasets_to_estimate_trees){
  # Ensure the folder for species trees data exists
  category_output_folder <- paste0(output_dirs[dataset], "species_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  
  # Get list of all files in that folder
  all_category_folder_files <- list.files(category_output_folder)
  # Filter into ASTRAL text files and IQ-Tree folders
  astral_files <- paste0(grep("\\.txt", grep("ASTRAL", all_category_folder_files, value = TRUE), value = TRUE))
  iqtree_files <- paste0(grep("\\.", grep("IQTREE", all_category_folder_files, value = TRUE), value = TRUE, invert = TRUE))
  # Contruct names of finished treefiles
  astral_files_finished_names <- paste0(category_output_folder, gsub(".txt", "_species.tre", astral_files))
  iqtree_files_finished_names <- paste0(category_output_folder, iqtree_files, ".contree")
  
  # Identify which ASTRAL analyses have not been run
  astral_files_to_run <- astral_files[!file.exists(astral_files_finished_names)]
  # Run remaining ASTRAL analyses
  if (length(astral_files_to_run) > 0){
    # Construct full file path
    astral_files_to_run <- paste0(category_output_folder, astral_files_to_run)
    # Estimate the species trees using ASTRAL
    lapply(astral_files_to_run, ASTRAL.wrapper, exec_paths["ASTRAL"])
  }
  
  # Identify which IQTREE analyses have not been run
  iqtree_files_to_run <- iqtree_files[!file.exists(iqtree_files_finished_names)]
  # Run remaining IQ-Tree analyses
  if (length(iqtree_files_to_run) > 0){
    # Construct full file path
    iqtree_files_to_run <- paste0(category_output_folder, iqtree_files_to_run)
    # Estimate the species trees using IQ-Tree
    if (estimate.species.trees.in.IQTREE == TRUE){
      if (partition.by.codon.position == TRUE){
        # Estimate the species tree on each folder of alignments using the partition file
        partitions_to_run <- paste0(dirname(iqtree_files_to_run), "/", basename(iqtree_files_to_run), "/", "partitions.nex")
        lapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"])
      } else if (partition.by.codon.position == FALSE){
        # If not partitioning data by codon position, simply call IQ-Tree on the folder of alignments to estimate a tree
        lapply(iqtree_files_to_run, estimate.IQTREE.species.tree, exec_paths["IQTree"])
      }
    }
  }
  
}
