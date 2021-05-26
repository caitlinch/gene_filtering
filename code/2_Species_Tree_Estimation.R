### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from loci with varying treelikeness
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
# Caitlin Cherryh 2021

##### Step 1: Set file paths and run variables #####
print("set filepaths")
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
  input_names <- c("1KP", "Strassert2021","Vanderpool2020")
  loci_to_remove <- list("Vanderpool2020" = "ORTHOMCL14552")
  number_of_taxa <- list("Vanderpool2020" = NA)
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  treelikeness_df_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/empiricalTreelikeness_Vanderpool2020_collated_results_20210526.csv"
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
  
  # Set number of replicates to do for the treelikeness category analysis
  number_of_category_replicates = 999
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("1KP", "Strassert2021","Vanderpool2020")
  datasets_to_estimate_trees <- c("1KP", "Strassert2021","Vanderpool2020")
  check.for.warnings = FALSE # check IQ-Tree .log file and .iqtree file output for each gene tree for warnings
  estimate.trees.in.IQTREE = FALSE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = TRUE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
  
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
  
  # Set number of replicates to do for the treelikeness category analysis
  number_of_category_replicates = 999
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("1KP", "Strassert2021","Vanderpool2020")
  datasets_to_estimate_trees <- c("1KP", "Strassert2021","Vanderpool2020")
  check.for.warnings = FALSE # check IQ-Tree .log file and .iqtree file output for each gene tree for warnings
  estimate.trees.in.IQTREE = FALSE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = TRUE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
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
# Save trimmed treelikeness_df
trimmed_treelikeness_df_file <- gsub(".csv", paste0("_trimmedLoci_trimmedTaxa_",format(Sys.Date(),"%Y%m%d"),".csv"), treelikeness_df_file)
write.csv(treelikeness_df, file = trimmed_treelikeness_df_file)



##### Step 4: Categorise loci by treelikeness test p-values (3seq and tree proportion) #####
# Iterate through each of the datasets
# Sort loci into four categories based on treelikeness p-values
#     * both 3seq and tree proportion significant
#     * neither 3seq nor tree proportion significant
#     * only 3seq significant
#     * only tree proportion significant
# Save the trees from a category into a separate folder and collect all trees from a category into a collated text file
for (dataset in datasets_to_copy_loci){
  # Create a new folder to put these 400 text files (per dataset) in
  category_output_folder <- paste0(output_dirs[dataset], "ASTRAL_category_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  # filter treelikeness_df by dataset
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  ### split loci into four groups (neither, 3seq, tp or both), then copy all loci alignments from each group into a new folder and all trees from each group into a new collated text file
  ### Note that these categories are based on the number of statistical tests that are statistically significant
  ### e.g. none = no loci in this category have a statistically significant p-value, which means we accept the null hypotheses of treelikeness and clonal evolution)
  # 3seq p-value and tree proportion p-value both >0.05 (not significant = do not reject null hypothesis of treelikeness/clonal evolution)
  # none = NO SIGNIFICANT P-VALUES
  cat_none <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
  copy.loci.trees(cat_none, dataset_df[dataset_df$loci %in% cat_none,]$tree, category_output_folder, 
                  "p-value_categories_none_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  # 3seq p-value and tree proportion p-value both <=0.05 (significant = reject null hypothesis of treelikeness/clonal evolution)
  # both = BOTH SIGNIFICANT P-VALUES
  cat_both <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_both, dataset_df[dataset_df$loci %in% cat_both,]$tree, category_output_folder, 
                  "p-value_categories_both_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  # Only 3SEQ p-value <=0.05 and significant, tree proportion p-value not significant  (reject null hypothesis of clonal evolution but do not reject null hypothesis of treelikeness)
  # 3seq only = treelike, not clonal
  cat_3seq <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
  copy.loci.trees(cat_3seq, dataset_df[dataset_df$loci %in% cat_3seq,]$tree, category_output_folder, 
                  "p-value_categories_3seq_only_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant (reject null hypothesis of treelikeness but do not reject null hypothesis of clonal evolution)
  # tree_proportion_only = not-treelike, clonal
  cat_tp   <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
  copy.loci.trees(cat_tp, dataset_df[dataset_df$loci %in% cat_tp,]$tree, category_output_folder, 
                  "p-value_categories_tree_proportion_only_ASTRAL", copy.all.individually = FALSE, copy.and.collate = TRUE)
  
  # Write out the loci from each category
  write(cat_none, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_none_loci_list.txt"))
  write(cat_both, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_both_loci_list.txt"))
  write(cat_3seq, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_3seq_only_loci_list.txt"))
  write(cat_tp, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_tree_proportion_only_loci_list.txt"))
  cat_op_df <- data.frame("dataset" = dataset, "n_loci_none" = length(cat_none), "n_plot_name_none" = "Treelike",
                          "n_loci_both" = length(cat_both), "n_plot_name_both" = "Non-treelike 1",
                          "n_loci_3seq_only" = length(cat_3seq), "n_plot_name_3seq_only" = "Non-treelike 3",
                          "n_loci_tree_proportion_only" = length(cat_tp), "n_plot_name_tree_proportion_only" = "Non-treelike 2")
  write.csv(cat_op_df, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_info.csv"))
  
  # Perform 99 replicates of each category with randomly sampled loci
  # Set information about categories
  category_name <- c("cat_none","cat_both","cat_3seq","cat_tp")
  category_list <- list("cat_none" = cat_none, "cat_both" = cat_both, "cat_3seq" = cat_3seq, "cat_tp" = cat_tp)
  category_output_names <- c("none","both","3seq_only","tree_proportion_only")
  # Iterate through the conditions
  for (c in 1:length(category_name)){
    # Select information about this category by index
    c_name <- category_name[c]
    c_loci <- category_list[[c]]
    c_op_name <- category_output_names[c]
    # Initialise list to store loci replicates
    rep_list <- list()
    # Perform 99 replicates of this condition
    for (i in 1:number_of_category_replicates){
      # Identify the number of loci included in this category
      n_loci <- length(c_loci)
      # Randomly sample the list of loci in the dataset_df
      replicate_loci <- sample(dataset_df$loci, n_loci)
      # Pad out the number for a nice output name
      rep_id <- sprintf("%04d",i)
      # Copy the trees into a separate file
      copy.loci.trees(replicate_loci, dataset_df[dataset_df$loci %in% replicate_loci,]$tree, category_output_folder, 
                      paste0("p-value_categories_", c_op_name,"_ASTRAL" ,"_replicate", rep_id), copy.all.individually = FALSE, copy.and.collate = TRUE)
      # Add loci selected for this replicate to the rep_list
      rep_list[[as.character(rep_id)]] <- replicate_loci
    }
    rep_df <- as.data.frame(do.call(rbind, rep_list))
    names(rep_df) <- paste0("sampled_loci_",sprintf("%04d", 1:ncol(rep_df)))
    write.csv(rep_df, file = paste0(output_dirs[dataset], dataset, "_p-value_categories_", c_op_name, "_ASTRAL_SampledReplicatesLoci.csv"))
  }
  
  # If estimating trees in IQ-Tree, copy the information you need into a new folder for each analysis
  if (estimate.trees.in.IQTREE == TRUE){
    if (partition.by.codon.position == FALSE){
      mclapply(cat_none, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_none_IQ-Tree/"), mc.cores = cores_to_use)
      mclapply(cat_both, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_both_IQ-Tree/"), mc.cores = cores_to_use)
      mclapply(cat_3seq, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_3seq_only_IQ-Tree/"), mc.cores = cores_to_use)
      mclapply(cat_tp, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_tree_proportion_only_IQ-Tree/"), mc.cores = cores_to_use)
      # If partitioning by codon position, repeat the IQ-Tree set-up
    } else if (partition.by.codon.position == TRUE){
      # filter treelikeness_df by dataset
      dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
      # split loci into four groups (neither, 3seq, tp or both), then copy all loci alignments from each group into a new folder and all trees from each group into a new collated text file
      # 3seq p-value and tree proportion p-value both >0.05 (not significant)
      cat_none <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
      mclapply(cat_none, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_none_IQ-Tree_partition/"), mc.cores = cores_to_use)
      # 3seq p-value and tree proportion p-value both <=0.05 (significant)
      cat_both <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
      mclapply(cat_both, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_both_IQ-Tree_partition/"), mc.cores = cores_to_use)
      # Only 3seq p-value <=0.05 and significant, tree proportion p-value not significant
      cat_3seq <- dataset_df[dataset_df$X3SEQ_p_value <= 0.05 & dataset_df$tree_proportion_p_value > 0.05,]$loci
      mclapply(cat_3seq, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_3seq_only_IQ-Tree_partition/"), mc.cores = cores_to_use)
      # Only tree proportion p-value <=0.05 and significant, 3seq p-value not significant
      cat_tp   <- dataset_df[dataset_df$X3SEQ_p_value > 0.05 & dataset_df$tree_proportion_p_value <= 0.05,]$loci
      mclapply(cat_tp, copy.loci.alignment, alignment_dir[dataset], paste0(category_output_folder,"p-value_categories_tree_proportion_only_IQ-Tree_partition/"), mc.cores = cores_to_use)
    }
  }
}

# Estimate a species tree for each of the four categories
for (dataset in datasets_to_estimate_trees){
  # Create a new folder to put these 400 text files (per dataset) in
  category_output_folder <- paste0(output_dirs[dataset], "ASTRAL_category_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  # Get list of all files in that folder
  all_category_folder_files <- list.files(category_output_folder)
  # Filter by text files
  text_category_folder_files <- grep(".txt", all_category_folder_files, value = TRUE)
  tree_category_folder_files <- grep("p-value_categories", text_category_folder_files, value = TRUE)
  run_astral_files <- paste0(category_output_folder, grep("ASTRAL", tree_category_folder_files, value = TRUE))
  # Calculate the species tree using ASTRAL for each of the four categories
  # estimate.ASTRAL.species.tree(astral_inputs[1], gsub(".txt","_ASTRAL_species.tre",astral_inputs[1]), gsub(".txt","_ASTRAL_species.log",astral_inputs[1]), exec_paths["ASTRAL"]
  lapply(run_astral_files, ASTRAL.wrapper, exec_paths["ASTRAL"])
  
  if (estimate.trees.in.IQTREE == TRUE){
    # Get list of all files in that folder
    all_category_folder_files <- list.files(category_output_folder)
    iqtree_files <- grep("IQ-Tree", all_category_folder_files, value = TRUE)
    partition_files <- grep("partition", iqtree_files, value = TRUE)
    unpartitioned_files <- grep("partition", iqtree_files, value = TRUE, invert = TRUE)
    if (partition.by.codon.position == FALSE){
      # Calculate the species tree using IQ-Tree for each of the four categories
      run_iqtree_files <- paste0(category_output_folder, unpartitioned_files, "/")
      lapply(run_iqtree_files, estimate.IQTREE.species.tree, exec_paths["IQTree"])
      # If partitioning by codon position, create a partition file for each folder and then estimate the tree using IQ-Tree
    } else if (partition.by.codon.position == TRUE){
      # Make list of folders for partition analysis
      run_iqtree_files <- paste0(category_output_folder, partition_files, "/")
      # Create the partition files
      mclapply(run_iqtree_files, make.partition.file, mc.cores = cores_to_use)
      # Run the analysis
      partitions_to_run <- paste0(run_iqtree_files, "partitions.nex")
      lapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"])
    }
  }
}



##### Step 5: Species tree with all loci #####
for (dataset in datasets_to_copy_loci){
  # Subset dataframe to only this dataset
  dataset_df <- treelikeness_df[treelikeness_df$dataset == dataset,]
  # Create name for ASTRAL run
  astral_name <- "all_loci_ASTRAL"
  # Copy all loci and all trees 
  all_loci <- dataset_df$loci
  all_trees <- dataset_df$tree
  copy.loci.trees(all_loci, all_trees, output_dirs[dataset], astral_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
  
  if (estimate.trees.in.IQTREE == TRUE){
    if (partition.by.codon.position == FALSE){
      # Create name for IQ-Tree run
      iqtree_name <- "all_loci_IQ-Tree/"
      # Copy all loci
      mclapply(all_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], iqtree_name), mc.cores = cores_to_use)
    } else if (partition.by.codon.position == TRUE){
      # Create name for IQ-Tree run
      partition_name <- "all_loci_IQ-Tree_partition/"
      # Copy all loci
      mclapply(all_loci, copy.loci.alignment, alignment_dir[dataset], paste0(output_dirs[dataset], partition_name), mc.cores = cores_to_use)
    }
  }
}

for (dataset in datasets_to_estimate_trees){
  # Calculate the species trees for all loci
  # estimate.ASTRAL.species.tree(astral_inputs[1], gsub(".txt","_ASTRAL_species.tre",astral_inputs[1]), gsub(".txt","_ASTRAL_species.log",astral_inputs[1]), exec_paths["ASTRAL"]
  ASTRAL.wrapper(paste0(output_dirs[[dataset]], astral_name,".txt"), exec_paths[["ASTRAL"]])
  
  if (estimate.trees.in.IQTREE == TRUE){
    if (partition.by.codon.position == FALSE){
      # Calculate the species tree using IQ-Tree for each of the four categories
      estimate.IQTREE.species.tree(paste0(output_dirs[dataset],iqtree_name), exec_paths[["IQTree"]])
    } else if (partition.by.codon.position == TRUE){
      partition_name <- "all_loci_IQTree_partition/"
      make.partition.file(paste0(output_dirs[dataset], partition_name))
      estimate.partitioned.IQTREE.species.tree(paste0(output_dirs[[dataset]], partition_name, "partitions.nex"), exec_paths[["IQTree"]])
    } 
  }
}


