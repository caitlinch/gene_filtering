### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from treelike or non-treelike loci
## Additional software packages are required:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (need version 2.0 or later)
# Caitlin Cherryh 2021

##### Step 1: Set file paths and run variables #####
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first and same dataset last)
# alignment_dir     <- the folder(s) containing the alignments for each loci
# csv_data_dir      <- directory containing the .csv file results from script 1_RecombinationDetection_empiricalTreelikeness.R
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# exec_paths        <- location to the software executables needed for analysis (ASTRAL and IQTREE)
# datasets_to_copy_loci   <- Out of the input names, select which datasets to copy loci trees for tree estimation
# datasets_to_estimate_trees <- Out of the input names, select which datasets to estimate species trees based on treelikeness results
# partition.by.codon.position <- Whether to run analysis partitioning by codon position
#                             <- set TRUE if you want to estimate species trees partitioning by codon position, and FALSE if you don't
#                             <- codon position simply counts every third base starting from 1st, 2nd, or 3rd base, and does not account for frame shift
# use.modelfinder.models.for.partitions <- can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. 
#                                       <- TRUE will use "-m MERGE" and include a charpartition with substitution models selected by ModelFinder in IQ-Tree

# # To run this program: 
# # 1. Delete the lines below that include Caitlin's paths/variables
# # 2. Uncomment lines 24 to 41 inclusive and fill with your own variable names
# input_names <- ""
# alignment_dir <- ""
# csv_data_dir <- ""
# output_dir <- ""
# maindir <- ""
# # Create a vector with all of the executable file paths:
# # exec_paths <- c("/path/to/ASTRAL_executable","/path/to/IQ-Tree_executable")
# # names(exec_paths) <- c("ASTRAL","IQTree")
# # To access a path: exec_paths[["name"]]
# exec_paths <- c()
# names(exec_paths) <- c("ASTRAL","IQTree")
# datasets_to_copy_loci <-  c()
# datasets_to_estimate_trees <- c()
# estimate.species.trees.in.IQTREE = TRUE
# partition.by.codon.position = FALSE

### Caitlin's paths ###
run_location = "server"

if (run_location == "local"){
  # Datasets/dataset information
  input_names <- c( "1KP", "Strassert2021", "Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths  in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("Astral/astral.5.7.5.jar","iqtree-2.0-rc1-MacOSX/bin/iqtree")
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # Select number of cores for parallelisation
  cores.to.use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c("Pease2016", "Vanderpool2020", "1KP", "Strassert2021")
  datasets_to_estimate_trees <- c("Pease2016", "Vanderpool2020", "1KP", "Strassert2021")
  estimate.species.trees.in.IQTREE = TRUE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd - based on position in alignment file) 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use "-m MERGE" and include a charpartition with substitution models
  
} else if (run_location=="server"){
  # Datasets/dataset information
  input_names <- c( "1KP", "Strassert2021", "Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/data/caitlin/empirical_treelikeness/Data_1KP/",
                     "/data/caitlin/empirical_treelikeness/Data_Strassert2021/",
                     "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/",
                     "/data/caitlin/empirical_treelikeness/Data_Pease2016/")
  csv_data_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/executables/ASTRAL/astral.5.7.5.jar",
                  "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # Select number of cores for parallelisation
  cores.to.use = 45
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci <-  c()
  datasets_to_estimate_trees <- c("Pease2016", "Vanderpool2020", "1KP", "Strassert2021")
  estimate.species.trees.in.IQTREE = TRUE # can be TRUE of FALSE - if TRUE, will run IQ-Tree analyses
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use "-m MERGE" and include a charpartition with substitution models
}
### End of Caitlin's paths ###



##### Step 2: Open packages and source files for functions #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
# Source the functions using the filepaths
source(paste0(maindir,"code/func_empirical.R"))
source(paste0(maindir,"code/func_analysis.R"))



##### Step 3: Set up input and output folders #####
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



##### Step 4: Assemble the treelikeness dataframe #####
if (length(datasets_to_copy_loci) > 0){
  # Check whether a collated, trimmed recombination detection results file exists
  trimmed_treelikeness_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci), collapse="_"), "_collated_RecombinationDetection_TrimmedLoci.csv")
  pass_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci), collapse="_"), "_RecombinationDetection_PassFail_record.csv")
  collated_exclude_file <- paste0(csv_data_dir, "01_IQ-Tree_warnings_",paste(sort(datasets_to_copy_loci), collapse="_"), "_LociToExclude.csv")
  
  if (file.exists(trimmed_treelikeness_df_file) & file.exists(pass_df_file)){
    treelikeness_df <- read.csv(trimmed_treelikeness_df_file, stringsAsFactors = TRUE)
    pass_df <- read.csv(pass_df_file, stringsAsFactors = TRUE)
  } else{
    # If the file doesn't exist, create it
    # Get a list of all the csv files in the csv_data_directory
    all_files <- list.files(csv_data_dir)
    # Get the results filenames for the datasets of interest
    all_results <- grep("RecombinationDetection_complete_collated_results", all_files, value = TRUE)
    all_results <- paste0(csv_data_dir, all_results)
    results <- c()
    for (dataset in datasets_to_copy_loci){
      f <- grep(dataset, all_results, value = TRUE)
      results <- c(results, f)
    }
    # Remove duplicates
    results <- unique(results)
    # Open and attach the datasets
    treelikeness_df <- as.data.frame(do.call(rbind, lapply(results, read.csv)))
    treelikeness_df$match <- paste0(treelikeness_df$dataset, ":", treelikeness_df$loci_name)
    # If the collated total file hasn't been saved, save it
    all_treelikeness_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci), collapse="_"), "_collated_RecombinationDetection.csv")
    if (file.exists(all_treelikeness_file) == FALSE){
      write.csv(treelikeness_df, all_treelikeness_file)
    }
    
    # Open the csv containing the list of loci to exclude from species tree analysis
    exclude_file <- paste0(csv_data_dir, grep("LociToExclude.csv", all_files, value = TRUE))
    # Open the csv files and bind into one dataframe
    exclude_df <- as.data.frame(do.call(rbind, lapply(exclude_file, read.csv)))
    # Collate and output if that file doesn't exist
    if (file.exists(collated_exclude_file) == FALSE){
      write.csv(exclude_df, file = collated_exclude_file, row.names = FALSE)
    }
    
    # Reduce down to the unique dataset/loci pairs to exclude from treelikeness_df
    exclude_df <- exclude_df[,c("dataset", "loci")]
    exclude_df <- exclude_df[!duplicated(exclude_df),]
    exclude_df$match <- paste0(exclude_df$dataset, ":", exclude_df$loci)
    
    # Identify the row numbers of the loci to exclude
    exclude_loci <- which(treelikeness_df$match %in% exclude_df$match)
    # Compare with the set of all loci in all datasets to get the loci to keep
    include_loci <- setdiff(1:nrow(treelikeness_df), exclude_loci)
    
    # Trim the treelikeness_df to remove the loci to be excluded
    treelikeness_df <- treelikeness_df[include_loci,]
    
    # Create new columns with pass/fail for each test
    treelikeness_df <- make.pass.fail.column("pass_3seq", "X3SEQ_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_phi", "PHI_normal_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_phi_permute", "PHI_permutation_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_maxchi", "max_chi_squared_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_NSS", "NSS_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_geneconv_inner", "geneconv_inner_fragment_simulated_p_value", treelikeness_df)
    treelikeness_df <- make.pass.fail.column("pass_geneconv_outer", "geneconv_outer_fragment_simulated_p_value", treelikeness_df)
    # Create a column for GeneConv where both inner and outer fragment p-value must be >0.05 to pass
    treelikeness_df$pass_geneconv <- "FALSE"
    treelikeness_df$pass_geneconv[((treelikeness_df$geneconv_outer_fragment_simulated_p_value > 0.05) & 
                                     (treelikeness_df$geneconv_inner_fragment_simulated_p_value > 0.05))] <- "TRUE"
    treelikeness_df$pass_geneconv <- as.logical(treelikeness_df$pass_geneconv)
    
    # Save the trimmed treelikeness_df
    trimmed_treelikeness_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci), collapse="_"), "_collated_RecombinationDetection_TrimmedLoci.csv")
    write.csv(treelikeness_df, file = trimmed_treelikeness_df_file, row.names = FALSE)
    # Save a df of just the pass/fail info
    pass_df <- treelikeness_df[,c("dataset", "loci_name", "alphabet", "n_taxa", "n_bp", "pass_3seq", "pass_phi", "pass_maxchi", 
                                  "pass_NSS", "pass_geneconv_inner", "pass_geneconv_outer", "pass_geneconv")]
    pass_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci), collapse="_"), "_RecombinationDetection_PassFail_record.csv")
    write.csv(pass_df, file = pass_df_file, row.names = FALSE)
  }
}



##### Step 5: Categorize loci by test results and prepare subsets of data for tree estimation #####
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
  
  ### Apply a single recombination detection test ###
  ## Iterate through each var and save the loci/trees for a tree made from only loci that pass the test
  # make a list of the variables on which to filter the loci - should be columns from the treelikeness_df
  vars <- c("pass_phi", "pass_maxchi", "pass_geneconv")
  # Assign output names for each of the variables in vars
  vars_names <- c("PHI", "maxchi", "geneconv")
  names(vars_names) <- vars
  # Iterate through each var: 
  for (v in vars){
    # Make a tree for all the loci that pass the test and all the loci that fail the test
    tree_type <- c("pass", "fail")
    for (tt in tree_type){
      print(paste0(dataset, " : ", v, " : ", tt))
      # Set boolean to collect either loci that passed or failed the test
      if (tt == "pass"){
        bool = TRUE
      } else if (tt == "fail"){
        bool = FALSE
      }
      
      # Get short version of name for output files
      v_name <- vars_names[v]
      # Make names for output files
      v_text_name <- paste0(text_records_dir, dataset, "_", v_name, "_", tt, "_loci_record.txt")
      v_ASTRAL_name <- paste0(dataset,"_",v_name,"_", tt, "_ASTRAL")
      if (estimate.species.trees.in.IQTREE == TRUE){
        if (partition.by.codon.position == TRUE){
          v_IQTree_name <- paste0(dataset,"_",v_name, "_", tt, "_IQTREE_partitioned")
        } else if (partition.by.codon.position == FALSE){
          v_IQTree_name <- paste0(dataset,"_",v_name, "_", tt, "_IQTREE") 
        }
      }
      
      # Remove any loci that have an NA result for this test
      na_test_inds <- which(is.na(dataset_df[[v]]))
      dataset_inds <- 1:nrow(dataset_df)
      keep_inds <- setdiff(dataset_inds, na_test_inds)
      v_df <- dataset_df[keep_inds,]
      
      # Break up dataframe into only loci that pass the test
      v_inds <- which(dataset_df[,c(v)] == bool)
      # Use the indexes to subset the dataframe to just loci that pass the test 
      # (i.e. have a non significant p-value, meaning the null hypothesis of treelikeness cannot be rejected)
      v_df <- dataset_df[v_inds,]
      
      # Copy trees of all loci that pass the test into one file that can be fed into ASTRAL
      copy.loci.trees(v_df$loci_name, v_df$tree, category_output_folder, v_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
      # If running IQ-Tree analysis, copy all loci into a separate folder that can be fed into IQ-Tree
      if (estimate.species.trees.in.IQTREE == TRUE){
        # create the partition file required to run this IQ-Tree analysis
        partition.file.from.loci.list(loci_list = v_df$loci_name, directory = paste0(category_output_folder, v_IQTree_name, "/"),
                                      original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                      substitution_models = v_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
      }
      
      # Create a record of which loci went into which analysis
      output_text <- v_df$loci_name
      write(output_text, file = v_text_name)
      # Add the number of loci in this category to the summary row
      summary_row <- c(summary_row, length(v_df$loci_name))
    }
  }
  
  
  
  ### Apply all three tests - get loci that pass all and those that fail all ###
  ## Remove loci that have an NA result for any test - do not include these in any analysis for comparative purposes
  # identify loci that have an NA for any one test
  na_phi_inds <- which(is.na(dataset_df$pass_phi))
  na_maxchi_inds <- which(is.na(dataset_df$pass_maxchi))
  na_geneconv_inds <- which(is.na(dataset_df$pass_geneconv))
  na_dataset_inds <- sort(unique(c(na_phi_inds, na_maxchi_inds, na_geneconv_inds)))
  dataset_inds <- 1:nrow(dataset_df)
  allTest_inds <- setdiff(dataset_inds, na_dataset_inds)
  # Reduce allTest_df to loci without na for any test
  allTest_df <- dataset_df[allTest_inds, ]
  
  # 
  tree_type <- c("pass", "fail")
  for (tt in tree_type){
    print(paste0(dataset, " : all tests : ", tt))
    
    # Create output names
    all_text_name <- paste0(text_records_dir, dataset, "_allTests_",tt ,"_loci_record.txt")
    all_ASTRAL_name <- paste0(dataset,"_allTests_", tt, "_ASTRAL")
    if (estimate.species.trees.in.IQTREE == TRUE){
      if (partition.by.codon.position == TRUE){
        all_IQTree_name <- paste0(dataset,"_allTests_", tt, "_IQTREE_partitioned")
      } else if (partition.by.codon.position == FALSE){
        all_IQTree_name <- paste0(dataset,"_allTests_", tt, "_IQTREE") 
      }
    }
    
    ## Subset the dataframe_df to loci that pass all three tests
    # Select loci that pass/fail all three tests
    if (tt == "pass"){
      all_df <- allTest_df[((allTest_df$pass_phi == TRUE) & (allTest_df$pass_maxchi == TRUE) & (allTest_df$pass_geneconv == TRUE)), ]
    } else if (tt == "fail"){
      all_df <- allTest_df[((allTest_df$pass_phi == FALSE) & (allTest_df$pass_maxchi == FALSE) & (allTest_df$pass_geneconv == FALSE)), ] 
    }
    # Copy loci trees for ASTRAL
    copy.loci.trees(all_df$loci_name, all_df$tree, category_output_folder, all_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
    # Copy loci alignments for IQ-Tree
    if (estimate.species.trees.in.IQTREE == TRUE){
      # create the partition file required to run this IQ-Tree analysis
      partition.file.from.loci.list(loci_list = all_df$loci_name, directory = paste0(category_output_folder, all_IQTree_name, "/"),
                                    original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                    substitution_models = all_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
    }
    # Create a record of which loci went into which analysis
    output_text <- all_df$loci_name
    write(output_text, file = all_text_name)
    # Add the number of loci in this category to the summary row
    summary_row <- c(summary_row, length(all_df$loci_name))
  }
  
  ### Apply no tests ###
  print(paste0(dataset, " : No tests -- include all loci from modified dataset_df"))
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
    partition.file.from.loci.list(loci_list = dataset_df$loci_name, directory = paste0(category_output_folder, NoTest_IQTree_name, "/"),
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = FALSE,
                                  substitution_models = dataset_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
  }
  
  ### Create a record of which loci went into which analysis ###
  output_text <- dataset_df$loci_name
  write(output_text, file = NoTest_text_name)
  # Add the number of loci in this category to the summary row
  summary_row <- c(summary_row, length(dataset_df$loci_name))
  
  # Expand the summary row to identify the number of loci that pass or fail the same tests
  n_pass_PHI_maxchi <- nrow(dataset_df[((dataset_df$pass_phi == TRUE) & (dataset_df$pass_maxchi == TRUE)), ])
  n_fail_PHI_maxchi <- nrow(dataset_df[((dataset_df$pass_phi == FALSE) & (dataset_df$pass_maxchi == FALSE)), ])
  n_pass_PHI_geneconv <- nrow(dataset_df[((dataset_df$pass_phi == TRUE) & (dataset_df$pass_geneconv == TRUE)), ])
  n_fail_PHI_geneconv <- nrow(dataset_df[((dataset_df$pass_phi == FALSE) & (dataset_df$pass_geneconv == FALSE)), ])
  n_pass_maxchi_geneconv <- nrow(dataset_df[((dataset_df$pass_maxchi == TRUE) & (dataset_df$pass_geneconv == TRUE)), ])
  n_fail_maxchi_geneconv <- nrow(dataset_df[((dataset_df$pass_maxchi == FALSE) & (dataset_df$pass_geneconv == FALSE)), ])
  summary_row <- c(summary_row, n_pass_PHI_maxchi, n_fail_PHI_maxchi, n_pass_PHI_geneconv, n_fail_PHI_geneconv, n_pass_maxchi_geneconv, n_fail_maxchi_geneconv)
  
  ### Write out the summary row as a dataframe ###
  names(summary_row) <- c("dataset", "n_pass_PHI", "n_fail_PHI", "n_pass_maxchi", "n_fail_maxchi", "n_pass_geneconv", "n_fail_geneconv",
                          "n_pass_allTests", "n_fail_allTests", "n_NoTest", "n_pass_PHI_maxchi", "n_fail_PHI_maxchi", "n_pass_PHI_geneconv", 
                          "n_fail_PHI_geneconv", "n_pass_maxchi_geneconv", "n_fail_maxchi_geneconv")
  summary_df <- data.frame(as.list(summary_row))
  summary_op_file <- paste0(output_dirs[dataset], dataset, "_species_tree_summary.csv")
  write.csv(summary_df, file = summary_op_file, row.names = FALSE)
}




##### Step 6: Estimate species trees #####
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
    mclapply(astral_files_to_run, ASTRAL.wrapper, exec_paths["ASTRAL"], mc.cores = cores.to.use)
  }
  
  # Identify which IQTREE analyses have not been run
  iqtree_files_to_run <- iqtree_files[!file.exists(iqtree_files_finished_names)]
  # Run remaining IQ-Tree analyses
  if (length(iqtree_files_to_run) > 0){
    # Construct full file path
    iqtree_files_to_run <- paste0(category_output_folder, iqtree_files_to_run)
    # Estimate the species trees using IQ-Tree
    if (estimate.species.trees.in.IQTREE == TRUE){
      # Estimate the species tree on each folder of alignments using the partition file
      partitions_to_run <- paste0(dirname(iqtree_files_to_run), "/", basename(iqtree_files_to_run), "/", "partitions.nex")
      if (dataset == "Vanderpool2020" | dataset == "Strassert2021"){
        mclapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"], 
                 IQTREE_model_command = "MFP+MERGE", mc.cores = cores.to.use)
      } else if (dataset = "1KP" | dataset = "Pease2016"){
        mclapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"],
                 IQTREE_model_command = "MERGE", mc.cores = cores.to.use)
      }
    }
  }
  
}
