### empirical_treelikeness/code/2_TreeEstimation_EmpiricalData.R
## R program to estimate trees from treelike or non-treelike loci
## Additional software packages are required:
##     - ASTRAL (Zhang et al 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/) (version 2.0 or later)
# Caitlin Cherryh 2021

##### Step 1: Set file paths and run variables #####
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first and same dataset last)
# alignment_dir     <- the folder(s) containing the alignments for each loci
# csv_data_dir      <- directory containing the .csv file results from script 1_RecombinationDetection_empiricalTreelikeness.R
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# exec_paths        <- location to the software executables needed for analysis (ASTRAL and IQTREE)
# datasets_to_copy_loci_ASTRAL_IQTREE <- Out of the input names, select which datasets to copy loci trees for tree estimation in ASTRAL or IQ-Tree
# datasets_to_copy_loci_RAxML <- Out of the input names, select which datasets to copy loci trees for tree estimation in RAxML-NG
# datasets_to_estimate_ASTRAL_trees <- Out of the input names, select which datasets to estimate species trees in ASTRAL based on recombination results
# datasets_to_estimate_IQTREE_trees <- Out of the input names, select which datasets to estimate species trees in IQ-Tree based on recombination results
# datasets_to_estimate_RAxML_trees <- Out of the input names, select which datasets to estimate species trees in IQ-Tree based on recombination results
# partition.by.codon.position <- Whether to run analysis partitioning by codon position
#                             <- set TRUE if you want to estimate species trees partitioning by codon position, and FALSE if you don't
#                             <- codon position simply counts every third base starting from 1st, 2nd, or 3rd base, and does not account for frame shift
# use.modelfinder.models.for.partitions <- can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. 
#                                       <- TRUE will use "-m MERGE" and include a charpartition with substitution models selected by ModelFinder in IQ-Tree

# We estimated species trees in ASTRAl for all four datasets (1KP, Strassert2021, Vanderpool2020 and Pease2016)
# We estimated concatenated trees in IQ-Tree for the two shallow datasets (Vanderpool2020 and Pease2016) as the deep datasets were too large to reasonably run
# We estimated concatenated trees in RAxML-NG for the two deep datasets (1KP and Strassert2021)

# # To run this program: 
# # 1. Delete the lines below that include Caitlin's paths/variables
# # 2. Uncomment lines 25 to 43 inclusive and fill with your own variable names
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
# datasets_to_copy_loci_ASTRAL_IQTREE <-  c()
# datasets_to_copy_loci_RAxML <- c()
# datasets_to_estimate_ASTRAL_trees <- c()
# datasets_to_estimate_IQTREE_trees <- c()
# partition.by.codon.position <- FALSE
# use.modelfinder.models.for.partitions <- TRUE

### Caitlin's paths ###
run_location = "server"

if (run_location == "local"){
  # Datasets/dataset information
  input_names <- c( "1KP", "Strassert2021", "Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes_renamed/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths  in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("Astral/astral.5.7.5.jar","iqtree-2.0-rc1-MacOSX/bin/iqtree", "raxml-ng_v1.0.3_macos_x86_64/raxml-ng")
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("ASTRAL","IQTree", "RAxML-NG")
  
  # Select number of cores for parallelisation
  cores.to.use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci_ASTRAL_IQTREE <-  c()
  datasets_to_copy_loci_RAxML <- c("1KP", "Strassert2021")
  datasets_to_estimate_ASTRAL_trees <- c()
  datasets_to_estimate_IQTREE_trees <- c()
  datasets_to_estimate_RAxML_trees <- c("1KP", "Strassert2021")
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd - based on position in alignment file) 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use substitution models from the gene trees
  
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
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical_treelikeness repository is
  
  # Create a vector with all of the executable file paths in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/executables/ASTRAL/astral.5.7.5.jar",
                  "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree")
  names(exec_paths) <- c("ASTRAL","IQTree")
  
  # Select number of cores for parallelisation
  cores.to.use = 45
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci_ASTRAL_IQTREE <-  c()
  datasets_to_copy_loci_RAxML <- c("1KP", "Strassert2021")
  datasets_to_estimate_ASTRAL_trees <- c()
  datasets_to_estimate_IQTREE_trees <- c()
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use partition file with substitution models specified
}
### End of Caitlin's paths ###



##### Step 2: Open packages and source files for functions #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
library(phylotools)
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
all_datasets_to_copy <- unique(c(datasets_to_copy_loci_ASTRAL_IQTREE, datasets_to_copy_loci_RAxML))


##### Step 4: Assemble the dataframe of gene recombination results #####
if (length(all_datasets_to_copy) > 0){
  # Check whether a collated, trimmed recombination detection results file exists
  trimmed_gene_result_df_file <- paste0(csv_data_dir, "02_",paste(sort(all_datasets_to_copy), collapse="_"), "_collated_RecombinationDetection_TrimmedLoci.csv")
  pass_df_file <- paste0(csv_data_dir, "02_",paste(sort(all_datasets_to_copy), collapse="_"), "_RecombinationDetection_PassFail_record.csv")
  collated_exclude_file <- paste0(csv_data_dir, "01_IQ-Tree_warnings_",paste(sort(all_datasets_to_copy), collapse="_"), "_LociToExclude.csv")
  
  if (file.exists(trimmed_gene_result_df_file) & file.exists(pass_df_file)){
    gene_result_df <- read.csv(trimmed_gene_result_df_file, stringsAsFactors = TRUE)
    pass_df <- read.csv(pass_df_file, stringsAsFactors = TRUE)
  } else{
    # If the file doesn't exist, create it
    # Get a list of all the csv files in the csv_data_directory
    all_files <- list.files(csv_data_dir)
    # Get the results filenames for the datasets of interest
    all_results <- grep("RecombinationDetection_complete_collated_results", all_files, value = TRUE)
    all_results <- paste0(csv_data_dir, all_results)
    results <- c()
    for (dataset in datasets_to_copy_loci_ASTRAL_IQTREE){
      f <- grep(dataset, all_results, value = TRUE)
      results <- c(results, f)
    }
    # Remove duplicates
    results <- unique(results)
    # Open and attach the datasets
    gene_result_df <- as.data.frame(do.call(rbind, lapply(results, read.csv)))
    gene_result_df$match <- paste0(gene_result_df$dataset, ":", gene_result_df$loci_name)
    # If the collated total file hasn't been saved, save it
    all_gene_result_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci_ASTRAL_IQTREE), collapse="_"), "_collated_RecombinationDetection.csv")
    if (file.exists(all_gene_result_file) == FALSE){
      write.csv(gene_result_df, all_gene_result_file)
    }
    
    # Open the csv containing the list of loci to exclude from species tree analysis
    exclude_file <- paste0(csv_data_dir, grep("LociToExclude.csv", all_files, value = TRUE))
    # Open the csv files and bind into one dataframe
    exclude_df <- as.data.frame(do.call(rbind, lapply(exclude_file, read.csv)))
    # Collate and output if that file doesn't exist
    if (file.exists(collated_exclude_file) == FALSE){
      write.csv(exclude_df, file = collated_exclude_file, row.names = FALSE)
    }
    
    # Reduce down to the unique dataset/loci pairs to exclude from gene_result_df
    exclude_df <- exclude_df[,c("dataset", "loci")]
    exclude_df <- exclude_df[!duplicated(exclude_df),]
    exclude_df$match <- paste0(exclude_df$dataset, ":", exclude_df$loci)
    
    # Identify the row numbers of the loci to exclude
    exclude_loci <- which(gene_result_df$match %in% exclude_df$match)
    # Compare with the set of all loci in all datasets to get the loci to keep
    include_loci <- setdiff(1:nrow(gene_result_df), exclude_loci)
    
    # Trim the gene_result_df to remove the loci to be excluded
    gene_result_df <- gene_result_df[include_loci,]
    
    # Create new columns with pass/fail for each test
    gene_result_df <- make.pass.fail.column("pass_3seq", "X3SEQ_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_phi", "PHI_normal_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_phi_permute", "PHI_permutation_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_maxchi", "max_chi_squared_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_NSS", "NSS_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_geneconv_inner", "geneconv_inner_fragment_simulated_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_geneconv_outer", "geneconv_outer_fragment_simulated_p_value", gene_result_df)
    # Create a column for GeneConv where both inner and outer fragment p-value must be >0.05 to pass
    gene_result_df$pass_geneconv <- "FALSE"
    gene_result_df$pass_geneconv[((gene_result_df$geneconv_outer_fragment_simulated_p_value > 0.05) & 
                                    (gene_result_df$geneconv_inner_fragment_simulated_p_value > 0.05))] <- "TRUE"
    gene_result_df$pass_geneconv <- as.logical(gene_result_df$pass_geneconv)
    
    # Save the trimmed gene_result_df
    trimmed_gene_result_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci_ASTRAL_IQTREE), collapse="_"), "_collated_RecombinationDetection_TrimmedLoci.csv")
    write.csv(gene_result_df, file = trimmed_gene_result_df_file, row.names = FALSE)
    # Save a df of just the pass/fail info
    pass_df <- gene_result_df[,c("dataset", "loci_name", "alphabet", "n_taxa", "n_bp", "pass_3seq", "pass_phi", "pass_maxchi", 
                                 "pass_NSS", "pass_geneconv_inner", "pass_geneconv_outer", "pass_geneconv")]
    pass_df_file <- paste0(csv_data_dir, "02_",paste(sort(datasets_to_copy_loci_ASTRAL_IQTREE), collapse="_"), "_RecombinationDetection_PassFail_record.csv")
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
for (dataset in datasets_to_copy_loci_ASTRAL_IQTREE){
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
  
  ## filter gene_result_df by dataset
  dataset_df <- gene_result_df[gene_result_df$dataset == dataset,]
  
  ### Apply a single recombination detection test ###
  ## Iterate through each var and save the loci/trees for a tree made from only loci that pass the test
  # Make a list of the variables on which to filter the loci
  # vars values are columns from the gene_result_df
  vars <- c("pass_phi", "pass_maxchi", "pass_geneconv")
  # Assign output names for each of the variables in vars
  vars_names <- c("PHI", "maxchi", "geneconv")
  names(vars_names) <- vars
  # Iterate through each var: 
  for (v in vars){
    # Make a tree for all the loci that pass the test and all the loci that fail the test
    tree_type <- c("pass", "fail")
    for (tt in tree_type){
      print(paste0(dataset, " : ", vars_names[v], " : ", tt))
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
      if (partition.by.codon.position == TRUE){
        v_IQTree_name <- paste0(dataset,"_",v_name, "_", tt, "_IQTREE_partitioned")
      } else if (partition.by.codon.position == FALSE){
        v_IQTree_name <- paste0(dataset,"_",v_name, "_", tt, "_IQTREE") 
      }
      
      # Remove any loci that have an NA result for this test
      na_test_inds <- which(is.na(dataset_df[[v]]))
      dataset_inds <- 1:nrow(dataset_df)
      keep_inds <- setdiff(dataset_inds, na_test_inds)
      v_df <- dataset_df[keep_inds,]
      
      # Break up dataframe into only loci that pass the test
      v_inds <- which(dataset_df[,c(v)] == bool)
      # Use the indexes to subset the dataframe to just loci that pass the test 
      # (i.e. have a non significant p-value, meaning the null hypothesis of non-recombination/treelikeness cannot be rejected)
      v_df <- dataset_df[v_inds,]
      
      # Copy trees of all loci that pass the test into one file that can be fed into ASTRAL
      copy.loci.trees(v_df$loci_name, v_df$tree, category_output_folder, v_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
      # If running IQ-Tree analysis, copy all loci into a separate folder that can be fed into IQ-Tree
      # create the partition file required to run this IQ-Tree analysis
      partition.file.from.loci.list(loci_list = v_df$loci_name, directory = paste0(category_output_folder, v_IQTree_name, "/"),
                                    original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                    substitution_models = v_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
      
      # Create a record of which loci went into which analysis
      output_text <- v_df$loci_name
      write(output_text, file = v_text_name)
      # Add the number of loci in this category to the summary row
      summary_row <- c(summary_row, length(v_df$loci_name))
    }
  }
  
  
  
  ### Apply all three tests - get loci that pass all, and those that fail any ###
  # Make a tree for all the loci that pass all tests and all the loci that fail one or more test
  tree_type <- c("pass", "fail")
  for (tt in tree_type){
    print(paste0(dataset, " : all tests : ", tt))
    
    # Create output names
    all_text_name <- paste0(text_records_dir, dataset, "_allTests_",tt ,"_loci_record.txt")
    all_ASTRAL_name <- paste0(dataset,"_allTests_", tt, "_ASTRAL")
    if (partition.by.codon.position == TRUE){
      all_IQTree_name <- paste0(dataset,"_allTests_", tt, "_IQTREE_partitioned")
    } else if (partition.by.codon.position == FALSE){
      all_IQTree_name <- paste0(dataset,"_allTests_", tt, "_IQTREE") 
    }
    
    ## Subset the dataframe_df to loci that pass all three tests and those that fail one or more test
    allTest_df <- dataset_df
    # Select loci that pass/fail all three tests
    pass_inds <- which((allTest_df$pass_phi == TRUE) & (allTest_df$pass_maxchi == TRUE) & (allTest_df$pass_geneconv == TRUE))
    # setdiff(x,y) = elements in x but not in y
    fail_inds <- setdiff(1:nrow(allTest_df), pass_inds)
    # Subset the dataframe
    if (tt == "pass"){
      all_df <- allTest_df[pass_inds, ]
    } else if (tt == "fail"){
      all_df <- allTest_df[fail_inds, ] 
    }
    # Copy loci trees for ASTRAL
    copy.loci.trees(all_df$loci_name, all_df$tree, category_output_folder, all_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
    # Create the partition file required to run this IQ-Tree analysis
    partition.file.from.loci.list(loci_list = all_df$loci_name, directory = paste0(category_output_folder, all_IQTree_name, "/"),
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                  substitution_models = all_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
    # Create a record of which loci went into which analysis
    output_text <- all_df$loci_name
    write(output_text, file = all_text_name)
    # Add the number of loci in this category to the summary row
    summary_row <- c(summary_row, length(all_df$loci_name))
  }
  
  ### Apply no tests ###
  print(paste0(dataset, " : No tests -- include all loci"))
  NoTest_text_name <- paste0(text_records_dir, dataset, "_NoTest_loci_record.txt")
  NoTest_ASTRAL_name <- paste0(dataset,"_NoTest_ASTRAL")
  if (partition.by.codon.position == TRUE){
    NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE_partitioned")
  } else if (partition.by.codon.position == FALSE){
    NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE") 
  }
  # Copy loci trees for ASTRAL
  copy.loci.trees(dataset_df$loci_name, dataset_df$tree, category_output_folder, NoTest_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
  # Copy loci alignments for IQ-Tree
  # create the partition file required to run this IQ-Tree analysis
  partition.file.from.loci.list(loci_list = dataset_df$loci_name, directory = paste0(category_output_folder, NoTest_IQTree_name, "/"),
                                original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                substitution_models = dataset_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
  
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
  n_na_PHI <- length(which(is.na(dataset_df$PHI_normal_p_value)))
  n_na_maxchi <- length(which(is.na(dataset_df$max_chi_squared_p_value)))
  n_na_geneconv <- length(unique(c(which(is.na(dataset_df$geneconv_outer_fragment_simulated_p_value)),
                                   which(is.na(dataset_df$geneconv_inner_fragment_simulated_p_value))) ))
  summary_row <- c(summary_row, n_pass_PHI_maxchi, n_fail_PHI_maxchi, n_pass_PHI_geneconv, n_fail_PHI_geneconv, 
                   n_pass_maxchi_geneconv, n_fail_maxchi_geneconv, n_na_PHI, n_na_maxchi, n_na_geneconv)
  
  ### Write out the summary row as a dataframe ###
  names(summary_row) <- c("dataset", "n_pass_PHI", "n_fail_PHI", "n_pass_maxchi", "n_fail_maxchi", "n_pass_geneconv", "n_fail_geneconv",
                          "n_pass_allTests", "n_fail_allTests", "n_NoTest", "n_pass_PHI_maxchi", "n_fail_PHI_maxchi", "n_pass_PHI_geneconv", 
                          "n_fail_PHI_geneconv", "n_pass_maxchi_geneconv", "n_fail_maxchi_geneconv", "n_na_PHI", "n_na_maxchi", "n_na_geneconv")
  summary_df <- data.frame(as.list(summary_row))
  summary_op_file <- paste0(output_dirs[dataset], dataset, "_species_tree_summary.csv")
  write.csv(summary_df, file = summary_op_file, row.names = FALSE)
}

if (length(datasets_to_copy_loci_ASTRAL_IQTREE) > 0 | length(datasets_to_copy_loci_RAxML) > 0){
  # Collate species_tree_summary.csvs
  all_files <- list.files(output_dir, recursive = TRUE)
  summary_csvs <- grep("species_tree_summary.csv", all_files, value = TRUE)
  csv_list <- lapply(paste0(output_dir, summary_csvs), read.csv)
  csv_df <- do.call(rbind, csv_list)
  write.csv(csv_df, paste0(output_dir, "02_species_tree_summary_numbers.csv"))
}




##### Step 6: Estimate species trees in ASTRAL and IQ-Tree #####
### Estimate a species tree for each of the five categories
for (dataset in datasets_to_estimate_ASTRAL_trees){
  # Ensure the folder for species trees data exists
  category_output_folder <- paste0(output_dirs[dataset], "species_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  
  # Get list of all files in that folder
  all_category_folder_files <- list.files(category_output_folder)
  # Filter into ASTRAL text files and IQ-Tree folders
  astral_files <- paste0(grep("\\.txt", grep("ASTRAL", all_category_folder_files, value = TRUE), value = TRUE))
  # Contruct names of finished treefiles
  astral_files_finished_names <- paste0(category_output_folder, gsub(".txt", "_species.tre", astral_files))
  # Identify which ASTRAL analyses have not been run
  astral_files_to_run <- astral_files[!file.exists(astral_files_finished_names)]
  # Run remaining ASTRAL analyses
  if (length(astral_files_to_run) > 0){
    # Construct full file path
    astral_files_to_run <- paste0(category_output_folder, astral_files_to_run)
    # Estimate the species trees using ASTRAL
    mclapply(astral_files_to_run, ASTRAL.wrapper, exec_paths["ASTRAL"], mc.cores = cores.to.use)
  }
}

for (dataset in datasets_to_estimate_IQTREE_trees){
  # Ensure the folder for species trees data exists
  category_output_folder <- paste0(output_dirs[dataset], "species_trees/")
  if (dir.exists(category_output_folder) == FALSE){
    dir.create(category_output_folder)
  }
  
  # Get list of all files in that folder
  all_category_folder_files <- list.files(category_output_folder, recursive = TRUE, full.names = TRUE)
  all_category_folder_files <- gsub("//", "/", all_category_folder_files)
  iqtree_files <- grep("partitions.nex",grep(".nex.", grep("IQTREE", all_category_folder_files, value = TRUE), value = TRUE, invert = TRUE), value = TRUE)
  # Contruct names of finished treefiles
  iqtree_files_finished_names <- paste0(iqtree_files, ".contree")
  # Identify which IQTREE analyses have not been run
  partitions_to_run <- iqtree_files[!file.exists(iqtree_files_finished_names)]
  # Run remaining IQ-Tree analyses
  if (length(partitions_to_run) > 0){
    # Estimate the species tree on each folder of alignments using the partition files
    if (use.modelfinder.models.for.partitions == TRUE){
      mclapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"], 
               IQTREE_model_command = NA, mc.cores = cores.to.use)
    } else if (use.modelfinder.models.for.partitions == FALSE){
      mclapply(partitions_to_run, estimate.partitioned.IQTREE.species.tree, exec_paths["IQTree"],
               IQTREE_model_command = "MFP+MERGE", mc.cores = cores.to.use)
    }
  }
}



##### Step 7: Prepare partition and supermatrix files for tree estimation in RAxML #####
for (dataset in datasets_to_copy_loci_RAxML){
  # filter the gene_result_df for this dataset
  dataset_df <- gene_result_df[(gene_result_df$dataset == dataset), ]
  dataset_df$ModelFinder_model <- as.character(dataset_df$ModelFinder_model)
  dataset_df$loci_name <- as.character(dataset_df$loci_name)
  
  # identify file extension for alignment files
  if (dataset == "1KP"){
    file_extension = "fasta"
  } else if (dataset == "Strassert2021"){
    file_extension = "fas"
  }
  
  ## Create the supermatrix and partition file for the NoTest tree (tree estimated from all genes)
  # Create a new directory for this analysis
  raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_NoTest_RAxML/")
  if (dir.exists(raxml_dir) == FALSE){
    dir.create(raxml_dir)
  }
  # Identify those loci from the alignment folder
  all_als <- list.files(alignment_dir[dataset], recursive = TRUE)
  all_als <- grep(file_extension, all_als, value = TRUE)
  # Extend the alignment names to be the full alignment folder paths
  all_als_alignment_dir <- paste0(alignment_dir[dataset], all_als)
  # Get list of loci names for all alignments
  if (dataset == "1KP"){
    all_als_loci_names <- list.files(alignment_dir[dataset])
  } else {
    all_als_file_names <- unlist(strsplit(all_als, "\\."))
    all_als_file_name_parts <- which(all_als_file_names %in% c("filtered", "ginsi", "bmge", "merged", "fa", "divvy", "trimal", "fas"))
    all_als_loci_name_indexes <- setdiff(1:length(all_als_file_names), all_als_file_name_parts)
    all_als_loci_names <- all_als_file_names[all_als_loci_name_indexes]
  }
  # Filter out alignments that are not in the dataset_df$loci_name
  keep_als <- which(all_als_loci_names %in% dataset_df$loci_name)
  keep_al_paths <- all_als_alignment_dir[keep_als]
  # Create the output names for the supermatrix and partition files
  supermatrix_file <- paste0(raxml_dir, dataset, "_NoTest_supermat.phy")
  partition_file <- paste0(raxml_dir, dataset, "_NoTest_partition.txt")
  # Build PHYLIP supermatrix and RAxML partition file using aligned FASTA files
  supermat(keep_al_paths, outfile = supermatrix_file, partition.file = partition_file)
  # Reset the models in the partition file one at a time
  fix.all.models.in.partition.file(locus_names = dataset_df$loci_name, locus_models = dataset_df$ModelFinder_model, 
                                   dataset = dataset, partition_file = partition_file)
  
  ## Create the supermatrix and partition file for the PHI,pass and MaxChi,pass trees
  tests_to_run = c("PHI", "maxchi")
  for (test in tests_to_run){
    raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_", test, "_pass_RAxML/")
    if (dir.exists(raxml_dir) == FALSE){
      dir.create(raxml_dir)
    }
    
    # Set column names
    col_names <- c("pass_phi", "pass_maxchi")
    names(col_names) <- c("PHI", "maxchi")
    col <- col_names[test]
    
    # Remove any loci that have an NA result for this test
    na_test_inds <- which(is.na(dataset_df[[col]]))
    dataset_inds <- 1:nrow(dataset_df)
    keep_inds <- setdiff(dataset_inds, na_test_inds)
    test_df <- dataset_df[keep_inds,]
    
    # Break up dataframe into only loci that pass the test
    col_inds <- which(test_df[,c(col)] == TRUE)
    # Use the indexes to subset the dataframe to just loci that pass the test 
    # (i.e. have a non significant p-value, meaning the null hypothesis of non-recombination/treelikeness cannot be rejected)
    test_df <- test_df[col_inds,]
    
    # Identify file paths for those loci
    test_al_inds <- which(all_als_loci_names %in% test_df$loci_name)
    test_al_paths <- all_als_alignment_dir[test_al_inds]
    
    # Create the output names for the supermatrix and partition files
    supermatrix_file <- paste0(raxml_dir, dataset, "_", test, "_pass_supermat.phy")
    partition_file <- paste0(raxml_dir, dataset, "_", test, "_pass_partition.txt")
    # Build PHYLIP supermatrix and RAxML partition file using aligned FASTA files
    supermat(test_al_paths, outfile = supermatrix_file, partition.file = partition_file)
    # Reset the models in the partition file one at a time
    fix.all.models.in.partition.file(locus_names = test_df$loci_name, locus_models = test_df$ModelFinder_model, 
                                     dataset = dataset, partition_file = partition_file)
  }
}



##### Step 8: Estimate trees in RAxML #####
for (dataset in datasets_to_estimate_RAxML_trees){
  species_tree_dir <- paste0(output_dirs[dataset], "species_trees/")
  all_dirs <- list.dirs(species_tree_dir)
  all_dirs <- gsub("//", "/", all_dirs) # remove double slashes from file names
  raxml_dirs <- grep("RAxML", all_dirs, value = TRUE)
  
  tests_to_run <- c("NoTest", "PHI_pass", "maxchi_pass")
  for (test in tests_to_run){
    test_dir <- grep(paste0(test, "_RAxML"), raxml_dirs, value = TRUE)
    test_files <- list.files(test_dir)
    # Remove any RAxML files from the list of test files (by removing anything that has extra filename parts after the .phy or .txt)
    test_files <- grep(".phy.", test_files, invert = TRUE, value = TRUE)
    test_files <- grep(".txt.", test_files, invert = TRUE, value = TRUE)
    # Identify the partition and supermatrix files
    partition_file <- paste0(test_dir, "/", grep("partition.txt", test_files, value = TRUE))
    supermatrix_file <- paste0(test_dir, "/", grep("supermat.phy", test_files, value = TRUE))
    # Assemble RAxML-NG command line 
    raxml_call <- paste0("raxml-ng --all --msa ", supermatrix_file, " --model ", partition_file, 
                         " --prefix ", dataset, "_", test, " --brlen scaled --bs-metric fbp,tbe --bs-trees 100 --lh-epsilon 1")
    print(raxml_call)
    system(raxml_call)
  }
}



