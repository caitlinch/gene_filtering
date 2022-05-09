### gene_filtering/code/2_Species_Tree_Estimation.R
## R program to estimate trees from recombinant or non-recombinant loci
# Caitlin Cherryh 2021

## Additional software packages are required:
##     - ASTRAL (Zhang et. al. 2019) (https://github.com/smirarab/ASTRAL)
##     - IQTREE2 (Minh et. al. 2020) (http://www.iqtree.org/)
##     - RAxML-ng (Kozlov et. al. 2019) (https://github.com/amkozlov/raxml-ng)

##### Step 1: Set file paths and run variables #####
# input_names       <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir 
#                      (e.g. for 2 datasets, put same dataset first and same dataset last)
# alignment_dir     <- the folder(s) containing the alignments for each loci (should be vector with same length as input_names, with directories in same order as for input_names)
# csv_data_dir      <- directory containing the .csv file results from script 1_RecombinationDetection.R
# output_dir        <- where the coalescent/concatenated trees and tree comparisons will be stored 
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/gene_filtering)
# exec_paths        <- location to the software executables needed for analysis (ASTRAL and IQTREE)
# datasets_to_copy_loci_ASTRAL_IQTREE <- Out of the input names, select which datasets to copy loci trees for tree estimation in ASTRAL or IQ-Tree
# datasets_to_copy_loci_RAxML <- Out of the input names, select which datasets to copy loci trees for tree estimation in RAxML-NG
# datasets_to_estimate_ASTRAL_trees <- Out of the input names, select which datasets to estimate species trees in ASTRAL based on recombination results
# datasets_to_estimate_IQTREE_trees <- Out of the input names, select which datasets to estimate species trees in IQ-Tree based on recombination results
# datasets_to_estimate_RAxML_trees <- Out of the input names, select which datasets to estimate species trees in IQ-Tree based on recombination results
# dataset_tests_to_run <- list of tests to apply to each dataset (only apply categories that contain 50 or more loci)
# dataset_trees_to_estimate <- list of trees to estimate for each dataset (only estimate trees that contain 50 or more loci)
# partition.by.codon.position <- whether to run analysis partitioning by codon position
#                             <- set TRUE if you want to estimate species trees partitioning by codon position, and FALSE if you don't
#                             <- codon position simply counts every third base starting from 1st, 2nd, or 3rd base, and does not account for frame shift
# use.modelfinder.models.for.partitions <- can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. 
#                                       <- TRUE will use "-m MERGE" and include a charpartition with substitution models selected by ModelFinder in IQ-Tree
# use.free.rate.models.for.deep.datasets <- during our tree estimation step we found estimating trees from the deep dataset (1KP) was extremely
#                                           time consuming. This argument excludes the free rate parameters for the deep datasets by instead using the best
#                                           model found by ModelFinder that doesn't include a free rate parameter. 
#                                           To find the correct models needed to set this value to FALSE, run file 2.5_ExtractingModels_DeepDatasets.R

# We estimated species trees in ASTRAl for all four datasets (1KP, Whelan2017, Vanderpool2020 and Pease2016)
# We estimated concatenated trees in IQ-Tree for the two shallow datasets (Vanderpool2020 and Pease2016) as the deep datasets were too large to reasonably run
# We estimated concatenated trees in RAxML-NG for the 1000 Plants deep dataset (1KP)

# # To run this program: 
# # 1. Delete the lines below that include Caitlin's paths/variables
# # 2. Uncomment lines 25 to 43 inclusive and fill with your own variable names
# input_names <- c( "1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
# alignment_dir <- ""
# csv_data_dir <- ""
# output_dir <- ""
# maindir <- "/path/to/gene_filtering/"
# # Create a vector with all of the executable file paths:
# # exec_paths <- c("/path/to/ASTRAL_executable","/path/to/IQ-Tree_executable")
# # names(exec_paths) <- c("ASTRAL","IQTree")
# # To access a path: exec_paths[["name"]]
# exec_paths <- c()
# names(exec_paths) <- c("ASTRAL","IQTree")
# datasets_to_copy_loci_ASTRAL_IQTREE <-  c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
# datasets_to_copy_loci_RAxML <- c("1KP")
# datasets_to_estimate_ASTRAL_trees <- c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
# datasets_to_estimate_IQTREE_trees <- c("Whelan2017", "Vanderpool2020", "Pease2016")
# datasets_to_estimate_RAxML_trees <- c("1KP")
# dataset_tests_to_run <- list("1KP" = c("PHI", "maxchi"),
#                              "Whelan2017" = c("PHI", "maxchi", "geneconv"),
#                              "Vanderpool2020" = c("PHI", "maxchi", "geneconv", "all"),
#                              "Pease2016" = c("PHI", "maxchi", "geneconv", "all"))
# # List of trees to estimate for each dataset
# dataset_trees_to_estimate <- list("1KP" = c("PHI,pass", "maxchi,pass"),
#                                   "Whelan2017" = c("PHI", "maxchi", "geneconv,pass", "geneconv,fail"),
#                                   "Vanderpool2020" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
#                                                        "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"),
#                                   "Pease2016" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
#                                                   "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"))
# partition.by.codon.position <- FALSE
# use.modelfinder.models.for.partitions <- TRUE
# use.free.rate.models.for.deep.datasets <- TRUE

### Caitlin's paths ###
run_location = "local"

if (run_location == "local"){
  # Datasets/dataset information
  input_names <- c( "1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes_renamed/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2017/genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/" # where the empirical treelikeness code is
  
  # Create a vector with all of the executable file paths  in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("Astral/astral.5.7.5.jar","iqtree-2.0-rc1-MacOSX/bin/iqtree", "raxml-ng_v1.0.3_macos_x86_64/raxml-ng")
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("ASTRAL","IQTree", "RAxML-NG")
  
  # Select number of cores for parallelisation
  cores.to.use = 1
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci_ASTRAL_IQTREE <-  c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_copy_loci_RAxML <- c("1KP")
  datasets_to_estimate_ASTRAL_trees <- c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_estimate_IQTREE_trees <- c("Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_estimate_RAxML_trees <- c("1KP")
  # List of tests to run for each dataset
  dataset_tests_to_run <- list("1KP" = c("PHI", "maxchi"),
                               "Whelan2017" = c("PHI", "maxchi", "geneconv"),
                               "Vanderpool2020" = c("PHI", "maxchi", "geneconv", "all"),
                               "Pease2016" = c("PHI", "maxchi", "geneconv", "all"))
  # List of trees to estimate for each dataset
  dataset_trees_to_estimate <- list("1KP" = c("PHI,pass", "maxchi,pass"),
                                    "Whelan2017" = c("PHI,pass", "maxchi,pass", "geneconv,pass"),
                                    "Vanderpool2020" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
                                                         "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"),
                                    "Pease2016" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
                                                    "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"))
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd - based on position in alignment file) 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use substitution models from the gene trees
  use.free.rate.models.for.deep.datasets = FALSE # whether to use modelFinder best model or best model that doesn't include a free rates parameter
  
} else if (run_location=="server"){
  # Datasets/dataset information
  input_names <- c( "1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/data/caitlin/empirical_treelikeness/Data_1KP/",
                     "/data/caitlin/empirical_treelikeness/Data_Whelan2017/",
                     "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/",
                     "/data/caitlin/empirical_treelikeness/Data_Pease2016/")
  csv_data_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical_treelikeness repository is
  
  # Create a vector with all of the executable file paths in this order: ASTRAL, IQ-Tree
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/executables/ASTRAL/astral.5.7.5.jar",
                  "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree",
                  "/data/caitlin/linux_executables/raxml-ng")
  names(exec_paths) <- c("ASTRAL","IQTree", "RAxML-NG")
  
  # Select number of cores for parallelisation
  cores.to.use = 45
  
  # Select datasets to run analysis and collect results
  datasets_to_copy_loci_ASTRAL_IQTREE <-  c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_copy_loci_RAxML <- c("1KP")
  datasets_to_estimate_ASTRAL_trees <- c("1KP", "Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_estimate_IQTREE_trees <- c("Whelan2017", "Vanderpool2020", "Pease2016")
  datasets_to_estimate_RAxML_trees <- c("1KP")
  # List of tests to run for each dataset
  dataset_tests_to_run <- list("1KP" = c("PHI", "maxchi"),
                               "Whelan2017" = c("PHI", "maxchi", "geneconv"),
                               "Vanderpool2020" = c("PHI", "maxchi", "geneconv", "all"),
                               "Pease2016" = c("PHI", "maxchi", "geneconv", "all"))
  # List of trees to estimate for each dataset
  dataset_trees_to_estimate <- list("1KP" = c("PHI,pass", "maxchi,pass"),
                                    "Whelan2017" = c("PHI,pass", "maxchi,pass", "geneconv,pass"),
                                    "Vanderpool2020" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
                                                         "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"),
                                    "Pease2016" = c("PHI,pass", "PHI,fail", "maxchi,pass", "maxchi,fail",
                                                    "geneconv,pass", "geneconv,fail", "all,pass", "all,fail"))
  partition.by.codon.position = FALSE # can be TRUE or FALSE: TRUE will partition by codon position (1st, 2nd and 3rd), FALSE will treat each gene homogeneously 
  use.modelfinder.models.for.partitions = TRUE # can be TRUE or FALSE. FALSE will use "-m MFP+MERGE" in IQ-Tree. TRUE will use partition file with substitution models specified
  use.free.rate.models.for.deep.datasets = FALSE # whether to use modelFinder best model or best model that doesn't include a free rates parameter
}
### End of Caitlin's paths ###



##### Step 2: Open packages and source files for functions #####
print("opening packages")
library(ape) # analyses of phylogenetics and evolution
library(parallel) # support for parallel computation
library(phangorn) # phylogenetic reconstruction and analysis
library(phytools) # tools for comparative biology and phylogenetics
if (length(datasets_to_copy_loci_RAxML) > 0){
  library(phylotools)
}
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



##### Step 4: Assemble the dataframe of gene recombination results #####
# This section:
#   1. Takes in the recombination detection results from file 1. 
#   2. Checks the LociToExclude file, and removes any loci that have been flagged (by IQ-Tree warnings)
#   3. Outputs a csv of the trimmed data
# All downstream analyses use the trimmed dataframe (i.e. the complete results, minus individual loci that were removed on the basis of IQ-Tree warnings)

if (length(input_names) > 0){
  # Check whether a collated, trimmed recombination detection results file exists
  trimmed_gene_result_df_file <- paste0(csv_data_dir, "02_AllDatasets_collated_RecombinationDetection_TrimmedLoci.csv")
  pass_df_file <- paste0(csv_data_dir, "02_AllDatasets_RecombinationDetection_PassFail_record.csv")
  collated_exclude_file <- paste0(csv_data_dir, "01_AllDatasets_IQ-Tree_warnings_LociToExclude.csv")
  
  if (file.exists(trimmed_gene_result_df_file) & file.exists(pass_df_file)){
    # If collated RecombinationDetection results file for all datasets has been created, open that file
    gene_result_df <- read.csv(trimmed_gene_result_df_file, stringsAsFactors = TRUE)
    pass_df <- read.csv(pass_df_file, stringsAsFactors = TRUE)
  } else{
    # If the collated, trimmed recombination detection results file doesn't exist, create it

    all_gene_result_file <- paste0(csv_data_dir, "01_AllDatasets_RecombinationDetection_complete_collated_results.csv")
    if (file.exists(all_gene_result_file)){
      gene_result_df <- read.csv(all_gene_result_file)
      gene_result_df$match <- paste0(gene_result_df$dataset, ":", gene_result_df$loci_name)
    } else {
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
      all_gene_result_file <- paste0(csv_data_dir, "01_AllDatasets_RecombinationDetection_complete_collated_results.csv")
      if (file.exists(all_gene_result_file) == FALSE){
        write.csv(gene_result_df, all_gene_result_file)
      }
    }
    
    # If collated LociToExclude file for all datasets has been collated, open that file
    # If not, create it.
    if (file.exists(collated_exclude_file) == TRUE){
      exclude_df <- read.csv(collated_exclude_file)
    } else {
      # Get a list of all the csv files in the csv_data_directory
      all_files <- list.files(csv_data_dir)
      # Open the csv containing the list of loci to exclude from species tree analysis
      exclude_file <- paste0(csv_data_dir, grep("LociToExclude.csv", all_files, value = TRUE))
      # Open the csv files and bind into one dataframe
      exclude_df <- as.data.frame(do.call(rbind, lapply(exclude_file, read.csv)))
      # Collate and output if that file doesn't exist
      if (file.exists(collated_exclude_file) == FALSE){
        write.csv(exclude_df, file = collated_exclude_file, row.names = FALSE)
      }
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
    gene_result_df <- make.pass.fail.column("pass_phi", "PHI_normal_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_phi_permute", "PHI_permutation_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_maxchi", "max_chi_squared_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_geneconv_inner", "geneconv_inner_fragment_simulated_p_value", gene_result_df)
    gene_result_df <- make.pass.fail.column("pass_geneconv_outer", "geneconv_outer_fragment_simulated_p_value", gene_result_df)
    # Create a column for GeneConv where both inner and outer fragment p-value must be >0.05 to pass
    gene_result_df$pass_geneconv <- "0"
    gene_result_df$pass_geneconv[((gene_result_df$geneconv_outer_fragment_simulated_p_value > 0.05) & 
                                    (gene_result_df$geneconv_inner_fragment_simulated_p_value > 0.05))] <- "TRUE"
    gene_result_df$pass_geneconv[((gene_result_df$geneconv_outer_fragment_simulated_p_value <= 0.05) | 
                                    (gene_result_df$geneconv_inner_fragment_simulated_p_value <= 0.05))] <- "FALSE"
    gene_result_df$pass_geneconv[(is.na(gene_result_df$geneconv_outer_fragment_simulated_p_value) |
                                    is.na(gene_result_df$geneconv_inner_fragment_simulated_p_value))] <- "NA"
    gene_result_df$pass_geneconv <- as.logical(gene_result_df$pass_geneconv)
    # Select columns to save for gene results df
    gene_result_df <- gene_result_df[,c("dataset", "loci_name", "alphabet", "best_model", "ModelFinder_model", "n_taxa", "n_bp",
                                        "analytical_PHI_mean_value",  "permutation_PHI_mean_value", "analytical_PHI_variance_value",
                                        "permutation_PHI_variance_value", "analytical_PHI_observed_value", "permutation_PHI_observed_value",
                                        "NSS_p_value", "max_chi_squared_p_value", "PHI_permutation_p_value", "PHI_normal_p_value",
                                        "geneconv_seed", "geneconv_num_global_inner_fragments", "geneconv_num_global_outer.sequence_fragments",
                                        "geneconv_num_pairwise_inner_fragments", "geneconv_num_pairwise_outer.sequence_fragments",
                                        "geneconv_inner_fragment_maximum_blast.like_score", "geneconv_inner_fragment_simulated_p_value",
                                        "geneconv_inner_fragment_sd_above_sim_mean", "geneconv_inner_fragment_sd_of_sim", 
                                        "geneconv_outer_fragment_maximum_blast.like_score", "geneconv_outer_fragment_simulated_p_value",
                                        "geneconv_outer_fragment_sd_above_sim_mean", "geneconv_outer_fragment_sd_of_sim", "tree",
                                        "pass_phi", "pass_phi_permute", "pass_maxchi", "pass_geneconv_inner", "pass_geneconv_outer",
                                        "pass_geneconv")]
    
    # Save the trimmed gene_result_df and pass_df_file
    write.csv(gene_result_df, file = trimmed_gene_result_df_file, row.names = FALSE)
    # Save a df of just the pass/fail info
    pass_df <- gene_result_df[,c("dataset", "loci_name", "alphabet", "n_taxa", "n_bp", "pass_phi", "pass_maxchi", 
                                 "pass_NSS", "pass_geneconv_inner", "pass_geneconv_outer", "pass_geneconv")]
    write.csv(pass_df, file = pass_df_file, row.names = FALSE)
  }
}



##### Step 5: Categorize loci by test results and prepare subsets of data for tree estimation #####
# Estimate trees from the putatively non-recombinant and the putatively recombinant loci for each dataset 
# Iterate through each dataset and:
#    For each recombination detection test:
#        1. Record which loci pass the test
#        2. Estimate a species tree from the loci that pass the test
#        3. Record which loci fail the test
#        4. Estimate a species tree from the loci that fail the test
#    Then:
#        1. Estimate a species tree from all loci
#        2. Estimate a species tree from loci that pass every recombination detection test

# This section of the code prepares all files for IQ-Tree and ASTRAL to estimate the trees described above.
# Trees are estimated in the next section (Section 6)

### Save the loci trees (for ASTRAL) and the loci alignment (for IQ-Tree)
for (dataset in datasets_to_copy_loci_ASTRAL_IQTREE){
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
  # Make a list of the three tests for recombination (used to index columns in gene_result_df)
  vars <- c("PHI" = "pass_phi", "maxchi" = "pass_maxchi", "geneconv" = "pass_geneconv")
  
  # Get list of vars to run and trees to estimate for this dataset
  dataset_vars <- vars[which(names(vars) %in% dataset_tests_to_run[[dataset]])]
  dataset_var_trees <- dataset_trees_to_estimate[[dataset]]
  
  # Iterate through each var: 
  for (v in dataset_vars){
    # Make a tree for all the loci that pass the test and all the loci that fail the test
    tree_type <- c("pass", "fail")
    for (tt in tree_type){
      print(paste0(dataset, " : ", names(dataset_vars)[which(dataset_vars == v)], " : ", tt))
      
      # Check whether tree is in list to estimate 
      v_tt <- paste0(names(v),",",tt)
      run_check <- "PHI,pass" %in% dataset_var_trees
      
      # If run_check = TRUE, this tree is in the list to estimate.
      # Collect the loci to estimate this tree 
      if (run_check == TRUE){
        # Set boolean to collect either loci that passed or failed the test
        if (tt == "pass"){
          bool = TRUE
        } else if (tt == "fail"){
          bool = FALSE
        }
        
        # Get short version of name for output files
        v_name <- names(dataset_vars)[which(dataset_vars == v)]
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
        
        # If there are more than 50 loci (50+ rows in v_df), prepare the files to run this analysis
        if (nrow(v_df) >= 50){
          # Copy trees of all loci that pass the test into one file that can be fed into ASTRAL
          copy.loci.trees(v_df$loci_name, v_df$tree, category_output_folder, v_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
          # If running IQ-Tree analysis, copy all loci into a separate folder that can be fed into IQ-Tree
          # create the partition file required to run this IQ-Tree analysis
          partition.file.from.loci.list(loci_list = v_df$loci_name, directory = paste0(category_output_folder, v_IQTree_name, "/"),
                                        original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                        substitution_models = v_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
        } else {
          print(paste0("ERROR: less than fifty loci in this category."))
          print(paste0("CATEGORY: Dataset = " ,dataset, ". Test = ", names(dataset_vars)[which(dataset_vars == v)], ". Pass/fail = ", tt))
          print("CATEGORY NOT PREPARED FOR TREE ESTIMATION")
        }
        
        # Create a record of which loci went into which analysis
        output_text <- v_df$loci_name
        write(output_text, file = v_text_name)
      }
    }
  }
  
  ### Apply all three tests - get loci that pass all, and those that fail any ###
  # Make a tree for all the loci that pass all tests and all the loci that fail one or more test
  if ("all,fail" %in% dataset_var_trees | "all,pass" %in% dataset_var_trees){
    # Collect the trees from the "all" test to estimate
    all_vars <- grep("all", dataset_var_trees, value = TRUE)
    # Iterate through each of the "all" tests in the dataset_var_trees vector (may be "all,pass", "all,fail", or both)
    for (all_var in all_vars){
      # Assign test type (tt) as "pass" or "fail" to match the current all_var
      if (all_var == "all,pass"){
        tt = "pass"
      } else if (all_var == "all,fail"){
        tt = fail
      }
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
      
      # If there are more than 50 loci (50+ rows in all_df), prepare the files to run this analysis
      if (nrow(all_df) >= 50){
        # Copy loci trees for ASTRAL
        copy.loci.trees(all_df$loci_name, all_df$tree, category_output_folder, all_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
        # Create the partition file required to run this IQ-Tree analysis
        partition.file.from.loci.list(loci_list = all_df$loci_name, directory = paste0(category_output_folder, all_IQTree_name, "/"),
                                      original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                      substitution_models = all_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
      }
      
      # Create a record of which loci went into which analysis
      output_text <- all_df$loci_name
      write(output_text, file = all_text_name)
    }
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
  
  # If there are more than 50 loci (50+ rows in dataset_df), prepare the files to run this analysis
  if (nrow(dataset_df) >= 50){
    # Copy loci trees for ASTRAL
    copy.loci.trees(dataset_df$loci_name, dataset_df$tree, category_output_folder, NoTest_ASTRAL_name, copy.all.individually = FALSE, copy.and.collate = TRUE)
    # Copy loci alignments for IQ-Tree
    # create the partition file required to run this IQ-Tree analysis
    partition.file.from.loci.list(loci_list = dataset_df$loci_name, directory = paste0(category_output_folder, NoTest_IQTree_name, "/"),
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                  substitution_models = dataset_df$ModelFinder_model, add.codon.positions = partition.by.codon.position)
  }
  
  ### Create a record of which loci went into which analysis ###
  output_text <- dataset_df$loci_name
  write(output_text, file = NoTest_text_name)
  
  # Expand the summary row to identify the number of loci that pass or fail the same tests
  n_total <- nrow(dataset_df)
  n_pass_PHI <- nrow(dataset_df[(dataset_df$pass_phi == TRUE),])
  n_fail_PHI <- nrow(dataset_df[(dataset_df$pass_phi == FALSE),])
  n_na_PHI <- length(which(is.na(dataset_df$PHI_normal_p_value)))
  n_pass_maxchi <- nrow(dataset_df[(dataset_df$pass_maxchi == TRUE),])
  n_fail_maxchi <- nrow(dataset_df[(dataset_df$pass_maxchi == FALSE),])
  n_na_maxchi <- length(which(is.na(dataset_df$max_chi_squared_p_value)))
  n_pass_geneconv <- nrow(dataset_df[(dataset_df$pass_geneconv == TRUE),])
  n_fail_geneconv <- nrow(dataset_df[(dataset_df$pass_geneconv == FALSE),])
  n_na_geneconv <- nrow(dataset_df[is.na(dataset_df$pass_geneconv),])
  n_pass_all <- length(which((dataset_df$pass_phi == TRUE) & (dataset_df$pass_maxchi == TRUE) & (dataset_df$pass_geneconv == TRUE)))
  n_fail_all <- length(setdiff(1:nrow(dataset_df), which((dataset_df$pass_phi == TRUE) & (dataset_df$pass_maxchi == TRUE) & (dataset_df$pass_geneconv == TRUE))))
  n_pass_PHI_maxchi <- nrow(dataset_df[((dataset_df$pass_phi == TRUE) & (dataset_df$pass_maxchi == TRUE)), ])
  n_fail_PHI_maxchi <- nrow(dataset_df[((dataset_df$pass_phi == FALSE) & (dataset_df$pass_maxchi == FALSE)), ])
  n_pass_PHI_geneconv <- nrow(dataset_df[((dataset_df$pass_phi == TRUE) & (dataset_df$pass_geneconv == TRUE)), ])
  n_fail_PHI_geneconv <- nrow(dataset_df[((dataset_df$pass_phi == FALSE) & (dataset_df$pass_geneconv == FALSE)), ])
  n_pass_maxchi_geneconv <- nrow(dataset_df[((dataset_df$pass_maxchi == TRUE) & (dataset_df$pass_geneconv == TRUE)), ])
  n_fail_maxchi_geneconv <- nrow(dataset_df[((dataset_df$pass_maxchi == FALSE) & (dataset_df$pass_geneconv == FALSE)), ])
  summary_row <- c(dataset, n_total, n_pass_PHI, n_fail_PHI, n_na_PHI, n_pass_maxchi, n_fail_maxchi, n_na_maxchi, n_pass_geneconv, n_fail_geneconv, n_na_geneconv,
                   n_pass_all, n_fail_all, n_pass_PHI_maxchi, n_fail_PHI_maxchi, n_pass_PHI_geneconv, n_fail_PHI_geneconv,
                   n_pass_maxchi_geneconv, n_fail_maxchi_geneconv)
  
  ### Write out the summary row as a dataframe ###
  names(summary_row) <- c("dataset", "n_total", "n_pass_PHI", "n_fail_PHI", "n_na_PHI", "n_pass_maxchi", "n_fail_maxchi", "n_na_maxchi", "n_pass_geneconv", "n_fail_geneconv", "n_na_geneconv",
                          "n_pass_all", "n_fail_all", "n_pass_PHI_maxchi", "n_fail_PHI_maxchi", "n_pass_PHI_geneconv", "n_fail_PHI_geneconv",
                          "n_pass_maxchi_geneconv", "n_fail_maxchi_geneconv")
  summary_df <- data.frame(as.list(summary_row))
  summary_op_file <- paste0(output_dirs[dataset], dataset, "_species_tree_summary.csv")
  write.csv(summary_df, file = summary_op_file, row.names = FALSE)
}

if (length(datasets_to_copy_loci_ASTRAL_IQTREE) > 0 | length(datasets_to_copy_loci_RAxML) > 0){
  # Collate species_tree_summary.csvs
  all_files <- list.files(output_dir, recursive = TRUE)
  summary_csvs <- grep("species_tree_summary.csv", all_files, value = TRUE)
  summary_csvs <- summary_csvs[grepl(paste0(input_names, collapse = "|"), summary_csvs)]
  csv_list <- lapply(paste0(output_dir, summary_csvs), read.csv)
  csv_df <- do.call(rbind, csv_list)
  write.csv(csv_df, paste0(output_dir, "02_species_tree_summary_numbers.csv"))
}

# If estimating deep datasets in IQ-Tree and using models without free rate parameters, create partition file for IQ-Tree for those analyses
if ((use.free.rate.models.for.deep.datasets == FALSE) & 
    ("1KP" %in% datasets_to_copy_loci_ASTRAL_IQTREE)){
  for (dataset in c("1KP")){
    # Create new folders to put these tree files/loci files and records in
    category_output_folder <- paste0(output_dirs[dataset], "species_trees/")
    if (dir.exists(category_output_folder) == FALSE){
      dir.create(category_output_folder)
    }
    
    # Filter gene_result_df by dataset
    dataset_df <- gene_result_df[gene_result_df$dataset == dataset,]
    
    # Open file containing best models without free rate parameters
    all_output_files <- list.files(csv_data_dir)
    noFreeRate_csvs <- grep("loci_models_noFreeRates", all_output_files, value = TRUE)
    noFreeRate_dataset_csv <- grep(dataset, noFreeRate_csvs, value = TRUE)
    noFreeRate_df <- read.csv(paste0(csv_data_dir, noFreeRate_dataset_csv), stringsAsFactors = FALSE)
    
    # Trim noFreeRate_df to only include rows that appear in the dataset_df$loci_name column
    keep_row_inds <- which(noFreeRate_df$loci %in% dataset_df$loci_name)
    noFreeRate_df <- noFreeRate_df[keep_row_inds, ]
    noFreeRate_df$loci <- as.character(noFreeRate_df$loci)
    # Order noFreeRate_df in same order as dataset_df
    noFreeRate_df <- noFreeRate_df[match(dataset_df$loci_name, noFreeRate_df$loci),]
    
    ### Create partition file to estimate tree for full set of loci: NoTests ###
    print(paste0(dataset, " : No tests -- include all loci"))
    if (partition.by.codon.position == TRUE){
      NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE_noFreeRates_partitioned")
    } else if (partition.by.codon.position == FALSE){
      NoTest_IQTree_name <- paste0(dataset,"_NoTest_IQTREE_noFreeRates") 
    }
    iqtree_dir <- paste0(category_output_folder, NoTest_IQTree_name, "/")
    if (dir.exists(iqtree_dir) == FALSE){
      dir.create(iqtree_dir)
    }
    
    # Create the partition file required to run this IQ-Tree analysis
    partition.file.from.loci.list(loci_list = noFreeRate_df$loci, directory = iqtree_dir,
                                  original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                  substitution_models = noFreeRate_df$no_freerates_best_model, add.codon.positions = partition.by.codon.position)
    
    ### Create partition file for PHI and MaxChi tree estimation runs ### 
    tests_to_run = c("PHI", "maxchi")
    for (test in tests_to_run){
      print(paste0("Dataset: ", dataset, " : ", test))
      # Create directory
      if (partition.by.codon.position == TRUE){
        test_IQTree_name <- paste0(dataset,"_", test, "_pass_IQTREE_noFreeRates_partitioned")
      } else if (partition.by.codon.position == FALSE){
        test_IQTree_name <- paste0(dataset,"_", test, "_pass_IQTREE_noFreeRates") 
      }
      iqtree_dir <- paste0(category_output_folder, test_IQTree_name, "/")
      if (dir.exists(iqtree_dir) == FALSE){
        dir.create(iqtree_dir)
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
      
      # Trim noFreeRate_df to only include rows that appear in the dataset_df$loci_name column
      keep_row_inds <- which(noFreeRate_df$loci %in% test_df$loci_name)
      test_noFreeRate_df <- noFreeRate_df[keep_row_inds, ]
      test_noFreeRate_df$loci <- as.character(test_noFreeRate_df$loci)
      # Order noFreeRate_df in same order as dataset_df
      test_noFreeRate_df <- test_noFreeRate_df[match(test_df$loci_name, test_noFreeRate_df$loci),]
      
      # Create the partition file required to run this IQ-Tree analysis
      partition.file.from.loci.list(loci_list = test_noFreeRate_df$loci, directory = iqtree_dir,
                                    original_alignment_folder = alignment_dir[[dataset]], add.charpartition.models = TRUE,
                                    substitution_models = test_noFreeRate_df$no_freerates_best_model, add.codon.positions = partition.by.codon.position)
      
    }
  }
}




##### Step 6: Estimate species trees in ASTRAL and IQ-Tree #####
### Estimate a species trees from the files prepared in Section 5.
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
# For the 1KP dataset, estimating a ML tree in IQ-Tree was not computationally feasible
# We estimated the tree in RAxML-NG instead

for (dataset in datasets_to_copy_loci_RAxML){
  print(paste0("Dataset: ", dataset))
  # filter the gene_result_df for this dataset
  dataset_df <- gene_result_df[(gene_result_df$dataset == dataset), ]
  dataset_df$ModelFinder_model <- as.character(dataset_df$ModelFinder_model)
  dataset_df$loci_name <- as.character(dataset_df$loci_name)
  
  # identify file extension for alignment files
  if (dataset == "1KP"){
    file_extension = "fasta"
  }
  
  ## Create the supermatrix and partition file for the NoTest tree (tree estimated from all genes)
  print(paste0("Dataset: ", dataset, " : No test"))
  # Create a new directory for this analysis
  if (use.free.rate.models.for.deep.datasets == TRUE){
    raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_NoTest_RAxML/")
  } else if(use.free.rate.models.for.deep.datasets == FALSE){
    raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_NoTest_RAxML_noFreeRates/")
  }
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
  if (use.free.rate.models.for.deep.datasets == TRUE){
    # Reset the models in the partition file one at a time
    fix.all.models.in.partition.file(locus_names = dataset_df$loci_name, locus_models = dataset_df$ModelFinder_model, 
                                     dataset = dataset, partition_file = partition_file)
  } else if (use.free.rate.models.for.deep.datasets == FALSE){
    # Open file containing best models without free rate parameters
    all_output_files <- list.files(csv_data_dir)
    noFreeRate_csvs <- grep("loci_models_noFreeRates", all_output_files, value = TRUE)
    noFreeRate_dataset_csv <- grep(dataset, noFreeRate_csvs, value = TRUE)
    noFreeRate_df <- read.csv(paste0(csv_data_dir, noFreeRate_dataset_csv), stringsAsFactors = FALSE)
    # Trim noFreeRate_df to only include rows that appear in the dataset_df$loci_name column
    keep_row_inds <- which(noFreeRate_df$loci %in% dataset_df$loci_name)
    noFreeRate_df <- noFreeRate_df[keep_row_inds, ]
    noFreeRate_df$loci <- as.character(noFreeRate_df$loci)
    # Order noFreeRate_df in same order as dataset_df
    noFreeRate_df <- noFreeRate_df[match(dataset_df$loci_name, noFreeRate_df$loci),]
    # Extract the names and models of the loci and feed them into the function fix.all.models.in.partition.file
    fix.all.models.in.partition.file(locus_names = noFreeRate_df$loci, locus_models = noFreeRate_df$no_freerates_best_model, 
                                     dataset = dataset, partition_file = partition_file)
  }
  
  ## Create the supermatrix and partition file for the PHI,pass and MaxChi,pass trees
  tests_to_run = c("PHI", "maxchi")
  for (test in tests_to_run){
    print(paste0("Dataset: ", dataset, " : ", test))
    # Create directory
    if (use.free.rate.models.for.deep.datasets == TRUE){
      raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_", test, "_pass_RAxML/")
    } else if(use.free.rate.models.for.deep.datasets == FALSE){
      raxml_dir <- paste0(output_dirs[dataset], "species_trees/", dataset, "_", test, "_pass_RAxML_noFreeRates/")
    }
    
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
    
    if (use.free.rate.models.for.deep.datasets == TRUE){
      # Reset the models in the partition file one at a time
      fix.all.models.in.partition.file(locus_names = test_df$loci_name, locus_models = test_df$ModelFinder_model, 
                                       dataset = dataset, partition_file = partition_file)
    } else if (use.free.rate.models.for.deep.datasets == FALSE){
      # Open file containing best models without free rate parameters
      all_output_files <- list.files(csv_data_dir)
      noFreeRate_csvs <- grep("loci_models_noFreeRates", all_output_files, value = TRUE)
      noFreeRate_dataset_csv <- grep(dataset, noFreeRate_csvs, value = TRUE)
      noFreeRate_df <- read.csv(paste0(csv_data_dir, noFreeRate_dataset_csv), stringsAsFactors = FALSE)
      # Trim noFreeRate_df to only include rows that appear in the test_df$loci_name column
      keep_row_inds <- which(noFreeRate_df$loci %in% test_df$loci_name)
      noFreeRate_df <- noFreeRate_df[keep_row_inds, ]
      noFreeRate_df$loci <- as.character(noFreeRate_df$loci)
      # Order noFreeRate_df in same order as test_df
      noFreeRate_df <- noFreeRate_df[match(test_df$loci_name, noFreeRate_df$loci),]
      # Extract the names and models of the loci and feed them into the function fix.all.models.in.partition.file
      fix.all.models.in.partition.file(locus_names = noFreeRate_df$loci, locus_models = noFreeRate_df$no_freerates_best_model, 
                                       dataset = dataset, partition_file = partition_file)
    }
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
    # Set working directory to test_dir to save output files with partition/supermatrix file
    setwd(test_dir)
    # Remove any RAxML files from the list of test files (by removing anything that has extra filename parts after the .phy or .txt)
    test_files <- grep(".phy.", test_files, invert = TRUE, value = TRUE)
    test_files <- grep(".txt.", test_files, invert = TRUE, value = TRUE)
    # Identify the partition and supermatrix files
    partition_file <- paste0(test_dir, "/", grep("partition.txt", test_files, value = TRUE))
    supermatrix_file <- paste0(test_dir, "/", grep("supermat.phy", test_files, value = TRUE))
    # Assemble RAxML-NG command line (to estimate tree with ML and perform bootstraps)
    # raxml_call <- paste0(exec_paths[["RAxML-NG"]], " --all --msa ", supermatrix_file, " --model ", partition_file, 
    #                      " --prefix ", dataset, "_", test, " --brlen scaled --bs-metric fbp,tbe --bs-trees 100 --lh-epsilon 1")
    # Assemble RAxML-NG command line (to estimate tree with ML)
    # Perform quick seach from single random starting tree using "--search1"
    # Perform search using n parsimony starting trees by adding "--tree pars{n}" (can also have "pars{n},rand{m}" to specify n parsimony and m random starting trees)
    raxml_call <- paste0(exec_paths[["RAxML-NG"]], " --search --msa ", supermatrix_file, " --model ", partition_file, 
                         " --prefix ", dataset, "_", test, " --brlen scaled --tree pars{1} --threads ", cores.to.use, " --lh-epsilon 1")
    print(raxml_call)
    system(raxml_call)
  }
}



