### empirical_treelikeness/code/3_Species_Tree_Comparison.R
## R program to identify the best fitting tree for a candidate dataset
## Additional software/packages required:
#     - IQ-Tree (http://www.iqtree.org/)
#     - Julia (https://julialang.org/)
#     - PhyloNetworks.jl (https://github.com/crsl4/PhyloNetworks.jl)
#     - QuartetNetworkGoodnessFit.jl (https://github.com/cecileane/QuartetNetworkGoodnessFit.jl)
# Caitlin Cherryh 2021


##### Step 1: Set file paths and run variables #####
# input_names               <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir and dataset_tree_roots
#                              (e.g. for 2 datasets, put same dataset first and same dataset last for each variable)
# dataset_tree_roots        <- set which taxa is outgroup for each dataset
# alignment_dir             <- the folder(s) containing the alignments for each loci

# compare_ASTRAL_trees      <- set which datasets you want to apply the QuartetNetworkGoF test to
# compare_IQTREE_trees      <- set which datasets you want to apply the AU test to
# tests_to_run              <- a list, with a vector for each dataset specifying which of the recombination detection methods should be tested 
#                              Options: "allTests", "PHI", "maxchi" and "geneconv"
# new.ASTRAL.terminal.branch.length <- ASTRAL does not estimate terminal branch lengths. Terminal branch lengths will be initially set to this value.
# n_julia_reps              <- Number of simulated data sets to generate for the QuartetNetworkGoF test.

# csv_data_dir              <- directory containing the .csv file results from script 1_RecombinationDetection_empiricalTreelikeness.R
# tree_dir                  <- directory containing species trees output from script 2_Species_Tree_Estimation.R
# output_dir                <- where the coalescent/concatenated trees and tree comparisons will be stored 
# main_dir                  <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# iqtree_path               <- location of IQ-Tree executable 

run = "local"

if (run == "local"){
  # Set datasets, which taxa to root the tree at for each dataset, and the location of alignments for each dataset
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  dataset_tree_roots <- c("BAJW", "Apusozoa_Apusozoa_N_A_N_A_N_A_Nutomonas_longa_SRR1617398", "Mus_musculus", "LA4116")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  
  # Set which datasets and which tests to run
  # Set which datasets and which tests to run
  compare_ASTRAL_trees <- c("Vanderpool2020")
  compare_IQTREE_trees <- c("Vanderpool2020")
  tests_to_run <- list("Vanderpool2020" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Pease2016" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Strassert2021" = c(),
                       "1KP" = c())
  new.ASTRAL.terminal.branch.length <- 0.1
  n_julia_reps <- 100
  
  # File and directory locations
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  tree_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
  output_dir <-"/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Software locations
  iqtree_path <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree"
} else if (run == "server"){
  # Set datasets, which taxa to root the tree at for each dataset, and the location of alignments for each dataset
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  dataset_tree_roots <- c("BAJW", "Apusozoa_Apusozoa_N_A_N_A_N_A_Nutomonas_longa_SRR1617398", "Mus_musculus", "LA4116")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  
  # Set which datasets and which tests to run
  compare_ASTRAL_trees <- c("Vanderpool2020")
  compare_IQTREE_trees <- c("Vanderpool2020")
  tests_to_run <- list("Vanderpool2020" = c("allTests", "PHI", "maxchi", "geneconv"),
                    "Pease2016" = c("allTests", "PHI", "maxchi", "geneconv"),
                    "Strassert2021" = c(),
                    "1KP" = c())
  new.ASTRAL.terminal.branch.length <- 0.1
  n_julia_reps <- 100
  
  # File and directory locations
  csv_data_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  tree_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_dataAnalysis/"
  maindir <- "/data/caitlin/empirical_treelikeness/"
  
  # Software locations
  iqtree_path <- "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree"
}



##### Step 2: Open packages and functions #####
# Open packages
library(ape)
# Source the functions using the filepaths
source(paste0(maindir,"code/func_comparison.R"))
source(paste0(maindir,"code/func_analysis.R"))



##### Step 3: Prepare for analysis #####
names(dataset_tree_roots) <- input_names
names(alignment_dir) <- input_names



##### Step 4: Compare ASTRAL  trees #####
# For each test, I have a fail tree, a pass tree. I also have an all tree - one that includes every loci from the dataset
# I want to take my three trees and modify the format to work for QGOF: add arbitrary terminal branches (1) and remove posterior probability values
# Then I need to convert the list of gene trees to a list of gene trees that work for QGOF
# Working in Julia:
#   - Convert the list of gene trees into a set of quartet CFs
#   - Apply the QGOF test and save the output values
#   - Repeat this for the other two trees
# Iterate through each dataset
for (dataset in compare_ASTRAL_trees){
  # Set a folder for the analyses for this dataset
  dataset_folder <- paste0(output_dir, dataset, "/")
  if (dir.exists(dataset_folder) == FALSE){dir.create(dataset_folder)}
  # Extract the names of the ASTRAL species trees and the .txt files containing the list of gene trees
  species_tree_folder <- paste0(tree_dir, dataset, "/", "species_trees", "/")
  all_species_trees_files <- list.files(species_tree_folder, recursive = TRUE)
  
  # Identify which tests to run for this dataset
  dataset_tests <- tests_to_run[[dataset]]
  # Iterate through each test
  for (test in dataset_tests){
    # Name the new folder for running this test
    new_folder <- paste0(dataset_folder, "quarnetGoFtest_", test, "/")
    if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
    
    # Check whether the results file for this test exists already
    quarnet_results_file <- paste0(new_folder, dataset, "_", test, "_QuarNetGoF_test_results.csv")
    # If the file does not exist, run the tests
    if (file.exists(quarnet_results_file) == FALSE){
      # Find all ASTRAL files from this test/dataset combination
      all_astral_files <- grep("ASTRAL", all_species_trees_files, value = TRUE)
      all_astral_trees <- grep("\\.tre", all_astral_files, value = TRUE)
      all_astral_gene_trees <- gsub("_species.tre", ".txt", all_astral_trees)
      
      # Collect files for this test
      files <- c(grep(test, all_astral_trees, value = TRUE), grep("NoTest", all_astral_trees, value = TRUE), 
                 grep("pass", grep(test, all_astral_gene_trees, value = TRUE), value = TRUE))
      
      # Print the dataset and test details
      print(paste0(dataset, " - ", test))
      
      # If all four files exist, continue the analysis
      if (length(files) == 4){
        # Move files into the new folder
        for (i in files){file.copy(from = paste0(species_tree_folder, i), to = paste0(new_folder, i))}
        # Give files their new full file paths
        files <- paste0(new_folder, files)
        
        # Provide extra parameters for the function, depending on dataset
        tree_root = dataset_tree_roots[[dataset]]
        
        # Rewrite ASTRAL species trees to match format for quarnetGoFtest (will all have .tre extension)
        lapply(grep(".tre", files, value = TRUE), reformat.ASTRAL.tree.for.Julia, add.arbitrary.terminal.branches = TRUE, 
               terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Extend the ASTRAL species trees to be ultrametric
        lapply(grep(".tre", files, value = TRUE), make.tree.ultrametric, root.tree = TRUE, outgroup = tree_root)
        # Rewrite IQ-Tree gene trees to match format for quarnetGoFTest (will have .txt extension)
        lapply(grep(".txt", files, value = TRUE), reformat.gene.tree.list.for.Julia, add.arbitrary.terminal.branches = FALSE)
        # Name the files
        noTest_tree_file <- grep("NoTest", grep("tre", files, value = TRUE), value = TRUE)
        pass_tree_file <- grep("pass", grep("tre", files, value = TRUE), value = TRUE)
        fail_tree_file <- grep("fail", grep("tre", files, value = TRUE), value = TRUE)
        gene_trees_file <- grep(".txt", files, value = TRUE)
        
        ### Apply the Quartet Network Goodness of Fit test in Julia ###
        # Write the Julia code into a file
        write.Julia.GoF.script(test_name = test, dataset = dataset, directory = new_folder, pass_tree = pass_tree_file, 
                               fail_tree = fail_tree_file, all_tree = noTest_tree_file, gene_trees = gene_trees_file, 
                               root.species.trees = FALSE, tree_root = tree_root, output_csv_file_path = quarnet_results_file,
                               number_of_simulated_replicates = n_julia_reps)
        # Run the script in Julia to calculate the adequacy of each tree for the quartet concordance factors calculated from the gene trees
        julia_command <- paste0("Julia ",new_folder, "apply_GoF_test.jl")
        system(julia_command)
      } else {
        # For some tests for some datasets, there were no loci that passed/failed that test and so there are not three trees, meaning this analysis
        # cannot be carried out as designed
        # Output a dataframe for this test with NA results
        df <- data.frame(dataset = rep(dataset, 3), concordance_factors = rep(NA,3), test = rep(test, 3), 
                         tree = c("test_pass", "test_fail", "no_test"), tree_root = rep(dataset_tree_roots[[dataset]], 3), 
                         p_value_overall_GoF_test = rep(NA, 3), uncorrected_z_value_test_statistic = rep(NA, 3), 
                         estimated_sigma_for_test_statistic_correction = rep(NA, 3))
        write.csv(df, file = quarnet_results_file)
      }
    }
    
  }
}

## Collate and output ASTRAL results files ##
# Find all QuarNetGoF_test_results.csv files
all_output_files <- list.files(output_dir, recursive = TRUE)
all_gof_results <- paste0(output_dir, grep("QuarNetGoF_test_results.csv", all_output_files, value = TRUE))
# Open and collate the csv files
gof_results_list <- lapply(all_gof_results, read.csv)
gof_results_df <- do.call(rbind, gof_results_list)
# Output compiled csv
gof_results_df_name <- paste0(output_dir, "03_collated_ComparisonTrees_QuarNetGoF_test_results.csv")
write.csv(gof_results_df, file = gof_results_df_name, row.names = FALSE)



##### Step 5: Compare IQ-Tree trees #####
# Iterate through each dataset
for (dataset in compare_IQTREE_trees){
  # Set a folder for the analyses for this dataset
  dataset_folder <- paste0(output_dir, dataset, "/")
  if (dir.exists(dataset_folder) == FALSE){dir.create(dataset_folder)}
  # Extract the names of the ASTRAL species trees and the .txt files containing the list of gene trees
  species_tree_folder <- paste0(tree_dir, dataset, "/", "species_trees", "/")
  all_species_trees_files <- list.files(species_tree_folder, recursive = TRUE)
  
  # Find all IQ-Tree trees and partition files for this dataset
  all_IQTree_files <- grep("IQTREE", all_species_trees_files, value = TRUE)
  all_IQTree_partitions <- grep(".nex.", grep("partitions.nex", all_IQTree_files, value = TRUE), value = TRUE, invert = TRUE)
  all_IQTree_trees <- grep("contree", all_IQTree_files, value = TRUE)
  
  # Identify which tests to run for this dataset
  dataset_tests <- tests_to_run[[dataset]]
  # Iterate through each test
  for (test in dataset_tests){
    # Name the new folder for running this test
    new_folder <- paste0(dataset_folder, "AUtest_", test, "/")
    if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
    
    # Check whether the results file for this test exists already
    au_results_file <- paste0(new_folder, dataset, "_", test, "_AU_test_results.csv")
    # If the file does not exist, run the tests
    if (file.exists(au_results_file) == FALSE){
      # Collect files for this test: three trees and the partition file (containing the loci that pass the test)
      test_IQTREE_trees <- grep(test, all_IQTree_trees, value = TRUE)
      # Sort the files into the following order: test pass, test fail, no test
      three_trees_location <- c(grep("pass", test_IQTREE_trees, value = TRUE), grep("fail", test_IQTREE_trees, value = TRUE), 
                                grep("NoTest", all_IQTree_trees, value = TRUE))
      three_trees_location <- paste0(species_tree_folder, three_trees_location)
      # Read in the three trees
      three_trees_text <- unlist(lapply(three_trees_location, readLines))
      # Write the three trees into one file, inside the new folder for the AU test
      three_trees_path <- paste0(new_folder, dataset, "_", test, "three_trees_Pass-Fail-NoTest.tree")
      write(three_trees_text, three_trees_path)
      
      # Find the partitions file containing the location of the loci that pass the test
      test_pass_partition_file <- grep("pass", grep(test, all_IQTree_partitions, value = TRUE), value = TRUE)
      # Write the new name for the partition file
      partition_path <- paste0(new_folder, "partitions.nex")
      # Copy the partition file
      file.copy(from = paste0(species_tree_folder, test_pass_partition_file), to = partition_path, overwrite = TRUE)
      
      ### Apply the AU tests in IQ-Tree ###
      # Apply the AU test
      au_test_df <- perform.partition.AU.test(partition_path, three_trees_path, iqtree_path)
      # Assemble the output dataframe
      au_results_df <- data.frame(dataset = rep(dataset, 3), test = rep(test, 3), tree = c("test_pass", "test_fail", "no_test"))
      au_results_df <- cbind(au_results_df, au_test_df)
      # Save the output dataframe
      write.csv(au_results_df, file = au_results_file, row.names = FALSE)
    }
  }
}

## Collate and output IQ-Tree results files ##
# Find all AU_test_results.csv files
all_output_files <- list.files(output_dir, recursive = TRUE)
all_au_results <- paste0(output_dir, grep("AU_test_results.csv", all_output_files, value = TRUE))
# Open and collate the csv files
au_results_list <- lapply(all_au_results, read.csv)
au_results_df <- do.call(rbind, au_results_list)
# Output compiled csv
au_results_df_name <- paste0(output_dir, "03_collated_ComparisonTrees_AU_test_results.csv")
write.csv(au_results_df, file = au_results_df_name, row.names = FALSE)



