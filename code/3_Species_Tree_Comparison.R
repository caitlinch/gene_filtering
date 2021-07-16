### empirical_treelikeness/code/3_Species_Tree_Comparison.R
## R program to identify the best fitting tree for a candidate dataset
## Additional software/packages required:
#     - IQ-Tree (http://www.iqtree.org/)
#     - Julia (https://julialang.org/)
#     - PhyloNetworks.jl (https://github.com/crsl4/PhyloNetworks.jl)
#     - QuartetNetworkGoodnessFit.jl (https://github.com/cecileane/QuartetNetworkGoodnessFit.jl)
# Caitlin Cherryh 2021


##### Step 1: Set file paths and run variables #####
run = "local"

if (run == "local"){
  # Set datasets, which taxa to root the tree at for each dataset, and the location of alignments for each dataset
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  dataset_tree_roots <- c("", "", "Mus_musculus", "LA4116")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  
  # Set which datasets and which tests to run
  datasets_to_run <- c("Vanderpool2020")
  tests_to_run <- c("allTests", "PHI", "maxchi", "geneconv")
  
  # File and directory locations
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  tree_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  output_dir <-c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  
  # Software locations
  iqtree_path <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree"
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
for (dataset in datasets_to_run){
  # Set a folder for the analyses for this dataset
  dataset_folder <- paste0(output_dir, dataset, "/")
  if (dir.exists(dataset_folder) == FALSE){dir.create(dataset_folder)}
  # Extract the names of the ASTRAL species trees and the .txt files containing the list of gene trees
  species_tree_folder <- paste0(tree_dir, dataset, "/", "species_trees", "/")
  all_species_trees_files <- list.files(species_tree_folder, recursive = TRUE)
  
  # Iterate through each test
  for (test in tests_to_run){
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
      # Move files into the new folder
      for (i in files){file.copy(from = paste0(species_tree_folder, i), to = paste0(new_folder, i))}
      # Give files their new full file paths
      files <- paste0(new_folder, files)
      
      # Rewrite ASTRAL trees to match format for quarnetGoFtest (will all have .tre extension)
      lapply(grep(".tre", files, value = TRUE), reformat.ASTRAL.tree.for.Julia)
      # Rewrite ASTRAL gene trees to match format for quarnetGoFTest (will have .txt extension)
      lapply(grep(".txt", files, value = TRUE), reformat.gene.tree.list.for.Julia, gene.tree.source = "IQ-TREE")
      # Name the files
      noTest_tree_file <- grep("NoTest", grep("tre", files, value = TRUE), value = TRUE)
      pass_tree_file <- grep("pass", grep("tre", files, value = TRUE), value = TRUE)
      fail_tree_file <- grep("fail", grep("tre", files, value = TRUE), value = TRUE)
      gene_trees_file <- grep(".txt", files, value = TRUE)
      # Provide extra parameters for the function, depending on dataset
      tree_root = dataset_tree_roots[[dataset]]
      
      ### Apply the Quartet Network Goodness of Fit test in Julia ###
      # Write the Julia code into a file
      write.Julia.GoF.script(test_name = test, dataset = dataset, directory = new_folder, pass_tree = pass_tree_file, fail_tree = fail_tree_file, 
                             all_tree = noTest_tree_file, gene_trees = gene_trees_file, tree_root = tree_root, output_csv_file_path = quarnet_results_file)
      # Run the script in Julia to calculate the adequacy of each tree for the quartet concordance factors calculated from the gene trees
      julia_command <- paste0("Julia ",new_folder, "apply_GoF_test.jl")
      system(julia_command)
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
gof_results_df_name <- paste0(output_dir, "03_",paste(sort(datasets_to_run), collapse = "_"), "_collated_ComparisonTrees_QuarNetGoF_test_results.csv")
write.csv(gof_results_df, file = gof_results_df_name, row.names = FALSE)



##### Step 5: Compare IQ-Tree trees #####
# Iterate through each dataset
for (dataset in datasets_to_run){
  # Set a folder for the analyses for this dataset
  dataset_folder <- paste0(output_dir, dataset, "/")
  if (dir.exists(dataset_folder) == FALSE){dir.create(dataset_folder)}
  # Extract the names of the ASTRAL species trees and the .txt files containing the list of gene trees
  species_tree_folder <- paste0(tree_dir, dataset, "/", "species_trees", "/")
  all_species_trees_files <- list.files(species_tree_folder, recursive = TRUE)
  
  for (test in tests_to_run){
    # Name the new folder for running this test
    new_folder <- paste0(dataset_folder, "AUtest_", test, "/")
    if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
    
    # Check whether the results file for this test exists already
    au_results_file <- paste0(new_folder, dataset, "_", test, "_AU_test_results.csv")
    # If the file does not exist, run the tests
    if (file.exists(au_results_file) == FALSE){
      # Find all IQ-Tree files from this test/dataset combination
      all_IQTree_files <- grep("IQTREE", all_species_trees_files, value = TRUE)
      all_IQTree_partitions <- grep(".nex.", grep("partitions.nex", all_IQTree_files, value = TRUE), value = TRUE, invert = TRUE)
      all_IQTree_trees <- gsub(".nex", ".nex.contree", all_IQTree_partitions)
      
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
      file.copy(from = paste0(species_tree_folder, test_pass_partition_file), to = partition_path)
      
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
au_results_df_name <- paste0(output_dir, "03_",paste(sort(datasets_to_run), collapse = "_"), "_collated_ComparisonTrees_AU_test_results.csv")
write.csv(au_results_df, file = au_results_df_name, row.names = FALSE)



