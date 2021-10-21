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
  input_names <- c("1KP", "Whelan2017","Vanderpool2020", "Pease2016")
  dataset_tree_roots <- c("BAJW", "Salpingoeca_pyxidium", "Mus_musculus", "LA4116")
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Whelan2017/genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  
  # Set which datasets and which tests to run
  compare_ASTRAL_trees <- c("Whelan2017")
  compare_IQTREE_trees <- c()
  tests_to_run <- list("Vanderpool2020" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Pease2016" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Whelan2017" = c("PHI", "maxchi", "geneconv"),
                       "1KP" = c("PHI", "maxchi"))
  new.ASTRAL.terminal.branch.length <- 0.1
  n_julia_reps <- 100
  
  # File and directory locations
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  tree_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
  output_dir <-"/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/" # where the empirical treelikeness code is
  
  # Software locations
  iqtree_path <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree"
} else if (run == "server"){
  # Set datasets, which taxa to root the tree at for each dataset, and the location of alignments for each dataset
  input_names <- c("1KP", "Whelan2017","Vanderpool2020", "Pease2016")
  dataset_tree_roots <- c("BAJW", "Salpingoeca_pyxidium", "Mus_musculus", "LA4116")
  alignment_dir <- c("/data/caitlin/empirical_treelikeness/Data_1KP/",
                     "/data/caitlin/empirical_treelikeness/Data_Whelan2017/",
                     "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/",
                     "/data/caitlin/empirical_treelikeness/Data_Pease2016/")
  
  # Set which datasets and which tests to run
  compare_ASTRAL_trees <- c()
  compare_IQTREE_trees <- c("1KP")
  tests_to_run <- list("Vanderpool2020" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Pease2016" = c("allTests", "PHI", "maxchi", "geneconv"),
                       "Whelan2017" = c("PHI", "maxchi", "geneconv"),
                       "1KP" = c("PHI", "maxchi"))
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
    print(paste0("Applying QuarNet GoF test to ", dataset, ": ", test))
    # Name the new folder for running this test
    new_folder <- paste0(dataset_folder, "quarnetGoFtest_", test, "/")
    if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
    
    # Find all ASTRAL files from this test/dataset combination
    all_astral_files <- grep("ASTRAL", all_species_trees_files, value = TRUE)
    all_astral_trees <- grep("\\.tre", all_astral_files, value = TRUE)
    all_astral_gene_trees <- gsub("_species.tre", ".txt", all_astral_trees)
    
    # Collect files for this test
    files <- c(grep(test, all_astral_trees, value = TRUE), grep("NoTest", all_astral_trees, value = TRUE), 
               grep("pass", grep(test, all_astral_gene_trees, value = TRUE), value = TRUE))
    old_files <- paste0(species_tree_folder, files)
    # Move files into the new folder
    for (i in files){file.copy(from = paste0(species_tree_folder, i), to = paste0(new_folder, i))}
    # Give files their new full file paths
    files <- paste0(new_folder, files)

    # Name the files
    noTest_tree_file <- grep("NoTest", grep("tre", files, value = TRUE), value = TRUE)
    pass_tree_file <- grep("pass", grep("tre", files, value = TRUE), value = TRUE)
    if (dataset == "Vanderpool2020" | dataset == "Pease2016" |(dataset == "Whelan2017" & test == "geneconv")){
      # Only collect fail tree for Vanderpool2020 and Pease2016 datasets
      # The fail trees for the other two datasets contained too few loci (due to tests not working for loci with few variable sites)
      fail_tree_file <- grep("fail", grep("tre", files, value = TRUE), value = TRUE)
    }
    gene_trees_file <- grep(".txt", files, value = TRUE)
    
    # Create a new vector of the trees used for this dataset
    if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
      dataset_tree_files <- c(pass_tree_file, fail_tree_file, noTest_tree_file)
    } else if (dataset == "Whelan2017" & test == "geneconv"){
      dataset_tree_files <- c(pass_tree_file, fail_tree_file, noTest_tree_file)
    } else if (dataset == "1KP" | (dataset == "Whelan2017" & test != "geneconv")){
      dataset_tree_files <- c(pass_tree_file, noTest_tree_file)
    }
    
    # Provide extra parameters for the function, depending on dataset
    tree_root = dataset_tree_roots[[dataset]]
    
    # Rewrite ASTRAL species trees to match format for quarnetGoFtest (will all have .tre extension)
    lapply(dataset_tree_files, reformat.ASTRAL.tree.for.Julia, add.arbitrary.terminal.branches = TRUE, 
           terminal.branch.length = new.ASTRAL.terminal.branch.length)
    # Extend the ASTRAL species trees to be ultrametric
    lapply(dataset_tree_files, make.tree.ultrametric, root.tree = TRUE, outgroup = tree_root)
    # Rewrite IQ-Tree gene trees to match format for quarnetGoFTest (will have .txt extension)
    lapply(gene_trees_file, reformat.gene.tree.list.for.Julia, add.arbitrary.terminal.branches = FALSE)
    
    # Assemble the output csv name for the Quartet Network Goodness of Fit results
    quarnet_results_file <- paste0(new_folder, dataset, "_", test, "_QuarNetGoF_test_results.csv")
    
    # Run slightly different versions of the code based on the dataset: shallow datasets compare 3 trees, and deeper datasets compare two trees
    if (dataset == "Vanderpool2020" | dataset == "Pease2016" | (dataset == "Whelan2017" & test == "geneconv")){
      # Assemble the output csv name for the Quartet Network Goodness of Fit results
      quarnet_results_file <- paste0(new_folder, dataset, "_", test, "_QuarNetGoF_test_results.csv")
      
      if (file.exists(quarnet_results_file) == FALSE){
        ### Run GoodnessOfFit test
        # Run Julia script with T_test,pass and T_test,fail and T_none
        write.Julia.GoF.script(test_name = test, dataset = dataset, directory = new_folder, pass_tree = pass_tree_file, 
                               fail_tree = fail_tree_file, all_tree = noTest_tree_file, gene_trees = gene_trees_file, 
                               root.species.trees = FALSE, tree_root = tree_root, output_csv_file_path = quarnet_results_file,
                               number_of_simulated_replicates = n_julia_reps)
        # Run the script in Julia to calculate the adequacy of each tree for the quartet concordance factors calculated from the gene trees
        julia_command <- paste0("Julia ",new_folder, "apply_GoF_test.jl")
        system(julia_command)
      }
      
      # Assemble the file name for the RF distances csv file
      rf_csv <- paste0(new_folder, dataset, "_", test, "_ASTRAL_tree_RF_distances.csv")
      if (file.exists(rf_csv) == FALSE){
        ### Calculate the RF and wRF distances between trees
        # Copy the test,pass tree across again and add terminal branches 
        t_test_pass_file <- paste0(new_folder, dataset, "_", test, "_pass_ASTRAL_terminalBranches0.1.tre")
        file.copy(from = grep("pass", grep("\\.tre", old_files, value = TRUE), value = TRUE), to = t_test_pass_file)
        reformat.ASTRAL.tree.for.Julia(t_test_pass_file, add.arbitrary.terminal.branches = TRUE, terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Copy the test,fail tree across again and add terminal branches 
        t_test_fail_file <- paste0(new_folder, dataset, "_", test, "_fail_ASTRAL_terminalBranches0.1.tre")
        file.copy(from = grep("fail", grep("\\.tre", old_files, value = TRUE), value = TRUE), to = t_test_fail_file)
        reformat.ASTRAL.tree.for.Julia(t_test_fail_file, add.arbitrary.terminal.branches = TRUE, terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Copy the NoTest tree across again and add terminal branches 
        t_none_file <- paste0(new_folder, dataset, "_noTest_ASTRAL_terminalBranches0.1.tre")
        file.copy(from = grep("NoTest", grep("\\.tre", old_files, value = TRUE), value = TRUE), to = t_none_file)
        reformat.ASTRAL.tree.for.Julia(t_none_file, add.arbitrary.terminal.branches = TRUE, terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Read in both trees
        t_test_pass <- read.tree(t_test_pass_file)
        t_test_fail <- read.tree(t_test_fail_file)
        t_none <- read.tree(t_none_file)
        # Calculate distances between the trees
        dist_df <- data.frame(dataset = rep(dataset, 3),
                              test = rep(test, 3),
                              tree = c("test_pass", "test_fail", "none"),
                              analysis = rep("ASTRAL", 3),
                              RF_dist_to_test_pass = c(RF.dist(t_test_pass, t_test_pass, check.labels = TRUE), RF.dist(t_test_pass, t_test_fail, check.labels = TRUE), RF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              RF_dist_to_test_fail = c(RF.dist(t_test_fail, t_test_pass, check.labels = TRUE), RF.dist(t_test_fail, t_test_fail, check.labels = TRUE), RF.dist(t_test_fail, t_none, check.labels = TRUE)),
                              RF_dist_to_test_none = c(RF.dist(t_none, t_test_pass, check.labels = TRUE),      RF.dist(t_none, t_test_fail, check.labels = TRUE),      RF.dist(t_none, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_pass = c(wRF.dist(t_test_pass, t_test_pass, check.labels = TRUE), wRF.dist(t_test_pass, t_test_fail, check.labels = TRUE), wRF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_fail = c(wRF.dist(t_test_fail, t_test_pass, check.labels = TRUE), wRF.dist(t_test_fail, t_test_fail, check.labels = TRUE), wRF.dist(t_test_fail, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_none = c(wRF.dist(t_none, t_test_pass, check.labels = TRUE),      wRF.dist(t_none, t_test_fail, check.labels = TRUE),      wRF.dist(t_none, t_none, check.labels = TRUE)))
        # Save RF distance as a dataframe
        write.csv(dist_df, file = rf_csv, row.names = FALSE)
      }
      
    } else if (dataset == "1KP" | (dataset == "Whelan2017" & test != "geneconv")){
      # Assemble the output csv name for the Quartet Network Goodness of Fit results
      quarnet_results_file <- paste0(new_folder, dataset, "_", test, "_QuarNetGoF_test_results.csv")
      
      if (file.exists(quarnet_results_file) == FALSE){
        ### Run GoodnessOfFit test
        # Run Julia script with T_test,pass and T_none
        write.Julia.GoF.script.two.trees(test_name = test, dataset = dataset, directory = new_folder, pass_tree = pass_tree_file, 
                                         all_tree = noTest_tree_file, gene_trees = gene_trees_file,root.species.trees = FALSE, tree_root = tree_root,
                                         output_csv_file_path = quarnet_results_file, number_of_simulated_replicates = n_julia_reps)
        # Run the script in Julia to calculate the adequacy of each tree for the quartet concordance factors calculated from the gene trees
        julia_command <- paste0("Julia ",new_folder, "apply_GoF_test.jl")
        system(julia_command)
      }
      
      # Assemble the file name for the RF distances csv file
      rf_csv <- paste0(new_folder, dataset, "_", test, "_ASTRAL_tree_RF_distances.csv")
      if (file.exists(rf_csv) == FALSE){
        ### Calculate the RF and wRF distances between trees
        # Copy the test,pass tree across again and add terminal branches 
        t_test_pass_file <- paste0(new_folder, dataset, "_", test, "_pass_ASTRAL_terminalBranches0.1.tre")
        file.copy(from = grep("pass", grep("\\.tre", old_files, value = TRUE), value = TRUE), to = t_test_pass_file)
        reformat.ASTRAL.tree.for.Julia(t_test_pass_file, add.arbitrary.terminal.branches = TRUE, terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Copy the NoTest tree across again and add terminal branches 
        t_none_file <- paste0(new_folder, dataset, "_noTest_ASTRAL_terminalBranches0.1.tre")
        file.copy(from = grep("NoTest", grep("\\.tre", old_files, value = TRUE), value = TRUE), to = t_none_file)
        reformat.ASTRAL.tree.for.Julia(t_none_file, add.arbitrary.terminal.branches = TRUE, terminal.branch.length = new.ASTRAL.terminal.branch.length)
        # Read in both trees
        t_test_pass <- read.tree(t_test_pass_file)
        t_none <- read.tree(t_none_file)
        # If the number of tips are different, drop the tips that aren't included in the test_pass tree
        if ( length(t_none$tip.label) != length(t_test_pass$tip.label)){
          keep_tips <- t_test_pass$tip.label
          t_none <- keep.tip(t_none, keep_tips)
        }
        # Calculate distances between the trees
        dist_df <- data.frame(dataset = rep(dataset, 2),
                              test = rep(test, 2),
                              analysis = rep("ASTRAL", 2),
                              tree = c("test_pass", "none"),
                              RF_dist_to_test_pass = c(RF.dist(t_test_pass, t_test_pass, check.labels = TRUE), RF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              RF_dist_to_test_fail = c(NA,NA),
                              RF_dist_to_test_none = c(RF.dist(t_none, t_test_pass, check.labels = TRUE), RF.dist(t_none, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_pass = c(wRF.dist(t_test_pass, t_test_pass, check.labels = TRUE), wRF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_fail = c(NA,NA),
                              wRF_dist_to_test_none = c(wRF.dist(t_none, t_test_pass, check.labels = TRUE),  wRF.dist(t_none, t_none, check.labels = TRUE)))
        # Save RF distance as a dataframe
        write.csv(dist_df, file = rf_csv, row.names = FALSE)
      } # end check for rf_csv file
    } # end else if dataset == deep
    
  } # end iterating through tests
  
} # end iterating through datasets




##### Collate and output ASTRAL results files #####
# Find all QuarNetGoF_test_results.csv files (and remove the collated results file from the list of files to combine)
all_output_files <- list.files(output_dir, recursive = TRUE)
all_gof_results <- grep("QuarNetGoF_test_results.csv", all_output_files, value = TRUE)
all_gof_results <- grep("collated", all_gof_results, value = TRUE, invert = TRUE)
if (length(all_gof_results) > 0){
  print("Collating QuarNet GoF test results")
  # Attach directory name to file names
  all_gof_results <- paste0(output_dir, all_gof_results)
  # Open and collate the csv files
  gof_results_list <- lapply(all_gof_results, read.csv)
  gof_results_df <- do.call(rbind, gof_results_list)
  # Output compiled csv
  gof_results_df_name <- paste0(output_dir, "03_AllDatasets_collated_ComparisonTrees_QuarNetGoF_test_results.csv")
  write.csv(gof_results_df, file = gof_results_df_name, row.names = FALSE)
}




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
  
  # If dataset was run using RAxML-NG, find the bestTree files (the trees) and the IQ-Tree partitions for the no free rate model runs
  if (dataset == "1KP"){
    # Find RAxML-NG trees
    all_raxml_trees <- grep("bestTree", all_species_trees_files, value = TRUE)
    # Replace IQ-Tree trees with RAxML-NG trees
    all_IQTree_trees <- all_raxml_trees
    # Find the partition files for the no free rates runs
    all_IQTree_partitions <- grep("noFreeRates", all_IQTree_partitions, value = TRUE)
  }
  
  # Identify which tests to run for this dataset
  dataset_tests <- tests_to_run[[dataset]]
  # Iterate through each test
  for (test in dataset_tests){
    print(paste0("Applying AU test to ", dataset, ": ", test))
    # Name the new folder for running this test
    new_folder <- paste0(dataset_folder, "AUtest_", test, "/")
    if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
    
    # Assemble the filename for the results file 
    au_results_file <- paste0(new_folder, dataset, "_", test, "_AU_test_results.csv")
    # If the file does not exist, run the tests
    if (dataset == "Vanderpool2020" | dataset == "Pease2016" | ((dataset == "Whelan2017") & (test == "geneconv"))){
      ### Apply the AU test
      # Check whether the results file for this test exists already
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
        
        # Apply the AU test
        au_test_df <- perform.partition.AU.test(partition_path, three_trees_path, iqtree_path)
        # Assemble the output dataframe
        au_results_df <- data.frame(dataset = rep(dataset, 3), test = rep(test, 3), tree = c("test_pass", "test_fail", "no_test"))
        au_results_df <- cbind(au_results_df, au_test_df)
        # Save the output dataframe
        write.csv(au_results_df, file = au_results_file, row.names = FALSE)
      }
      
      ### Calculate the RF and wRF distances between trees
      # Read in trees
      rf_csv <- paste0(new_folder, dataset, "_", test, "_IQTREE_tree_RF_distances.csv")
      if (file.exists(rf_csv) == FALSE){
        three_trees <- read.tree(file = three_trees_path)
        # Calculate distances between the trees
        t_test_pass <- three_trees[[1]]
        t_test_fail <- three_trees[[2]]
        t_none <- three_trees[[3]]
        # Calculate RF/wRF distances
        dist_df <- data.frame(dataset = rep(dataset, 3), 
                              test = rep(test, 3), 
                              tree = c("test_pass", "test_fail", "no_test"), 
                              analysis = rep("IQ-Tree", 3),
                              RF_dist_to_test_pass = c(RF.dist(t_test_pass, t_test_pass, check.labels = TRUE), RF.dist(t_test_pass, t_test_fail, check.labels = TRUE), RF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              RF_dist_to_test_fail = c(RF.dist(t_test_fail, t_test_pass, check.labels = TRUE), RF.dist(t_test_fail, t_test_fail, check.labels = TRUE), RF.dist(t_test_fail, t_none, check.labels = TRUE)),
                              RF_dist_to_test_none = c(RF.dist(t_none, t_test_pass, check.labels = TRUE),      RF.dist(t_none, t_test_fail, check.labels = TRUE),      RF.dist(t_none, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_pass = c(wRF.dist(t_test_pass, t_test_pass, check.labels = TRUE), wRF.dist(t_test_pass, t_test_fail, check.labels = TRUE), wRF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_fail = c(wRF.dist(t_test_fail, t_test_pass, check.labels = TRUE), wRF.dist(t_test_fail, t_test_fail, check.labels = TRUE), wRF.dist(t_test_fail, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_none = c(wRF.dist(t_none, t_test_pass, check.labels = TRUE),      wRF.dist(t_none, t_test_fail, check.labels = TRUE),      wRF.dist(t_none, t_none, check.labels = TRUE)))
        # Save the output dataframe
        write.csv(dist_df, file = rf_csv, row.names = FALSE)
      }
      
    } else if (dataset == "1KP" | (dataset == "Whelan2017" & test != "geneconv")){
      ### Apply the AU test
      # Check whether the results file for this test exists already
      if (file.exists(au_results_file) == FALSE){
        # Collect files for this test: two trees and the partition file (containing the loci that pass the test)
        test_IQTREE_trees <- grep(test, all_IQTree_trees, value = TRUE)
        # Sort the files into the following order: test pass, no test
        two_trees_location <- c(grep("pass", test_IQTREE_trees, value = TRUE), grep("NoTest", all_IQTree_trees, value = TRUE))
        two_trees_location <- paste0(species_tree_folder, two_trees_location)
        # Read in the two trees
        two_trees_text <- unlist(lapply(two_trees_location, readLines))
        # Write the two trees into one file, inside the new folder for the AU test
        two_trees_path <- paste0(new_folder, dataset, "_", test, "_two_trees_Pass-NoTest.tree")
        write(two_trees_text, two_trees_path)
        
        # Find the partitions file containing the location of the loci that pass the test
        test_pass_partition_file <- grep("pass", grep(test, all_IQTree_partitions, value = TRUE), value = TRUE)
        # Write the new name for the partition file
        partition_path <- paste0(new_folder, "partitions.nex")
        # Copy the partition file
        file.copy(from = paste0(species_tree_folder, test_pass_partition_file), to = partition_path, overwrite = TRUE)
        
        # Apply the AU test
        au_test_df <- perform.partition.AU.test.two.trees(partition_path, two_trees_path, iqtree_path)
        # Assemble the output dataframe
        au_results_df <- data.frame(dataset = rep(dataset, 2), test = rep(test, 2), tree = c("test_pass", "no_test"))
        au_results_df <- cbind(au_results_df, au_test_df)
        # Save the output dataframe
        write.csv(au_results_df, file = au_results_file, row.names = FALSE)
      }
      
      ### Calculate the RF and wRF distances between trees
      rf_csv <- paste0(new_folder, dataset, "_", test, "_IQTREE_tree_RF_distances.csv")
      if (file.exists(rf_csv) == FALSE){
        # Read in trees
        two_trees <- read.tree(file = two_trees_path)
        # Calculate distances between the trees
        t_test_pass <- two_trees[[1]]
        t_none <- two_trees[[2]]
        # Calculate RF/wRF distances
        dist_df <- data.frame(dataset = rep(dataset, 2), 
                              test = rep(test, 2), 
                              tree = c("test_pass", "no_test"), 
                              analysis = rep("IQ-Tree", 2),
                              RF_dist_to_test_pass = c(RF.dist(t_test_pass, t_test_pass, check.labels = TRUE), RF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              RF_dist_to_test_fail = c(NA,NA),
                              RF_dist_to_test_none = c(RF.dist(t_none, t_test_pass, check.labels = TRUE), RF.dist(t_none, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_pass = c(wRF.dist(t_test_pass, t_test_pass, check.labels = TRUE), wRF.dist(t_test_pass, t_none, check.labels = TRUE)),
                              wRF_dist_to_test_fail = c(NA,NA),
                              wRF_dist_to_test_none = c(wRF.dist(t_none, t_test_pass, check.labels = TRUE), wRF.dist(t_none, t_none, check.labels = TRUE)))
        # Save the output dataframe
        write.csv(dist_df, file = rf_csv, row.names = FALSE)
      } # end check for RF file
    } #end deep dataset
    
  } # end iterating through tests
  
} # end iterating through datasets




##### Collate and output IQ-Tree results files #####
# Find all AU_test_results.csv files (and remove the collated results file from the list of files to combine)
all_output_files <- list.files(output_dir, recursive = TRUE)
all_au_results <- grep("AU_test_results.csv", all_output_files, value = TRUE)
all_au_results <- grep("collated", all_au_results, value = TRUE, invert = TRUE)
if (length(all_au_results) > 0){
  print("Collating AU test results")
  # Attach directory name to file paths
  all_au_results <- paste0(output_dir, all_au_results)
  # Open and collate the csv files
  au_results_list <- lapply(all_au_results, read.csv)
  au_results_df <- do.call(rbind, au_results_list)
  # Output compiled csv
  au_results_df_name <- paste0(output_dir, "03_AllDatasets_collated_ComparisonTrees_AU_test_results.csv")
  write.csv(au_results_df, file = au_results_df_name, row.names = FALSE)
}




##### Collate and output all RF/wRF distance results files #####
# Find all tree_RF_distances.csv files (and remove the collated results file from the list of files to combine)
all_output_files <- list.files(output_dir, recursive = TRUE)
all_rf_results <- grep("tree_RF_distances.csv", all_output_files, value = TRUE)
all_rf_results <- grep("collated", all_rf_results, value = TRUE, invert = TRUE)
all_rf_results <- all_rf_results[grepl(paste(input_names, collapse = "|"), all_rf_results)]
if (length(all_rf_results) > 0){
  print("Collating RF/wRF distance results")
  # Attach directory name to file names
  all_rf_results <- paste0(output_dir, all_rf_results)
  # Open and collate the csv files
  rf_results_list <- lapply(all_rf_results, read.csv)
  rf_results_df <- do.call(rbind, rf_results_list)
  # Output compiled csv
  rf_results_df_name <- paste0(output_dir, "03_AllDatasets_collated_RF_wRF_distances_results.csv")
  write.csv(rf_results_df, file = rf_results_df_name, row.names = FALSE)
}

