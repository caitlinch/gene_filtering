### empirical_treelikeness/3_DataAnalysis.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
library(ape)
library(distory)
#library(ggplot2)
#library(adegenet)
#library(treespace)
#library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
# these methods all assume an unrooted tree so trees can be used as is for this analysis



##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
print("set filepaths")
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# tree_data_dir     <- Location of the gene trees
# test_data_dir     <- Location of the results from the AU test and QuartetNetwork Goodness of Fit tests
# output_dir        <- for saving collated output and results from treelikeness analysis.
#
# input_names               <- set name(s) for the dataset(s) - make sure input_names is in same order as alignment_dir and dataset_tree_roots
#                              (e.g. for 2 datasets, put same dataset first and same dataset last for each variable)
# dataset_tree_roots        <- set which taxa is outgroup for each dataset
# alignment_dir             <- the folder(s) containing the alignments for each loci
# tests_to_run              <- a list, with a vector for each dataset specifying which of the recombination detection methods should be tested 
#                              Options: "allTests", "PHI", "maxchi" and "geneconv"

maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
tree_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
test_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"

input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
dataset_tree_roots <- c("BAJW", "Apusozoa_Apusozoa_N_A_N_A_N_A_Nutomonas_longa_SRR1617398", "Mus_musculus", "LA4116")
alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                   "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                   "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                   "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
tests_to_run <- list("Vanderpool2020" = c("allTests", "PHI", "maxchi", "geneconv"),
                     "Pease2016" = c("allTests", "PHI", "maxchi", "geneconv"),
                     "Strassert2021" = c("PHI", "maxchi"),
                     "1KP" = c("PHI", "maxchi"))

datasets_to_identify_distinct_edges <- c("Pease2016", "Vanderpool2020")


##### Step 3: Source function files and prepare variables for analysis #####
print("Source function files")
source(paste0(maindir, "code/func_analysis.R"))

print("Prepare variables for indexing")
# Name vectors for tree roots and alignment locations so they can be accessed via index
names(dataset_tree_roots) <- input_names
names(alignment_dir) <- input_names



##### Step 4: Compare the posterior probabilities/ bootstraps of the trees #####
node_output_dir <- paste0(output_dir, "node_comparisons/")
if (dir.exists(node_output_dir) == FALSE){
  dir.create(node_output_dir)
}

for (dataset in datasets_to_identify_distinct_edges){
  # Identify file containing species trees
  dataset_tree_dir <- paste0(tree_data_dir, dataset, "/species_trees/")
  all_files <- list.files(dataset_tree_dir, recursive = TRUE)
  
  dataset_tests <- tests_to_run[[dataset]]
  # Iterate through each test and identify that information
  for (test in dataset_tests){
    print(paste0("Processing ", dataset, ": ", test))
    ## IQ-Tree trees: Create dataframe detailing differences in posterior probabilities between the two trees
    # Get the list of trees estimated in IQ-Tree for this dataset
    test_trees <- grep(test, all_files, value = TRUE)
    iq_trees <- grep(".contree", test_trees, value = TRUE)
    # Make the full filepaths for each of the three trees (test pass, test fail, and no test)
    none_tree_file <- paste0(dataset_tree_dir, grep(".contree", grep("NoTest", all_files, value = TRUE), value = TRUE))
    pass_tree_file <- paste0(dataset_tree_dir, grep("pass", iq_trees, value = TRUE))
    fail_tree_file <- paste0(dataset_tree_dir, grep("fail", iq_trees, value = TRUE))
    # Create dataframes
    # Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
    test_df_pass_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, tree1_name = "Pass", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "BS")
    test_df_fail_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, tree1_name = "Fail", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "BS")
    
    ## ASTRAL trees: Create dataframe detailing differences in posterior probabilities between the two trees
    # Get the list of trees estimated in ASTRAL for this dataset
    test_trees <- grep(test, all_files, value = TRUE)
    astral_trees <- grep(".tre", grep(".ASTRAL", test_trees, value = TRUE), value = TRUE)
    pass_tree_file <- paste0(dataset_tree_dir, grep("pass", astral_trees, value = TRUE))
    fail_tree_file <- paste0(dataset_tree_dir, grep("fail", astral_trees, value = TRUE))
    none_tree_file <- paste0(dataset_tree_dir, grep("NoTest", grep(".tre", grep(".ASTRAL", all_files, value = TRUE), value = TRUE), value = TRUE))
    # Create dataframes
    # Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
    test_df_pass_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, tree1_name = "Pass", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
    test_df_fail_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, tree1_name = "Fail", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
    
    # Combine all four dataframes into one
    test_df <- rbind(test_df_pass_iq, test_df_fail_iq, test_df_pass_astral, test_df_fail_astral)
    # Save dataset
    test_df_filename <- paste0(node_output_dir, dataset, "_", test, "_ExtractDistinctEdges.csv")
    write.csv(test_df, file = test_df_filename)
  }
}

# Collate all dataframes
all_csvs <- list.files(node_output_dir)
all_csvs <- grep("Collated", all_csvs, value = TRUE, invert = TRUE)
all_csvs <- paste0(node_output_dir, all_csvs)
all_csv_dfs <- lapply(all_csvs, read.csv)
node_df <- as.data.frame(do.call(rbind, all_csv_dfs))
node_df_filename <- paste0(node_output_dir, "Collated_ExtractDistinctEdges.csv")
write.csv(node_df, file = node_df_filename)



