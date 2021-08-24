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
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
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

treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
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



##### Step 3: Prepare variables for analysis #####
# Name vectors for tree roots and alignment locations so they can be accessed via index
names(dataset_tree_roots) <- input_names
names(alignment_dir) <- input_names



##### Step 4: Compare the posterior probabilities of the trees #####
dataset = "Pease2016"
test = "PHI"

dataset_tree_dir <- paste0(tree_data_dir, dataset, "/species_trees/")
all_files <- list.files(dataset_tree_dir, recursive = TRUE)
# Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
# Iterate through each test and identify that information

# Get the list of trees estimated in IQ-Tree for this dataset
test_trees <- grep(test, all_files, value = TRUE)
iq_trees <- grep(".contree", test_trees, value = TRUE)
# Make the full filepaths for each of the three trees (test pass, test fail, and no test)
none_tree_file <- paste0(dataset_tree_dir, grep(".contree", grep("NoTest", all_files, value = TRUE), value = TRUE))
pass_tree_file <- paste0(dataset_tree_dir, grep("pass", iq_trees, value = TRUE))
fail_tree_file <- paste0(dataset_tree_dir, grep("fail", iq_trees, value = TRUE))
# Read each in as a tree
t_none <- read.tree(none_tree_file)
t_pass <- read.tree(pass_tree_file)
t_fail <- read.tree(fail_tree_file)
# Compare splits on test trees with t_none
# Determine which edges are contained in one tree and not the other: returns a numeric vector of edge ids for the first tree
e_pass_none <- distinct.edges(t_pass, t_none)
e_none_pass <- distinct.edges(t_none, t_pass)
e_fail_none <- distinct.edges(t_fail, t_none)
e_none_fail <- distinct.edges(t_none, t_fail)
# Collect information about each of those different edges
if (length(e_pass_none) > 0){
  e_pass_none_list <- lapply(e_pass_none, get.edge.details, t_pass)
  e_pass_none_df <- as.data.frame(do.call(rbind, e_pass_none_list))
  e_pass_none_df$tree1 <- "Pass"
  e_pass_none_df$tree2 <- "None"
  e_pass_none_df$test <- test
  e_pass_none_df$dataset <- dataset
  e_pass_none_df$support_value_type <- "UFB"
} else {
  e_pass_none_df <- data.frame()
}

if (length(e_none_pass) > 0){
  e_none_pass_list <- lapply(e_none_pass, get.edge.details, t_none)
  e_none_pass_df <- as.data.frame(do.call(rbind, e_none_pass_list))
  e_none_pass_df$tree1 <- "None"
  e_none_pass_df$tree2 <- "Pass"
  e_none_pass_df$test <- test
  e_none_pass_df$dataset <- dataset
  e_none_pass_df$support_value_type <- "UFB"
} else {
  e_none_pass_df <- data.frame()
}

if (length(e_fail_none) > 0){
  e_fail_none_list <- lapply(e_fail_none, get.edge.details, t_fail)
  e_fail_none_df <- as.data.frame(do.call(rbind, e_fail_none_list))
  e_fail_none_df$tree1 <- "Fail"
  e_fail_none_df$tree2 <- "None"
  e_fail_none_df$test <- test
  e_fail_none_df$dataset <- dataset
  e_fail_none_df$support_value_type <- "UFB"
}  else {
  e_fail_none_df <- data.frame()
}

if (length(e_none_fail) > 0){
  e_none_fail_list <- lapply(e_none_fail, get.edge.details, t_none)
  e_none_fail_df <- as.data.frame(do.call(rbind, e_none_fail_list))
  e_none_fail_df$tree1 <- "None"
  e_none_fail_df$tree2 <- "Fail"
  e_none_fail_df$test <- test
  e_none_fail_df$dataset <- dataset
  e_none_fail_df$support_value_type <- "UFB"
}  else {
  e_none_fail_df <- data.frame()
}

# Assemble information into a dataframe
e_df <- rbind(e_pass_none_df, e_none_pass_df, e_fail_none_df, e_none_fail_df)
names(e_df) <- c("support_value", "edge_length", "node1", "node2", "tree1", "tree2", "test", "dataset", "support_value_type")
e_df <- e_df[, c("tree1", "tree2", "test", "dataset", "support_value_type", "support_value", "edge_length", "node1", "node2")]


# Write a function that given an edge and a tree, will go into the tree and get for that edge the support value and the branch length
get.edge.details <- function(edge, tree){
  edge_length <- tree$edge.length[edge]
  edge_support <- tree$node.label[edge]
  edge_node1 <- tree$edge[edge,][1]
  edge_node2 <- tree$edge[edge,][2]
  edge_info <- c(edge_support, edge_length, edge_node1, edge_node2)
  return(edge_info)
}

# Write a function that given two trees will return the dataframe of the information about the distinct edges
compare.distinct.edges.of.two.trees <- function(tree_file_1, tree_file_2, tree1_name, tree2_name, test_name, dataset_name, support_value_type_name){
  # Read in trees
  t_1 <- read.tree(tree_file_1)
  t_2 <- read.tree(tree_file_2)
  # Determine which edges are contained in one tree and not the other: returns a numeric vector of edge ids for the first tree
  e_1_2 <- distinct.edges(t_1, t_2)
  e_2_1 <- distinct.edges(t_2, t_1)
  
  # Collect information about each of those different edges
  if (length(e_1_2) > 0){
    e_1_2_list <- lapply(e_1_2, get.edge.details, t_1)
    e_1_2_df <- as.data.frame(do.call(rbind, e_pass_none_list))
    e_1_2_df$tree1 <- tree1_name
    e_1_2_df$tree2 <- tree2_name
    e_1_2_df$test <- test_name
    e_1_2_df$dataset <- dataset_name
    e_1_2_df$support_value_type <- support_value_type_name
  } else {
    e_1_2_df <- data.frame()
  }
  
  if (length(e_none_pass) > 0){
    e_2_1_list <- lapply(e_2_1, get.edge.details, t_2)
    e_2_1_df <- as.data.frame(do.call(rbind, e_none_pass_list))
    e_2_1_df$tree1 <- tree1_name
    e_2_1_df$tree2 <- tree2_name
    e_2_1_df$test <- test_name
    e_2_1_df$dataset <- dataset_name
    e_none_pass_df$support_value_type <- support_value_type_name
  } else {
    e_2_1_df <- data.frame()
  }
  
  
}








