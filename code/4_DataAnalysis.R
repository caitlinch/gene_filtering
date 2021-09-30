### empirical_treelikeness/4_DataAnalysis.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
library(ape)
library(distory)
library(ggplot2)
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

datasets_to_identify_distinct_edges <- c()


##### Step 3: Source function files and prepare variables for analysis #####
print("Source function files")
source(paste0(maindir, "code/func_analysis.R"))

print("Prepare variables for indexing")
# Name vectors for tree roots and alignment locations so they can be accessed via index
names(dataset_tree_roots) <- input_names
names(alignment_dir) <- input_names



##### Step 4: Extract the posterior probabilities/ bootstraps of the trees #####
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
    
    test_df_filename <- paste0(node_output_dir, dataset, "_", test, "_ExtractDistinctEdges.csv")
    
    # If results csv does not exist, calculate results
    if (file.exists(test_df_filename) == FALSE){
      
      # Run IQ-Tree for shallow datasets with completed IQ-Tree species trees
      if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
        print("Compare IQ-Tree trees")
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
        test_df_pass_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, 
                                                               tree1_name = "Pass", tree2_name = "None", test_name = test, 
                                                               dataset_name = dataset, support_value_type_name = "BS")
        test_df_fail_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, 
                                                               tree1_name = "Fail", tree2_name = "None", test_name = test, 
                                                               dataset_name = dataset, support_value_type_name = "BS")
      }
      
      print("Compare ASTRAL trees")
      ## ASTRAL trees: Create dataframe detailing differences in posterior probabilities between the two trees
      # Get the list of trees estimated in ASTRAL for this dataset
      test_trees <- grep(test, all_files, value = TRUE)
      astral_trees <- grep(".tre", grep(".ASTRAL", test_trees, value = TRUE), value = TRUE)
      pass_tree_file <- paste0(dataset_tree_dir, grep("pass", astral_trees, value = TRUE))
      none_tree_file <- paste0(dataset_tree_dir, grep("NoTest", grep(".tre", grep(".ASTRAL", all_files, value = TRUE), value = TRUE), value = TRUE))
      # Create dataframes
      # Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
      test_df_pass_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, tree1_name = "Pass", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
      # Run fail tree for shallow datasets only
      if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
        fail_tree_file <- paste0(dataset_tree_dir, grep("fail", astral_trees, value = TRUE))
        test_df_fail_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, tree1_name = "Fail", tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
      }
      
      print("Saving dataframe")
      # Combine all four dataframes into one
      if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
        test_df <- rbind(test_df_pass_iq, test_df_fail_iq, test_df_pass_astral, test_df_fail_astral)
      } else {
        test_df <- rbind(test_df_pass_astral)
      }
      # Save dataset
      write.csv(test_df, file = test_df_filename)
    }
  }
}


# Collate all dataframes
node_df_filename <- paste0(node_output_dir, "Collated_ExtractDistinctEdges.csv")
if (file.exists(node_df_filename) == FALSE){
  all_csvs <- list.files(node_output_dir)
  all_csvs <- grep("Collated", all_csvs, value = TRUE, invert = TRUE)
  all_csvs <- grep("plots", all_csvs, value = TRUE, invert = TRUE)
  all_csvs <- paste0(node_output_dir, all_csvs)
  all_csv_dfs <- lapply(all_csvs, read.csv)
  node_df <- as.data.frame(do.call(rbind, all_csv_dfs))
  write.csv(node_df, file = node_df_filename)
} else {
  node_df <- read.csv(node_df_filename)
}



##### Step 5: Compare the posterior probabilities/ bootstraps of the trees #####
#### Create new folder for plots ####
if (dir.exists(paste0(node_output_dir, "plots/")) == FALSE){
  dir.create(paste0(node_output_dir, "plots/"))
}

#### Create new columns to facilitate plotting ####
# Create a new factored dataset columns for nice plotting
node_df$dataset_fac <- factor(node_df$dataset, levels = c("Vanderpool2020", "Pease2016", "Strassert2021", "1KP"), 
                              labels = c("Vanderpool et al. (2020)", "Pease et al. (2016)", "Strassert et al. (2020)", "1000 Plants"), 
                              ordered = TRUE)
node_df$edge_type_fac <- factor(node_df$edge_type, levels = c("Congruent", "Conflicting"), ordered = TRUE)
node_df$test_fac <- factor(node_df$test, levels = c("PHI", "maxchi", "geneconv", "allTests"), labels = c("PHI", "MaxChi", "GENECONV", "All tests"), ordered = TRUE)

# Create a new comparison tree column 
pass_inds <- which(node_df$tree1 == "Pass" | node_df$tree2 == "Pass")
fail_inds <- which(node_df$tree1 == "Fail" | node_df$tree2 == "Fail")
node_df$comparison_tree <- NA
node_df$comparison_tree[pass_inds] <- "Pass"
node_df$comparison_tree[fail_inds] <- "Fail"
# Create a column combining the comparison tree and the edge type columns
node_df$boxplot_groups <- paste0(node_df$comparison_tree, ",\n ", node_df$edge_type_fac)
node_df$boxplot_fac <- factor(node_df$boxplot_groups, levels = c("Pass,\n Congruent", "Pass,\n Conflicting", "Fail,\n Congruent", "Fail,\n Conflicting"), ordered = TRUE)

#### Plot the shallow datasets ####
# Separate out the two shallow datasets
shallow_df <- node_df[(node_df$dataset == "Vanderpool2020" | node_df$dataset == "Pease2016"),]

# All trees estimated from loci that pass test have 0 conflicting branches. Add dummy data so that label appears in the plot
dummy_df <- data.frame(boxplot_fac = rep("Pass,\n Conflicting", 16), test_fac = rep(c("PHI", "MaxChi", "GENECONV", "All tests"),4),
                       dataset_fac = c(rep("Vanderpool et al. (2020)",8), rep("Pease et al. (2016)",8)), 
                       edge_length = rep(NA,16), support_value = rep(NA,16), support_value_type = c(rep("BS",4), rep("PP",4), rep("BS",4), rep("PP",4)))
# Make sure bs_df is formatted
shallow_bs_df <- bs_df[, c("boxplot_fac", "test_fac", "dataset_fac", "edge_length", "support_value" , "support_value_type")]
# Combine two datasets
shallow_bs_df <- rbind(shallow_bs_df, dummy_df)

# Separate into ASTRAL and IQ-Tree data frames
pp_df <- shallow_df[(shallow_df$support_value_type == "PP"),]
bs_df <- shallow_df[(shallow_df$support_value_type == "BS"),]

# Make a nice plot of the shallow datasets, faceted by test and dataset
# For ASTRAL trees (with posterior probabilities)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Posterior Probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05), limits = c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
      axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
      axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20),
      strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_ASTRAL_posteriorProbability_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_ASTRAL_posteriorProbability_conflicting_branches.png"), plot = p, device = "png")

# For ASTRAL trees (with branch lengths)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Branch length (coalescent units)", breaks = seq(0,6,2),  labels = seq(0,6,2), minor_breaks = seq(0,7,0.5), limits = c(0,7)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_ASTRAL_branchLength_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_ASTRAL_branchLength_conflicting_branches.png"), plot = p, device = "png")

# For IQ-Tree trees (with ultrafast bootstrap suport values)
p = ggplot(data = shallow_bs_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Ultrafast bootstrap support value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_IQTree_UFBootstraps_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_IQTree_UFBootstraps_conflicting_branches.png"), plot = p, device = "png")

# For IQ-Tree trees (with branch lengths)
p = ggplot(data = shallow_bs_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac, scales = "free_y") + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Branch length (substitutions per site)") +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_IQTree_branchLength_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/ShallowDatasets_IQTree_branchLength_conflicting_branches.png"), plot = p, device = "png")

#### Plot the deep datasets ####
# Separate out the two shallow datasets
deep_df <- node_df[(node_df$dataset == "Strassert2021" | node_df$dataset == "1KP"),]

# Separate into ASTRAL and IQ-Tree data frames
pp_df <- deep_df[(deep_df$support_value_type == "PP"),]
bs_df <- deep_df[(deep_df$support_value_type == "BS"),]

# Make a nice plot of the deep datasets, faceted by test and dataset
# For ASTRAL trees (with posterior probabilities)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Posterior Probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05), limits = c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22),
        axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/DeepDatasets_ASTRAL_posteriorProbability_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/DeepDatasets_ASTRAL_posteriorProbability_conflicting_branches.png"), plot = p, device = "png")

# For ASTRAL trees (with branch lengths)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Branch length (coalescent units)", breaks = seq(0,6,2),  labels = seq(0,6,2), minor_breaks = seq(0,6,0.5), limits = c(0,6)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 25), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 22),
        axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 25))
ggsave(filename = paste0(node_output_dir, "plots/DeepDatasets_ASTRAL_branchLength_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/DeepDatasets_ASTRAL_branchLength_conflicting_branches.png"), plot = p, device = "png")


#### For whole dataset ####
node_df2 <- node_df
node_df2$dataset_fac <- factor(node_df2$dataset, levels = c("Vanderpool2020", "Pease2016", "Strassert2021", "1KP"), 
                              labels = c("Vanderpool2020", "Pease2016", "Strassert2020", "1000 Plants"), 
                              ordered = TRUE)

pp_df <- node_df2[(node_df2$support_value_type == "PP"),]
bs_df <- node_df2[(node_df2$support_value_type == "BS"),]

# Make a nice plot of the deep datasets, faceted by test and dataset
# For ASTRAL trees (with posterior probabilities)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Posterior Probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05), limits = c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 22), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title.y = element_text(size = 22), axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 20))
ggsave(filename = paste0(node_output_dir, "plots/AllDatasets_ASTRAL_posteriorProbability_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/AllDatasets_ASTRAL_posteriorProbability_conflicting_branches.png"), plot = p, device = "png")

# For ASTRAL trees (with branch lengths)
p = ggplot(data = pp_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
  facet_grid(dataset_fac ~ test_fac) + 
  theme_bw() +
  scale_x_discrete(name = "Tree and edge type") +
  scale_y_continuous(name = "Branch length (coalescent units)", breaks = seq(0,6,2),  labels = seq(0,6,2), minor_breaks = seq(0,6,0.5), limits = c(0,6)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title.x = element_text(size = 22), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
        axis.title.y = element_text(size = 22), axis.text.y = element_text(size = 18),
        strip.text = element_text(size = 20))
ggsave(filename = paste0(node_output_dir, "plots/AllDatasets_ASTRAL_branchLength_conflicting_branches.pdf"), plot = p, device = "pdf")
ggsave(filename = paste0(node_output_dir, "plots/AllDatasets_ASTRAL_branchLength_conflicting_branches.png"), plot = p, device = "png")




