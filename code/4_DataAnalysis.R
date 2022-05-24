### gene_filtering/4_DataAnalysis.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Caitlin Cherryh 2022

## This script:
# 1. Creates a variety of plots and figures



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir             <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# tree_data_dir       <- Location of the gene trees
# test_data_dir       <- Location of the results from the AU test and QuartetNetwork Goodness of Fit tests
# output_dir          <- for saving collated output and results from treelikeness analysis.
# primate_data_dir    <- directory containing alignments for individual loci from the Vanderpool et. al. (2020) Primates dataset
# cebidae_trees       <- text file containing the three possible topologies for the Cebidae clade in the Primates dataset
# comparison_trees    <- text file containing the three possible topologies around a deep split within the Primates dataset

# iqtree_path         <- location of IQ-Tree2 executable

# input_names                         <- set name(s) for the dataset(s)
# dataset_tree_roots                  <- set which taxa is outgroup for each dataset
# tests_to_run                        <- a list, with a vector for each dataset specifying which of the recombination detection methods should be tested 
#                                           Options: "allTests", "PHI", "maxchi" and "geneconv"

# datasets_to_identify_distinct_edges <- which datasets to run the distinct.edges function on (to investigate the branches that appear in one tree but not the other)
#                                           To run all, set to = input_names OR to run none, set to = c()
# plotting                            <- whether to plot figures (TRUE = yes, FALSE = no)
# check_primate_loci                  <- whether to apply the AU test to each loci within the Primates dataset (TRUE = yes, FALSE = no)
# plot_primate_loci                   <- whether to plot results of the AU test from loci within the Primates dataset (TRUE = yes, FALSE = no)

### Caitlin's paths ###
location = "server"
if (location == "local"){
  maindir             <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
  tree_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
  test_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
  output_dir          <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"
  primate_data_dir    <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/"
  cebidae_trees       <- paste0(maindir, "primate_tree_topologies/Cebidae_three_possible_topologies.txt")
  comparison_trees    <- paste0(maindir, "primate_tree_topologies/ComparisonTrees_three_possible_topologies.txt")
  
  iqtree_path       <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree"
  
  num_threads       <- 1
  
} else if (location == "server"){
  maindir             <- "/data/caitlin/empirical_treelikeness/"
  tree_data_dir       <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  test_data_dir       <- "/data/caitlin/empirical_treelikeness/Output_dataAnalysis/"
  output_dir          <- "/data/caitlin/empirical_treelikeness/Output/"
  primate_data_dir    <- "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/"
  cebidae_trees       <- paste0(maindir, "primate_tree_topologies/Cebidae_three_possible_topologies.txt")
  comparison_trees    <- paste0(maindir, "primate_tree_topologies/ComparisonTrees_three_possible_topologies.txt")
  
  iqtree_path       <- "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree"
  
  num_threads       <- 20
}

input_names <- c("Vanderpool2020", "Pease2016", "Whelan2017", "1KP")
dataset_tree_roots <- list("1KP" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                     "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                     "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                           "Whelan2017" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                           "Vanderpool2020" = c("Mus_musculus"), 
                           "Pease2016" = c("LA4116", "LA2951", "LA4126"))
tests_to_run <- list("Vanderpool2020" = c("PHI", "maxchi", "geneconv", "allTests"),
                     "Pease2016" = c("PHI", "maxchi", "geneconv", "allTests"),
                     "Whelan2017" = c("PHI", "maxchi", "geneconv"),
                     "1KP" = c("PHI", "maxchi"))

datasets_to_identify_distinct_edges <- c()
plotting = FALSE
check_primate_loci = TRUE
plot_primate_loci = FALSE
### End of Caitlin's paths ###



##### Step 2: Open packages and set directories #####
# Open packages
if ( (length(datasets_to_identify_distinct_edges) > 0) | (plotting == TRUE) ){
  library(ape)
  library(distory)
  library(ggplot2)
  library(patchwork)
}
library(parallel)


##### Step 3: Source function files and prepare variables for analysis #####
source(paste0(maindir, "code/func_analysis.R"))



##### Step 4: Extract the posterior probabilities/ bootstraps of the trees #####
node_output_dir <- paste0(output_dir, "node_comparisons/")
if (dir.exists(node_output_dir) == FALSE){
  dir.create(node_output_dir)
}

for (dataset in datasets_to_identify_distinct_edges){
  # Identify file containing species trees
  dataset_tree_dir <- paste0(tree_data_dir, dataset, "/species_trees/")
  # Collect all files
  dataset_tree_dir_files <- list.files(dataset_tree_dir, recursive = TRUE)
  # Remove any trees or files from old GENECONV/All tests runs
  all_files <- grep("old_geneconv|Old_geneconv|Old_Geneconv|00_|zz_", dataset_tree_dir_files, invert = TRUE, value = TRUE)
  
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
        # Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
        test_df_pass_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, 
                                                               tree1_name = "Pass", tree2_name = "None", test_name = test, 
                                                               dataset_name = dataset, support_value_type_name = "BS")
        test_df_fail_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, 
                                                               tree1_name = "Fail", tree2_name = "None", test_name = test, 
                                                               dataset_name = dataset, support_value_type_name = "BS")
      } else if (dataset == "Whelan2017" | dataset == "1KP"){
        ## ML trees: Create dataframe detailing differences in posterior probabilities between the two trees
        # If dataset is 1KP, collect RAxML trees estimated with no free rate models
        # For other datasets, collect IQ-Tree trees (no restrictions on models)
        if (dataset == "1KP"){
          print("Compare ML trees")
          # Get the list of trees estimated in RAxML-NG for this dataset
          # Get the list of trees estimated in IQ-Tree for this dataset
          raxml_trees <- grep("bestTree", all_files, value = TRUE)
          raxml_trees <- grep("noFreeRates", raxml_trees, value = TRUE)
          test_trees <- grep(test, raxml_trees, value = TRUE)
          # Make the full filepaths for each of the three trees (test pass, test fail, and no test)
          none_tree_file <- paste0(dataset_tree_dir, grep("NoTest", raxml_trees, value = TRUE))
          pass_tree_file <- paste0(dataset_tree_dir, grep("pass", test_trees, value = TRUE))
        } else {
          print("Compare IQ-Tree trees")
          # Get the list of trees estimated in IQ-Tree for this dataset
          test_trees <- grep(test, all_files, value = TRUE)
          iq_trees <- grep(".contree", test_trees, value = TRUE)
          # Make the full filepaths for each of the three trees (test pass, test fail, and no test)
          none_tree_file <- paste0(dataset_tree_dir, grep(".contree", grep("NoTest", all_files, value = TRUE), value = TRUE))
          pass_tree_file <- paste0(dataset_tree_dir, grep("pass", iq_trees, value = TRUE))
        }
        # Want to collect the information about splits present in one tree but not the other (i.e. in T_all,pass vs T_None)
        test_df_pass_iq <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, 
                                                               tree1_name = "Pass", tree2_name = "None", test_name = test, 
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
      test_df_pass_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = pass_tree_file, tree_file_2 = none_tree_file, tree1_name = "Pass", 
                                                                 tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
      # Run fail tree for shallow datasets only
      if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
        fail_tree_file <- paste0(dataset_tree_dir, grep("fail", astral_trees, value = TRUE))
        test_df_fail_astral <- compare.distinct.edges.of.two.trees(tree_file_1 = fail_tree_file, tree_file_2 = none_tree_file, tree1_name = "Fail", 
                                                                   tree2_name = "None", test_name = test, dataset_name = dataset, support_value_type_name = "PP")
      }
      
      print("Saving dataframe")
      # Combine all four dataframes into one
      if (dataset == "Vanderpool2020" | dataset == "Pease2016"){
        test_df <- rbind(test_df_pass_iq, test_df_fail_iq, test_df_pass_astral, test_df_fail_astral)
      } else {
        test_df <- rbind(test_df_pass_iq, test_df_pass_astral)
      }
      # Save dataset
      write.csv(test_df, file = test_df_filename)
    }
  }
}



##### Step 5: Compare the posterior probabilities/ bootstraps of the trees #####
if (plotting == TRUE){
  #### Open the .csv file containing branch lengths/support for all analyses for all four datasets
  node_df_filename <- paste0(node_output_dir, "04_AllDatasets_Collated_ExtractDistinctEdges.csv")
  if (file.exists(node_df_filename) == FALSE){
    # Collate dataframe of conflicting/congruent branches from all four datasets
    all_csvs <- list.files(node_output_dir)
    all_csvs <- grep("\\.csv", all_csvs, value = TRUE)
    all_csvs <- grep("Collated", all_csvs, value = TRUE, invert = TRUE)
    all_csvs <- paste0(node_output_dir, all_csvs)
    all_csv_dfs <- lapply(all_csvs, read.csv)
    node_df <- as.data.frame(do.call(rbind, all_csv_dfs))
    write.csv(node_df, file = node_df_filename)
  } else {
    # Open collated dataframe
    node_df <- read.csv(node_df_filename)
  }
  
  #### Create new folder for plots ####
  plot_dir <- paste0(output_dir, "plots/")
  if (dir.exists(plot_dir) == FALSE){
    dir.create(plot_dir)
  }
  
  #### Create new columns to facilitate plotting ####
  # Create a new factored dataset columns for nice plotting
  node_df$dataset_fac <- factor(node_df$dataset, levels = c("Vanderpool2020", "Pease2016", "Whelan2017", "1KP"), 
                                labels = c("Primates", "Tomatoes", "Metazoans", "Plants"), 
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
  
  # All trees estimated from Pease2016/Vanderpool2020 loci that pass tests have 0 conflicting branches. 
  # Add dummy data to add an empty boxplot for the conflicting branches to the plot 
  # This makes it clear that there were no conflicting branches (as opposed to excluding the conflicting branches)
  dummy_df <- data.frame(boxplot_fac = rep("Pass,\n Conflicting", 16), test_fac = rep(c("PHI", "MaxChi", "GENECONV", "All tests"),4),
                         dataset_fac = c(rep("Primates",8), rep("Tomatoes",8)), 
                         edge_length = rep(NA,16), support_value = rep(NA,16), support_value_type = c(rep("BS",4), rep("PP",4), rep("BS",4), rep("PP",4)))
  collate_df <- node_df[, c("boxplot_fac", "test_fac", "dataset_fac", "edge_length", "support_value" , "support_value_type")]
  collate_df <- rbind(collate_df, dummy_df)
  
  # Separate into concatenated/coalescent data by support value type (posterior probability or bootstrap)
  pp_df <- collate_df[(collate_df$support_value_type == "PP"), ]
  bs_df <- collate_df[(collate_df$support_value_type == "BS"), ]
  
  
  #### Create a lovely plot of posterior probability support values for the ASTRAL trees
  # Break pp_df into three sections to plot: one for Tomatoes/Primates, one for Metazoans, and one for Plants
  pp1_df <- pp_df[((pp_df$dataset_fac == "Primates") | (pp_df$dataset_fac == "Tomatoes")),]
  pp2_df <- pp_df[(pp_df$dataset_fac == "Metazoans"),]
  pp3_df <- pp_df[(pp_df$dataset_fac == "Plants"),]
  
  # Plot each of the three sections, faceted by test and dataset
  p1 = ggplot(data = pp1_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Posterior probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
    labs(title = "a.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p2 <- ggplot(data = pp2_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Posterior probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
    labs(title = "b.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p3 <- ggplot(data = pp3_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Posterior probability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
    labs(title = "c.") +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  # Combine each plot into one big plot
  p = p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(8, 4, 4))
  ggsave(filename = paste0(plot_dir, "ASTRAL_posteriorProbability_conflicting_branches.pdf"), plot = p, device = "pdf", units = "in", width = 8, height = 10)
  
  
  #### Create a lovely plot of posterior probability support values for the ASTRAL trees
  # Break pp_df into three sections to plot: one for Tomatoes/Primates, one for Metazoans, and one for Plants
  pp1_df <- pp_df[((pp_df$dataset_fac == "Primates") | (pp_df$dataset_fac == "Tomatoes")),]
  pp2_df <- pp_df[(pp_df$dataset_fac == "Metazoans"),]
  pp3_df <- pp_df[(pp_df$dataset_fac == "Plants"),]
  
  # Plot each of the three sections, faceted by test and dataset
  p1 = ggplot(data = pp1_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,6,2),  labels = seq(0,6,2), minor_breaks = seq(0,7,1), limits = c(0,7)) +
    labs(title = "a.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p2 <- ggplot(data = pp2_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,4,2),  labels = seq(0,4,2), minor_breaks = seq(0,4,0.5), limits = c(0,4)) +
    labs(title = "b.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p3 <- ggplot(data = pp3_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,6,2),  labels = seq(0,6,2), minor_breaks = seq(0,6,1), limits = c(0,6)) +
    labs(title = "c.") +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  # Combine each plot into one big plot
  p = p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(8, 4, 4))
  ggsave(filename = paste0(plot_dir, "ASTRAL_edgeLength_conflicting_branches.pdf"), plot = p, device = "pdf", units = "in", width = 8, height = 10)
  
  
  #### Create a lovely plot of ultrafast bootstrap support values for the maximum likelihood trees
  # Break pp_df into three sections to plot: one for Tomatoes/Primates, one for Metazoans, and one for Plants
  bs1_df <- bs_df[((bs_df$dataset_fac == "Primates") | (bs_df$dataset_fac == "Tomatoes")),]
  bs2_df <- bs_df[(bs_df$dataset_fac == "Metazoans"),]
  
  # Plot each of the three sections, faceted by test and dataset
  p1 = ggplot(data = bs1_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Ultrafast bootstrap support value", breaks = seq(0,100,20),  labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
    labs(title = "a.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p2 <- ggplot(data = bs2_df, aes(x = boxplot_fac, y = support_value)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Ultrafast bootstrap support value", breaks = seq(0,100,20),  labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
    labs(title = "b.") +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  # Combine each plot into one big plot
  p = p1 + p2 + plot_layout(ncol = 1, heights = c(8, 4))
  ggsave(filename = paste0(plot_dir, "ML_ultrafastBootstrapSupport_conflicting_branches.pdf"), plot = p, device = "pdf", units = "in", width = 8, height = 8)
  
  
  #### Create a lovely plot of edge lengths for the maximum likelihood trees
  # Break pp_df into three sections to plot: one for Tomatoes/Primates, one for Metazoans, and one for Plants
  bs1_df <- bs_df[((bs_df$dataset_fac == "Primates") | (bs_df$dataset_fac == "Tomatoes")),]
  bs2_df <- bs_df[(bs_df$dataset_fac == "Metazoans"),]
  bs3_df <- bs_df[(bs_df$dataset_fac == "Plants"),]
  
  # Plot each of the three sections, faceted by test and dataset
  p1 = ggplot(data = bs1_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,0.06,0.02),  labels = seq(0,0.06,0.02), minor_breaks = seq(0,0.06,0.01), limits = c(0,0.06)) +
    labs(title = "a.") +
    theme(plot.title = element_text(size = 20), 
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p2 <- ggplot(data = bs2_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,0.6,0.2),  labels = seq(0,0.6,0.2), minor_breaks = seq(0,0.6,0.1), limits = c(0,0.6)) +
    labs(title = "b.") +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  p3 <- ggplot(data = bs3_df, aes(x = boxplot_fac, y = edge_length)) + geom_boxplot() +
    facet_grid(dataset_fac ~ test_fac) + 
    theme_bw() +
    scale_x_discrete(name = "Tree and edge type") +
    scale_y_continuous(name = "Edge length", breaks = seq(0,1.2,0.4),  labels = seq(0,1.2,0.4), minor_breaks = seq(0,1.2,0.2), limits = c(0,1.2)) +
    labs(title = "c.") +
    theme(plot.title = element_text(size = 20),
          axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 15))
  
  # Combine each plot into one big plot
  p = p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(8, 4, 4))
  ggsave(filename = paste0(plot_dir, "ML_edgeLength_conflicting_branches.pdf"), plot = p, device = "pdf", units = "in", width = 8, height = 10)
}



#### Step 6: Investigate topology of primates dataset ####
if (check_primate_loci == TRUE){
  ## Make a new folder to save all the output files
  check_dir <- paste0(output_dir, "check_primates/")
  if (dir.exists(check_dir) == FALSE){ dir.create(check_dir) }
  
  ## Assemble a list of all loci names
  all_loci <- gsub("\\.fa", "", list.files(primate_data_dir))
  
  ## Set the details for each run
  AU_test_ids <- c("Cebidae", "Comparison")
  AU_test_details <- list("Cebidae" = c(name = "Cebidae", tree_path = cebidae_trees, output_directory = paste0(check_dir, "AU_cebidae_trees/")),
                          "Comparison" = c(name = "Comparison", tree_path = comparison_trees, output_directory = paste0(check_dir, "AU_comparison_trees/")) )
  
  # Have to run the AU test for each loci twice: one for the three trees in cebidae_trees, and once for the three trees in comparison_trees
  # Now, iterate through the datasets
  for (id in AU_test_ids){
    ## Set folder for this set of AU test runs
    id_dir <- AU_test_details[[id]][["output_directory"]]
    if (dir.exists(id_dir) == FALSE){ dir.create(id_dir) }
    
    ## Set folders and parameters for running AU test
    # Set output folder
    output_folder <- paste0(AU_test_details[[id]][["output_directory"]], "output/")
    if (dir.exists(output_folder) == FALSE){ dir.create(output_folder) }
    # Set folder for saving CSV files
    csv_folder <- paste0(AU_test_details[[id]][["output_directory"]], "csvs/")
    if (dir.exists(csv_folder) == FALSE){ dir.create(csv_folder) }
    # Set tree topologies to compare
    three_trees_path <- AU_test_details[[id]][["tree_path"]]
    
    ## Perform AU test
    mclapply(all_loci, perform.AU.test, primate_data_dir, output_folder, csv_folder, three_trees_path, iqtree_path, trim.taxa = TRUE, mc.cores = num_threads)
    
    ## Read in all the csv files for this AU_test_id and combine them
    all_AU_csvs <- paste0(csv_folder,list.files(csv_folder))
    all_csvs <- lapply(all_AU_csvs, read.csv, stringsAsFactors = FALSE, row.names = 1)
    AU_df <- as.data.frame(do.call(rbind, all_csvs))
    AU_df_name <- paste0(check_dir, "Primates_", id,"_AU_test_collated.csv")
    write.csv(AU_df, file = AU_df_name)
  }
}

if (plot_primate_loci == TRUE){
  ## Make a new folder to save all the output files
  check_plots_dir <- paste0(output_dir, "check_primates/plots/")
  if (dir.exists(check_plots_dir) == FALSE){ dir.create(check_plots_dir) }
}










