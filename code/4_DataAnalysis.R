### gene_filtering/code/4_DataAnalysis.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Caitlin Cherryh 2023

## This script:
# 1. Extracts information about the branches within the subset trees for each dataset
# 2. Creates a variety of plots and figures



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                   <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# tree_data_dir             <- Location of the gene trees
# test_data_dir             <- Location of the results from the AU test and QuartetNetwork Goodness of Fit tests
# output_dir                <- for saving collated output and results from treelikeness analysis.
# plot_dir                  <- for saving outlier branch plots
# primate_data_dir          <- directory containing alignments for individual loci from the Vanderpool et. al. (2020) Primates dataset
# cebidae_trees             <- text file containing the three possible topologies for the Cebidae clade in the Primates dataset
# comparison_trees          <- text file containing the three possible topologies around a deep split within the Primates dataset
# iqtree_path               <- location of IQ-Tree2 executable
# recombination_output_file <- location of results of recombination detection tests
# dataset_tree_roots        <- set which taxa is outgroup for each dataset
# check_primate_loci        <- whether to apply the AU test to each loci within the Primates dataset (TRUE = yes, FALSE = no)
# plot_primate_loci         <- whether to plot results of the AU test from loci within the Primates dataset (TRUE = yes, FALSE = no)

### Caitlin's paths ###
location = "local"
if (location == "local"){
  maindir             <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
  tree_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/Species_trees_for_publication/"
  test_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
  output_dir          <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"
  plot_dir            <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
  primate_data_dir    <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/"
  
  iqtree_path       <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.0-rc1-MacOSX/bin/iqtree"
  
  recombination_output_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/01_AllDatasets_RecombinationDetection_complete_collated_results.csv"
  
  num_threads       <- 1
  
} else if (location == "server"){
  maindir             <- "/data/caitlin/empirical_treelikeness/"
  tree_data_dir       <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  test_data_dir       <- "/data/caitlin/empirical_treelikeness/Output_dataAnalysis/"
  output_dir          <- "/data/caitlin/empirical_treelikeness/Output/"
  plot_dir            <- output_dir
  primate_data_dir    <- "/data/caitlin/empirical_treelikeness/Data_Vanderpool2020/"
  
  iqtree_path       <- "/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree"
  
  num_threads       <- 20
}

cebidae_trees       <- paste0(maindir, "primate_tree_topologies/Cebidae_three_possible_topologies.txt")
comparison_trees    <- paste0(maindir, "primate_tree_topologies/ComparisonTrees_three_possible_topologies.txt")

dataset_tree_roots <- list("1KP" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                     "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                     "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                           "Whelan2017" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                           "Vanderpool2020" = c("Mus_musculus"), 
                           "Pease2016" = c("LA4116", "LA2951", "LA4126"))

check_primate_loci = FALSE
plot_primate_loci = FALSE
### End of Caitlin's paths ###



##### Step 2: Open packages and set directories #####
# Open packages
if ( (length(datasets_to_identify_distinct_edges) > 0) | identify_outlier_edges == TRUE){
  library(ape)
  library(distory)
} 
if (plot_distinct_edges == TRUE | identify_outlier_edges == TRUE){
  library(ggplot2)
  library(ggtree)
  library(patchwork)
}
library(parallel)



##### Step 3: Source function files and prepare variables for analysis #####
source(paste0(maindir, "code/func_analysis.R"))
source(paste0(maindir, "code/func_comparison.R"))
source(paste0(maindir, "code/func_plots.R"))



##### Step 4: Extract the posterior probabilities/ bootstraps of the trees #####
# Create output file names
tree_comp_df_file <- paste0(output_dir, "BranchSupport_input_trees.csv")
branch_support_df_file <- paste0(output_dir, "BranchSupport_values.csv")

# Open input dataframe
if (file.exists(tree_comp_df_file)){
  tree_comp_df <- read.csv(tree_comp_df_file, stringsAsFactors = FALSE)
} else {
  # Extract all tree files
  all_files <- list.files(tree_data_dir)
  tree_files <- grep(".tre|.treefile", all_files, value = T)
  # Create dataframe for testing every pair of trees
  tree_comp_df      <- data.frame(dataset = c(rep("Tomatoes", 8),
                                              rep("Primates",  8),
                                              rep("Metazoan", 3),
                                              rep("Plants", 2),
                                              rep("Tomatoes", 8),
                                              rep("Primates",  8),
                                              rep("Metazoan", 3),
                                              rep("Plants", 2)),
                                  clean_tree = paste0(rep(c("Tomatoes_allTests_pass", "Tomatoes_allTests_pass",
                                                            "Tomatoes_geneconv_pass", "Tomatoes_geneconv_pass",
                                                            "Tomatoes_maxchi_pass", "Tomatoes_maxchi_pass",
                                                            "Tomatoes_PHI_pass", "Tomatoes_PHI_pass",
                                                            "Primates_allTests_pass", "Primates_allTests_pass",
                                                            "Primates_geneconv_pass", "Primates_geneconv_pass",
                                                            "Primates_maxchi_pass", "Primates_maxchi_pass",
                                                            "Primates_PHI_pass", "Primates_PHI_pass",
                                                            "Metazoan_geneconv_pass", "Metazoan_maxchi_pass",
                                                            "Metazoan_PHI_pass",
                                                            "Plants_maxchi_pass", "Plants_PHI_pass"), 2), 
                                                      c(rep("_ASTRAL_species.tre", 21), rep("_CONCAT_IQTREE.treefile", 21))), 
                                  comparison_tree = paste0(rep(c("Tomatoes_allTests_fail",  "Tomatoes_NoTest",
                                                                 "Tomatoes_geneconv_fail",  "Tomatoes_NoTest",
                                                                 "Tomatoes_maxchi_fail",  "Tomatoes_NoTest",
                                                                 "Tomatoes_PHI_fail",  "Tomatoes_NoTest",
                                                                 "Primates_allTests_fail", "Primates_NoTest",
                                                                 "Primates_geneconv_fail", "Primates_NoTest",
                                                                 "Primates_maxchi_fail", "Primates_NoTest",
                                                                 "Primates_PHI_fail" , "Primates_NoTest",
                                                                 "Metazoan_NoTest", "Metazoan_NoTest",
                                                                 "Metazoan_NoTest",
                                                                 "Plants_NoTest", "Plants_NoTest"), 2), 
                                                           c(rep("_ASTRAL_species.tre", 21), rep("_CONCAT_IQTREE.treefile", 21))), 
                                  tree_directory = tree_data_dir,
                                  comparison_id = rep(c("allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                        "allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                        "noTest", "noTest", "noTest",
                                                        "noTest", "noTest"), 2),
                                  analysis_method = c(rep("ASTRAL", 21), rep("IQTREE", 21)) )
  tree_comp_df$clean_id <- unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){paste0(x[[2]], "_", x[[3]])}))
  tree_comp_df$recombination_test <-  unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){x[[2]]}))
  tree_comp_df$comparison_gene_status <- unlist(lapply(strsplit(tree_comp_df$comparison_tree, "_"), function(x){x[[3]]}))
  tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "ASTRAL")] <- "unfiltered"
  tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "CONCAT")] <- "unfiltered"
  tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "fail")] <- "recombinant"
  # Remove PLANT rows for CONCAT method (no support values) plants file names for CONCAT
  
  tree_comp_df$clean_tree[which(tree_comp_df$dataset == "Plants" & tree_comp_df$analysis_method == "IQTREE")] <- gsub("_CONCAT_IQTREE.treefile", "_CONCAT_RAxML_noFreeRates.raxml.bestTree",
                                                                                                                      tree_comp_df$clean_tree[which(tree_comp_df$dataset == "Plants" & tree_comp_df$analysis_method == "IQTREE")])
  tree_comp_df$comparison_tree[which(tree_comp_df$dataset == "Plants" & tree_comp_df$analysis_method == "IQTREE")] <- gsub("_CONCAT_IQTREE.treefile", "_CONCAT_RAxML_noFreeRates.raxml.bestTree",
                                                                                                                           tree_comp_df$comparison_tree[which(tree_comp_df$dataset == "Plants" & tree_comp_df$analysis_method == "IQTREE")])
  # Save csv
  write.csv(tree_comp_df, file = tree_comp_df_file, row.names = FALSE)
}

# Open BranchSupport dataframe
if (file.exists(branch_support_df_file)){
  branch_support_df <- read.csv(branch_support_df_file, stringsAsFactors = FALSE)
} else {
  # Extract qCF values
  branch_support_list <- lapply(1:nrow(tree_comp_df), compare.branch.wrapper, df = tree_comp_df)
  branch_support_df <- as.data.frame(do.call(rbind, branch_support_list))
  branch_support_df <- branch_support_df[ , c("dataset", "clean_tree", "comparison_tree", "clean_id", "comparison_id", "recombination_test", "comparison_gene_status",
                                              "analysis_method", "tree", "split_type", "confidence", "weights")]
  # Save csv
  write.csv(branch_support_df, file = branch_support_df_file, row.names = FALSE)
}



##### Step 5: Compare the posterior probabilities of the trees #####
# Separate posterior probabilities
pp_df <- branch_support_df[which(branch_support_df$analysis_method == "ASTRAL"),]
# Change split_type name
pp_df$split_type[which(pp_df$split_type == "Concordant")] <- "Congruent"
# Create a nicely formatted dataset column
pp_df$dataset_formatted <- factor(pp_df$dataset,
                                  levels = c("Tomatoes", "Primates", "Plants", "Metazoan"),
                                  labels = c("Tomatoes", "Primates", "Plants", "Metazoan"),
                                  ordered = T)
# Add a new columns for faceting/grouping
pp_df$gene_tree_formatted <- factor(pp_df$tree,
                                    levels = c("geneconv_pass", "geneconv_fail", "maxchi_pass", "maxchi_fail",
                                               "PHI_pass", "PHI_fail", "allTests_pass", "allTests_fail", 
                                               "noTest"),
                                    labels = c("Clean", "Recombinant", "Clean", "Recombinant",
                                               "Clean", "Recombinant", "Clean", "Recombinant",
                                               "Unfiltered"),
                                    ordered = TRUE)
pp_df$recombination_test_formatted <- factor(pp_df$recombination_test,
                                             levels = c("geneconv", "maxchi", "PHI", "allTests"),
                                             labels = c("GENCONV", "MaxChi", "PHI", "All tests"),
                                             ordered = T)
# Remove any test, fail rows
pp_df <- pp_df[grep("fail", pp_df$comparison_tree, invert = T), ]
# Convert confidence to numeric
pp_df$confidence <- as.numeric(pp_df$confidence)
# Separate into dataframes for the different datasets
tomatoes_pp <- pp_df[which(pp_df$dataset == "Tomatoes"), ]
primates_pp <- pp_df[which(pp_df$dataset == "Primates"), ]
plant_pp   <- pp_df[which(pp_df$dataset == "Plants"), ]
metazoan_pp <- pp_df[which(pp_df$dataset == "Metazoan"), ]
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes
tomatoes_pp_plot <- ggplot(tomatoes_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Posterior\nprobability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_pp_plot <- ggplot(primates_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Posterior\nprobability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_pp_plot <- ggplot(metazoan_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Posterior\nprobability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_pp_plot <- ggplot(plant_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Posterior\nprobability", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "d.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_pp_plot + primates_pp_plot + met_pp_plot + plants_pp_plot + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "PosteriorProbabilities_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "PosteriorProbabilities_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")



##### Step 6: Plot branch lengths for ASTRAL trees #####
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes
tomatoes_bla_plot <- ggplot(tomatoes_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_bla_plot <- ggplot(primates_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,4.5)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_bla_plot <- ggplot(metazoan_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5, 1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,3.5)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_bla_plot <- ggplot(plant_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,5.5)) +
  labs(title = "d.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_bla_plot + primates_bla_plot + met_bla_plot + plants_bla_plot + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "BranchLengths_ASTRAL_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 13)
quilt_png <- paste0(plot_dir, "BranchLengths_ASTRAL_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 13)



##### Step 7: Compare the bootstraps of the trees #####
# Separate posterior probabilities
bs_df <- branch_support_df[which(branch_support_df$analysis_method == "IQTREE"),]
# Remove any splits that are not either "Concordant" or "Conflicting"
bs_df$split_type[which(bs_df$split_type == "Concordant")] <- "Congruent"
# Create a nicely formatted dataset column
bs_df$dataset_formatted <- factor(bs_df$dataset,
                                  levels = c("Tomatoes", "Primates", "Plants", "Metazoan"),
                                  labels = c("Tomatoes", "Primates", "Plants", "Metazoan"),
                                  ordered = T)
# Add a new columns for faceting/grouping
bs_df$gene_tree_formatted <- factor(bs_df$tree,
                                    levels = c("geneconv_pass", "geneconv_fail", "maxchi_pass", "maxchi_fail",
                                               "PHI_pass", "PHI_fail", "allTests_pass", "allTests_fail", 
                                               "noTest"),
                                    labels = c("Clean", "Recombinant", "Clean", "Recombinant",
                                               "Clean", "Recombinant", "Clean", "Recombinant",
                                               "Unfiltered"),
                                    ordered = TRUE)
bs_df$recombination_test_formatted <- factor(bs_df$recombination_test,
                                             levels = c("geneconv", "maxchi", "PHI", "allTests"),
                                             labels = c("GENCONV", "MaxChi", "PHI", "All tests"),
                                             ordered = T)
# Remove any test, fail rows
bs_df <- bs_df[grep("fail", bs_df$comparison_tree, invert = T), ]
# Convert confidence to numeric
bs_df$confidence <- as.numeric(bs_df$confidence)
# Separate into dataframes for the different datasets
shallow_bs <- bs_df[which(bs_df$dataset %in% c("Tomatoes", "Primates")), ]
metazoan_bs <- bs_df[which(bs_df$dataset == "Metazoan"), ]
plant_bs <-  bs_df[which(bs_df$dataset == "Plants"), ]
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes
shallow_bs_plot <- ggplot(shallow_bs, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot()  +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering")+
  scale_y_continuous(name = "UFB value", breaks = seq(0,120,20),  labels = seq(0,120,20), minor_breaks = seq(0,110,10), limits = c(0,110)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_bs_plot <- ggplot(metazoan_bs, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "UFB value", breaks = seq(0,120,20),  labels = seq(0,120,20), minor_breaks = seq(0,110,10), limits = c(0,110)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- shallow_bs_plot + met_bs_plot + plot_layout(ncol = 1, heights = c(2,1))
quilt_pdf <- paste0(plot_dir, "UltrafastBootstrap_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt)
quilt_png <- paste0(plot_dir, "UltrafastBootstrap_quilt.png")
ggsave(filename = quilt_png, plot = quilt)



##### Step 8: Plot branch lengths for CONCAT trees #####
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes and primates
shallow_blc_plot <- ggplot(shallow_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0, 0.04, 0.01),  labels = seq(0, 0.04, 0.01), minor_breaks = seq(0, 0.04, 0.005), limits = c(0,0.04)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_blc_plot <- ggplot(metazoan_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0, 0.5, 0.1),  labels = seq(0, 0.5, 0.1), minor_breaks = seq(0, 0.5, 0.05), limits = c(0,0.5)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_blc_plot <- ggplot(plant_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,2,0.2),  labels = seq(0,2,0.2), minor_breaks = seq(0,2,0.1), limits = c(0,1.2)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- shallow_blc_plot + met_blc_plot + plants_blc_plot + plot_layout(ncol = 1, heights = c(2,1,1))
quilt_pdf <- paste0(plot_dir, "BranchLengths_CONCAT_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12)
quilt_png <- paste0(plot_dir, "BranchLengths_CONCAT_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 12)



##### Step 6: Identify outlier branches (ones with high support) #####
if (identify_outlier_edges == TRUE){
  ## Identify outlier edges
  # Assemble the filename for and open the completed "04_AllDatasets_Collated_ExtractDistinctEdges.csv" file created in Step 5
  node_output_dir <- paste0(output_dir, "node_comparisons/")
  node_df_filename <- paste0(node_output_dir, "04_AllDatasets_Collated_ExtractDistinctEdges.csv")
  node_df <- read.csv(node_df_filename, stringsAsFactors = FALSE)
  node_df <- node_df[, c("test", "dataset", "support_value_type", "tree1", "tree2", "edge_type", "edge_presence", "support_value", "edge_length", "node1", "node2")]
  
  # Create empty dataframe to append information to
  outlier_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 11))
  names(outlier_df) <- c("test", "dataset", "support_value_type", "tree1", "tree2", "edge_type", "edge_presence", "support_value", "edge_length", "node1", "node2")
  
  # Iterate through the datasets one at a time
  for (dataset in input_names){
    print(dataset)
    # Extract rows for only this dataset
    dataset_df <- node_df[node_df$dataset == dataset,]
    
    # For IQ-Tree2 runs, identify conflicting branches with length more than 1 standard deviations above the mean (on a test by test basis)
    tests <- unique(dataset_df$test)
    for (t in tests){
      print(t)
      # Separate out rows from concatenated trees (IQ-Tree runs)
      test_bs_df <- node_df[node_df$dataset == dataset & node_df$support_value_type == "BS" & node_df$test == t,]
      # Calculate mean and standard deviation
      mean_bs_edge_length <- mean(test_bs_df$edge_length)
      sd_bs_edge_length <- sd(test_bs_df$edge_length)
      # Identify outlier branches
      outlier_bs_df <- test_bs_df[test_bs_df$edge_length > (mean_bs_edge_length + sd_bs_edge_length) & test_bs_df$edge_type == "Conflicting", ]
      # If there are any outlier branches, add them to the outlier_df
      if (nrow(outlier_bs_df) > 0){
        outlier_df <- rbind(outlier_df, outlier_bs_df)
      } # End rbind outliers
    } # End iterating through tests for ultrafast bootstrap results
    
    # For ASTRAL-III runs, identify conflicting branches with length more than 1 standard deviations above the mean (on a test by test basis)
    for (t in tests){
      print(t)
      # Separate out rows from coalescent trees (ASTRAL-III runs)
      test_pp_df <- node_df[node_df$dataset == dataset & node_df$support_value_type == "PP" & node_df$test == t,]
      # Calculate mean and standard deviation
      mean_pp_edge_length <- mean(test_pp_df$edge_length)
      sd_pp_edge_length <- sd(test_pp_df$edge_length)
      # Identify outlier branches
      outlier_pp_df <- test_pp_df[test_pp_df$edge_length > (mean_pp_edge_length + sd_pp_edge_length) & test_pp_df$edge_type == "Conflicting", ]
      # If there are any outlier branches, add them to the outlier_df
      if (nrow(outlier_pp_df) > 0){
        outlier_df <- rbind(outlier_df, outlier_pp_df)
      } # End rbind outliers
    } # End iterating through tests for posterior probability results
    
  } # End iterating through dataset
  
  # Save outlier branches in new csv
  outlier_df_filename <- paste0(node_output_dir, "04_AllDatasets_Collated_ExtractDistinctEdges_OnlyOutliers.csv")
  write.csv(outlier_df, file = outlier_df_filename, row.names = FALSE)
  write.csv(outlier_df, file = paste0(node_output_dir, "outlier_branch_dataframe.csv"))
  
  ## Plot outlier edges
  # Create new folder to plot in
  outlier_plot_dir <- paste0(node_output_dir, "outlier_plots/")
  if (dir.exists(outlier_plot_dir) == FALSE){dir.create(outlier_plot_dir)}
  
  # Get list of trees
  all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
  
  # Identify each outlier branch within the tree
  for (i in 1:nrow(outlier_df)){
    # Get the row for this branch
    row <- outlier_df[i, ]
    # Identify and open the correct tree
    tree_file <- determine.outlier.tree.file(i, outlier_df, all_trees)
    tree <- read.tree(tree_file)
    # Subset the tree at the nodes on either end of the branch
    big_clade <- extract.clade(tree, row$node1)
    little_clade <- extract.clade(tree, row$node2)
    # Assemble file names for plotting
    plot_id <- determine.outlier.plot.name(i, outlier_df)
    big_clade_plot_path <- paste0(outlier_plot_dir, sprintf("%03d_", i), plot_id, "node", row$node1, "_bigClade_OutlierBranch_plot")
    little_clade_plot_path <- paste0(outlier_plot_dir, sprintf("%03d_", i), plot_id, "node", row$node2, "_littleClade_OutlierBranch_plot")
    # Save big clade plot
    pdf(file = paste0(big_clade_plot_path, ".pdf"))
    plot.phylo(big_clade)
    dev.off()
    # Save little clade plot
    pdf(file = paste0(little_clade_plot_path, ".pdf"))
    plot(little_clade)
    dev.off()
  }
  
} # End identify_outlier_edges code

## Pretty plotting for Outlier Branch 1
# Get list of trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Open the outlier dataframe
outlier_df <- read.csv(outlier_df_filename)

# Identify tree files
notest_tree_file <- grep("ASTRAL", grep("NoTest", grep("Metazoan", all_trees, value = TRUE), value = TRUE), value = TRUE)
test_tree_file <- grep("ASTRAL", grep("geneconv_pass", grep("Metazoan", all_trees, value = TRUE), value = TRUE), value = TRUE)
# Open trees and drop taxa
keep_tips <- c("Ocyropsis_sp_Florida_USA", "Bolinopsis_ashleyi", "Beroe_ovata", 
               "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Beroe_abyssicola", 
               "Beroe_sp_Antarctica","Beroe_sp_Queensland_Australia", "Beroe_forskalii",
               "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Lobata_sp_Punta_Arenas_Argentina",
               "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA")
notest_tree <- read.tree(notest_tree_file)
notest_tree <- keep.tip(notest_tree, keep_tips)
test_tree <- read.tree(test_tree_file)
test_tree <- keep.tip(test_tree, keep_tips)
# Change NaN edge lengths to 0.1
notest_tree$edge.length[is.nan(notest_tree$edge.length)] <- 0.1
test_tree$edge.length[is.nan(test_tree$edge.length)] <- 0.1
# Prepare labels
lab_df <- data.frame(taxa = keep_tips,
                     clean_taxa = gsub("_", " ", keep_tips),
                     color = c(rep("A", 1), rep("B", 2), (rep("C", length(keep_tips) - 3)) ) )
lab_df <- dplyr::mutate(lab_df, 
                        lab = glue('italic("{clean_taxa}")'))
# Plot each tree as a ggtree 
p1 <- ggtree(notest_tree)  %<+% lab_df + 
  geom_tiplab(aes(label = lab, color = color), size = 4, parse = T, show.legend = F) + 
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(subset = !isTip, label=label), nudge_x = -0.09, nudge_y = 0.20, color = "Gray55") +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) +
  theme(axis.text.x = element_text(size = 13, color = "White"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_color_manual(values = c(A = "Red", B = "Gray50", C = "Black"))

p2 <- ggtree(test_tree)  %<+% lab_df + 
  geom_tiplab(aes(label = lab, color = color), size = 4, parse = T, show.legend = F) + 
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(subset = !isTip, label=label), nudge_x = -0.08, nudge_y = 0.20, color = "Gray55") +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) +
  theme(axis.text.x = element_text(size = 13, color = "White"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_color_manual(values = c(A = "Red", B = "Gray50", C = "Black"))

quilt <- (p1 / p2) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
plot_id <- determine.outlier.plot.name(1, outlier_df)
quilt_name <- paste0(plot_dir, sprintf("%03d_", 1), plot_id, "OutlierBranch_plot")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 12, width = 10, units = "in")


