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
branch_support_df$confidence <- as.numeric(branch_support_df$confidence)
outlier_df <- branch_support_df[c(which(branch_support_df$analysis_method == "ASTRAL" & branch_support_df$confidence > 0.9 & branch_support_df$split_type == "Conflicting"), 
                                  which(branch_support_df$analysis_method == "IQTREE" & branch_support_df$confidence > 90 & branch_support_df$split_type == "Conflicting")), ]
outlier_df_path <- paste0(output_dir, "Outlier_Branches.csv")
write.csv(outlier_df, file = outlier_df_path)



