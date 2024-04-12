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
# annotations_csv_file      <- location of misc/annotations.csv file from the Leebens-Mack (2019) "Data from 1000 Plants Transcriptomes" data repository - used to assign taxa names and clades

# roots_by_group            <- set which taxa is outgroup for each dataset

### Caitlin's paths ###
location = "local"
if (location == "local"){
  maindir             <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
  tree_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/Species_trees_for_publication/"
  test_data_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/"
  output_dir          <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"
  plot_dir            <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
  annotation_csv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/misc/annotations.csv"
}


# Add roots for each dataset
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))
### End of Caitlin's paths ###



##### Step 2: Open packages and set directories #####
# Open packages
library(ape)
library(distory)
library(ggplot2)
library(ggtree)
library(ggpubr)
library(patchwork)
library(colorBlindness)



##### Step 3: Source function files and prepare variables for analysis #####
# Open annotations for Plants dataset
annotation_df <- read.csv(annotation_csv_file, stringsAsFactors = F)

# Create color palettes
primate_colour_palette <-  c("Non-primate"="#332288", "Strepsirrhini"="#117733", "Tarsiiformes"="#44AA99",
                             "Cebidae"="#88CCEE", "Hylobatidae"="#DDCC77", "Hominidae"="#CC6677",
                             "Colobinae"="#AA4499", "Cercopithecinae"="#882255")
tomato_colour_palette <- c("Esculentum" = "firebrick3", "Arcanum" = "goldenrod3", 
                           "Peruvianum" = "darkgreen", "Hirsutum" = "navy", 
                           "Outgroup" = "black")
metazoan_colour_palette <- c("Bilateria" = "#CC79A7", "Cnidaria" = "#009E73", "Ctenophora" = "#56B4E9",
                             "Porifera" = "#E69F00", "Outgroup" = "#999999", "Placozoa" = "#000000")
plants_color_palette <- c(SteppedSequential5Steps[c(1,3,5,6,8,10,11,13,15,16,18,20,21,23,25)], "black", "grey40", "grey70")
plant_classifications <-  c("Chromista", "Rhodophyta", "Glaucophyta", "Chlorophyta", "Streptophyte algae", 
                            "Hornworts", "Liverworts", "Mosses", "Lycophytes", "Monilophytes", "Gymnos",
                            "ANAGrade", "Monocots", "Magnoliids", "Chloranthales", "Eudicots", "Ceratophyllales")
names(plants_color_palette) <- plant_classifications

# Open functions
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
                                             labels = c("GENECONV", "MaxChi", "PHI", "All tests"),
                                             ordered = T)
# Remove any test, fail rows
pp_df <- pp_df[grep("fail", pp_df$comparison_tree, invert = T), ]
# Convert confidence to numeric
pp_df$confidence <- as.numeric(pp_df$confidence)
# Separate into dataframes for the different datasets
shallow_pp <- pp_df[which(pp_df$dataset %in% c("Tomatoes", "Primates")), ]
tomatoes_pp <- pp_df[which(pp_df$dataset == "Tomatoes"), ]
primates_pp <- pp_df[which(pp_df$dataset == "Primates"), ]
plant_pp   <- pp_df[which(pp_df$dataset == "Plants"), ]
metazoan_pp <- pp_df[which(pp_df$dataset == "Metazoan"), ]
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot shallow datasets
shallow_pp_plot <- ggplot(shallow_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot tomatoes
tomatoes_pp_plot <- ggplot(tomatoes_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_pp_plot <- ggplot(primates_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_pp_plot <- ggplot(metazoan_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_pp_plot <- ggplot(plant_pp, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "d.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_pp_plot + primates_pp_plot + met_pp_plot + plants_pp_plot + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "PosteriorProbabilities_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "PosteriorProbabilities_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")

## Add p-values
# Plot tomatoes
tomatoes_pp_plot <- ggplot(tomatoes_pp, aes(x = split_type, y = confidence, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
tomatoes_pp_pvals <- tomatoes_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 4, color = "darkgrey")
# Plot primates
primates_pp_plot <- ggplot(primates_pp, aes(x = split_type, y = confidence, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
primates_pp_pvals <- primates_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 4, color = "darkgrey")
# Plot metazoans
met_pp_plot <- ggplot(metazoan_pp, aes(x = split_type, y = confidence, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
met_pp_pvals <- met_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 4, color = "darkgrey")
# Plot plants
plants_pp_plot <- ggplot(plant_pp, aes(x = split_type, y = confidence, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Lpp", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
plants_pp_pvals <- plants_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 4, color = "darkgrey")
# Save
quilt <- tomatoes_pp_pvals + primates_pp_pvals + met_pp_pvals + plants_pp_pvals + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "PosteriorProbabilities_quilt_stats.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "PosteriorProbabilities_quilt_stats.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")



##### Step 6: Plot branch lengths for ASTRAL trees #####
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes
tomatoes_bla_plot <- ggplot(tomatoes_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_bla_plot <- ggplot(primates_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,4.5)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_bla_plot <- ggplot(metazoan_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5, 1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,3.5)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_bla_plot <- ggplot(plant_pp, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,5.5)) +
  labs(title = "d.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_bla_plot + primates_bla_plot + met_bla_plot + plants_bla_plot + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "BranchLengths_ASTRAL_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 13)
quilt_png <- paste0(plot_dir, "BranchLengths_ASTRAL_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 13)

## Add p-values
# Plot tomatoes
tomatoes_pp2 <- tomatoes_pp[which(tomatoes_pp$recombination_test == "PHI" & tomatoes_pp$split_type == "Conflicting"), ]
tomatoes_pp_plot <- ggplot(tomatoes_pp, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
tomatoes_pp_pvals <- tomatoes_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(6.2), size = 4, color = "darkgrey")
# Plot primates
primates_pp_plot <- ggplot(primates_pp, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
primates_pp_pvals <- primates_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(6.2), size = 4, color = "darkgrey")
# Plot metazoans
met_pp_plot <- ggplot(metazoan_pp, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
met_pp_pvals <- met_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(6.2), size = 4, color = "darkgrey")
# Plot plants
plants_pp_plot <- ggplot(plant_pp, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,6.5,1),  labels = seq(0,6.5,1), minor_breaks = seq(0,6.5,0.5), limits = c(0,6.5)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
plants_pp_pvals <- plants_pp_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(6.2), size = 4, color = "darkgrey")
# Save
quilt <- tomatoes_pp_pvals + primates_pp_pvals + met_pp_pvals + plants_pp_pvals + plot_layout(ncol = 1) 
quilt_pdf <- paste0(plot_dir, "ASTRAL_branchLength_quilt_stats.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "ASTRAL_branchLength_quilt_stats.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")



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
                                             labels = c("GENECONV", "MaxChi", "PHI", "All tests"),
                                             ordered = T)
# Remove any test, fail rows
bs_df <- bs_df[grep("fail", bs_df$comparison_tree, invert = T), ]
# Convert confidence to numeric
bs_df$confidence <- as.numeric(bs_df$confidence)
# Separate into dataframes for the different datasets
shallow_bs  <- bs_df[which(bs_df$dataset %in% c("Tomatoes", "Primates")), ]
tomatoes_bs <- bs_df[which(bs_df$dataset == "Tomatoes"), ]
primates_bs <- bs_df[which(bs_df$dataset == "Primates"), ]
metazoan_bs <- bs_df[which(bs_df$dataset == "Metazoan"), ]
plant_bs <-  bs_df[which(bs_df$dataset == "Plants"), ]
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes 
tomatoes_bs_plot <- ggplot(tomatoes_bs, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot()  +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset")+
  scale_y_continuous(name = "UFB", breaks = seq(0,120,20),  labels = seq(0,120,20), minor_breaks = seq(0,110,10), limits = c(0,110)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_bs_plot <- ggplot(primates_bs, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot()  +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset")+
  scale_y_continuous(name = "UFB", breaks = seq(0,120,20),  labels = seq(0,120,20), minor_breaks = seq(0,110,10), limits = c(0,110)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_bs_plot <- ggplot(metazoan_bs, aes(x = gene_tree_formatted, y = confidence, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "UFB", breaks = seq(0,120,20),  labels = seq(0,120,20), minor_breaks = seq(0,110,10), limits = c(0,110)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_bs_plot + primates_bs_plot + met_bs_plot + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "UltrafastBootstrap_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 10, width = 8)
quilt_png <- paste0(plot_dir, "UltrafastBootstrap_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 10, width = 8)



##### Step 8: Plot branch lengths for CONCAT trees #####
# Save theming as object 
theming <- theme_bw() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
# Plot tomatoes
tomatoes_blc_plot <- ggplot(tomatoes_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,0.02,0.005),  labels = seq(0,0.02,0.005), minor_breaks = seq(0,0.02,0.0025), limits = c(0,0.02)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot primates
primates_blc_plot <- ggplot(primates_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,0.04,0.01),  labels = seq(0,0.04,0.01), minor_breaks = seq(0,0.04,0.005), limits = c(0,0.04)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot metazoans
met_blc_plot <- ggplot(metazoan_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,0.5,0.1),  labels = seq(0,0.5,0.1), minor_breaks = seq(0,0.5,0.10), limits = c(0,0.5)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Plot plants
plants_blc_plot <- ggplot(plant_bs, aes(x = gene_tree_formatted, y = weights, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset") +
  scale_y_continuous(name = "Branch length", breaks = seq(0,1.20,0.20),  labels = seq(0,1.20,0.20), minor_breaks = seq(0,1.20,0.10), limits = c(0,1.20)) +
  labs(title = "d.") +
  scale_fill_manual(name = "Branch type", values = c("Congruent" = "#a6cee3", "Conflicting" = "#1f78b4")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
# Save
quilt <- tomatoes_blc_plot + primates_blc_plot + met_blc_plot + plants_blc_plot + plot_layout(ncol = 1, heights = c(1,1,1,1))
quilt_pdf <- paste0(plot_dir, "BranchLengths_CONCAT_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12)
quilt_png <- paste0(plot_dir, "BranchLengths_CONCAT_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 12)

## Add p-values
# Separate into dataframes for the different datasets
tomatoes_bs <- bs_df[which(bs_df$dataset == "Tomatoes"), ]
primates_bs <- bs_df[which(bs_df$dataset == "Primates"), ]
metazoan_bs <- bs_df[which(bs_df$dataset == "Metazoan"), ]
plant_bs <- bs_df[which(bs_df$dataset == "Plants"), ]
# Plot tomatoes
tomatoes_bs_plot <- ggplot(tomatoes_bs, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot()  +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Subset")+
  scale_y_continuous(name = "Branch length") +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
tomatoes_bs_pvals <- tomatoes_bs_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.02), size = 4, color = "darkgrey")
# Plot primates
primates_bs_plot <- ggplot(primates_bs, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length") +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
primates_bs_pvals <- primates_bs_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.035), size = 4, color = "darkgrey")
# Plot metazoans
met_bs_plot <- ggplot(metazoan_bs, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length") +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
met_bs_pvals <- met_bs_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(0.5), size = 4, color = "darkgrey")
# Plot metazoans
plant_bs_plot <- ggplot(plant_bs, aes(x = split_type, y = weights, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "Branch length") +
  labs(title = "a.") +
  scale_fill_manual(name = "Subset", values = c("Clean" = "#0571b0", "Unfiltered" = "#f7f7f7")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theming
plant_bs_pvals <- plant_bs_plot  + stat_compare_means(method = "t.test", paired = T, aes(label = paste0("p=", ..p.format..)), label.y = c(1.2), size = 4, color = "darkgrey")
# Save
quilt <- tomatoes_bs_pvals + primates_bs_pvals + met_bs_pvals + + plot_layout(ncol = 1) 
quilt_pdf <- paste0(plot_dir, "CONCAT_UFB_quilt_stats.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "CONCAT_UFB_quilt_stats.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")



##### Step 6: Identify outlier branches (ones with high support) #####
branch_support_df$confidence <- as.numeric(branch_support_df$confidence)
outlier_df <- branch_support_df[c(which(branch_support_df$analysis_method == "ASTRAL" & branch_support_df$confidence > 0.9 & branch_support_df$split_type == "Conflicting"), 
                                  which(branch_support_df$analysis_method == "IQTREE" & branch_support_df$confidence > 90 & branch_support_df$split_type == "Conflicting")), ]
outlier_df_path <- paste0(output_dir, "Outlier_Branches_raw.csv")
write.csv(outlier_df, file = outlier_df_path)



##### Step 7: Highlight outlier branches #####
old.par <- par(mar = c(0, 0, 0, 0))
par(mfrow = c(1, 1))
# Remove any "fail" rows
outlier_df <- outlier_df[which(outlier_df$comparison_id == "noTest"), ]
# Add new column of which tree to open
outlier_df$open_tree <- ""
outlier_df$open_tree[which(outlier_df$tree == "noTest")] <- outlier_df$comparison_tree[which(outlier_df$tree == "noTest")]
outlier_df$open_tree[which(outlier_df$tree != "noTest")] <- outlier_df$clean_tree[which(outlier_df$tree != "noTest")]
# Extract unique trees
input_trees <-  unique(outlier_df$open_tree)
# Need i = 8
for (i in 1:length(input_trees)){
  # Open tree file
  i_tree_file <- input_trees[i]
  i_tree <- read.tree(paste0(tree_data_dir, i_tree_file))
  # Identify dataset and plot with highlighted outlier branches
  if (i_rows$dataset[1] == "Tomatoes"){
    # Extract rows for this tree
    i_rows <- outlier_df[which(outlier_df$open_tree == i_tree_file), ]
    # Identify branches in tree
    i_outlier_branches <- which(round(i_tree$edge.length, digits = 5) %in% round(i_rows$weights, digits = 5))
    # Update node labels
    branch_labels = rep(NA, length(i_tree$edge.length))
    branch_labels[i_outlier_branches] <- "Conflicting"
    # Rename tip labels to have scientific names (not just numbers)
    i_tree$tip.label <- rename.tomato.tips(i_tree$tip.label)
    # Create plot
    plot(i_tree)
    title(main = "ASTRAL - Tomatoes - Pass PHI", sub = "Outlier branches")
    edgelabels(branch_labels)
    # Create plot name
    p_name <- paste0(plot_dir, i_tree_file, "_OutlierBranches")
    print(p_name)
    # Save plot
    ggsave(filename = paste0(p_name, ".pdf"), plot = plot, device = "pdf")
  } else if (i_rows$dataset[1] == "Metazoan"){
    # Reroot
    i_tree <- root(i_tree, outgroup = roots_by_group[["Metazoan"]], resolve.root = TRUE)
    # Extract rows for this tree
    i_rows <- outlier_df[which(outlier_df$open_tree == i_tree_file), ]
    # Identify branches in tree
    i_outlier_branches <- unlist(lapply(round(i_rows$weights, digits = 5), function(x){which(round(i_tree$edge.length, digits = 5) == x)}))
    # Update node labels
    branch_labels = rep(NA, length(i_tree$edge.length))
    branch_labels[i_outlier_branches] <- "Conflicting"
    # Create plot
    plot(i_tree, cex = 0.5)
    title(main = "CONCAT - Metazoan - Pass PHI", sub = "Outlier branches")
    edgelabels(branch_labels)
    # Create plot name
    p_name <- paste0(plot_dir, i_tree_file, "_OutlierBranches")
    print(p_name)
    # Save plot
    ggsave(filename = paste0(p_name, ".pdf"), plot = plot, device = "pdf")
  }
}
# For Plants: i = 8
i = 8
i_tree_file <- input_trees[i]
i_tree <- read.tree(paste0(tree_data_dir, i_tree_file))
# Extract rows for this tree
i_rows <- outlier_df[which(outlier_df$open_tree == i_tree_file), ]
# Identify branches in tree
i_outlier_branches <- which(round(i_tree$edge.length, digits = 5) %in% round(i_rows$weights, digits = 5))
# Extract clade
clade <- extract.clade(i_tree, 1306)
clade_outlier_branch <- which(round(clade$edge.length, digits = 5) %in% round(i_rows$weights, digits = 5))
# Update node labels
branch_labels = rep(NA, length(clade$edge.length))
branch_labels[clade_outlier_branch] <- "NoTest"
# Create labels
i_labs <- color.plants.by.clades(tree = i_tree, color_palette = plants_color_palette, clade_df = annotation_df)
# Trim to only labels in the clade
i_labs <- i_labs[which(i_labs$Code %in% clade$tip.label), ]
clade$tip.label <- unlist(lapply(clade$tip.label, function(x){i_labs[which(i_labs$Code == x), ]$Species}))
# Create plot
plot(clade, cex = 1)
title(main = "ASTRAL - Plants - Pass PHI\nClade in: CoreEudicots/Rosids", sub = "Outlier branches")
edgelabels(branch_labels)
# Create plot name
p_name <- paste0(plot_dir, i_tree_file, "_OutlierBranches")
print(p_name)
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = plot, device = "pdf")




