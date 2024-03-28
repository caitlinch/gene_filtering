### gene_filtering/code/6_Quartet_Concordance.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2024

## This script processes and plots ASTRAL trees with quartet concordance factors


##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir         <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir        <- for saving plots and analyses.
# qcf_dir         <- location of ASTRAL trees with qCFs
# output_files    <- location of csv files

### Caitlin's paths ###
# Folders and filepaths
maindir       <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
plot_dir      <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
qcf_dir       <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_qCF/"
output_dir    <- qcf_dir
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Source packages
library(ggplot2)
library(ggpubr) # To add significance values to plots
library(patchwork) # For assembling multi-panel plots

# Source functions
source(paste0(maindir, "code/func_comparison.R")) # includes ape, phangorn, phytools, dplyr

# Dataset information
dataset_names <- c("Plants" = "1KP", "Tomatoes" = "Pease2016", "Metazoa" = "Whelan2017", "Primates" = "Vanderpool2020")
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoa" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))



#### Step 3: Open quartet concordance factors ####
# Create output file names
tree_comp_df_file <- paste0(output_dir, "qCF_input_trees.csv")
qcf_df_file <- paste0(output_dir, "qCF_values.csv")

# Open input dataframe
if (file.exists(tree_comp_df_file)){
  tree_comp_df <- read.csv(tree_comp_df_file, stringsAsFactors = FALSE)
} else {
  # Extract all tree files
  all_files <- list.files(qcf_dir)
  tree_files <- grep(".tre", all_files, value = T)
  # Create dataframe for testing every pair of trees
  clean_trees       <- grep("pass", tree_files, value = T)
  recomb_trees      <- grep("fail", tree_files, value = T)
  alignment_trees   <- grep("NoTest", tree_files, value = T)
  tree_comp_df      <- data.frame(dataset = c(rep("Pease2016", 8),
                                              rep("Vanderpool2020",  8),
                                              rep("Whelan2017", 3),
                                              rep("1KP", 8)),
                                  clean_tree = c("Pease2016_allTests_pass_ASTRAL_qCF.tre", "Pease2016_allTests_pass_ASTRAL_qCF.tre",
                                                 "Pease2016_geneconv_pass_ASTRAL_qCF.tre", "Pease2016_geneconv_pass_ASTRAL_qCF.tre",
                                                 "Pease2016_maxchi_pass_ASTRAL_qCF.tre", "Pease2016_maxchi_pass_ASTRAL_qCF.tre",
                                                 "Pease2016_PHI_pass_ASTRAL_qCF.tre", "Pease2016_PHI_pass_ASTRAL_qCF.tre",
                                                 "Vanderpool2020_allTests_pass_ASTRAL_qCF.tre", "Vanderpool2020_allTests_pass_ASTRAL_qCF.tre",
                                                 "Vanderpool2020_geneconv_pass_ASTRAL_qCF.tre", "Vanderpool2020_geneconv_pass_ASTRAL_qCF.tre",
                                                 "Vanderpool2020_maxchi_pass_ASTRAL_qCF.tre", "Vanderpool2020_maxchi_pass_ASTRAL_qCF.tre",
                                                 "Vanderpool2020_PHI_pass_ASTRAL_qCF.tre", "Vanderpool2020_PHI_pass_ASTRAL_qCF.tre",
                                                 "Whelan2017_geneconv_pass_ASTRAL_qCF.tre", "Whelan2017_maxchi_pass_ASTRAL_qCF.tre",
                                                 "Whelan2017_PHI_pass_ASTRAL_qCF.tre",
                                                 "1KP_allTests_pass_ASTRAL_qCF.tre", "1KP_allTests_pass_ASTRAL_qCF.tre",
                                                 "1KP_geneconv_pass_ASTRAL_qCF.tre", "1KP_geneconv_pass_ASTRAL_qCF.tre",
                                                 "1KP_maxchi_pass_ASTRAL_qCF.tre", "1KP_maxchi_pass_ASTRAL_qCF.tre",
                                                 "1KP_PHI_pass_ASTRAL_qCF.tre" , "1KP_PHI_pass_ASTRAL_qCF.tre"),
                                  comparison_tree = c("Pease2016_allTests_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                      "Pease2016_geneconv_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                      "Pease2016_maxchi_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                      "Pease2016_PHI_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                      "Vanderpool2020_allTests_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                      "Vanderpool2020_geneconv_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                      "Vanderpool2020_maxchi_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                      "Vanderpool2020_PHI_fail_ASTRAL_qCF.tre" , "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                      "Whelan2017_NoTest_ASTRAL_qCF.tre", "Whelan2017_NoTest_ASTRAL_qCF.tre",
                                                      "Whelan2017_NoTest_ASTRAL_qCF.tre",
                                                      "1KP_allTests_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                      "1KP_geneconv_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                      "1KP_maxchi_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                      "1KP_PHI_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre"),
                                  tree_directory = qcf_dir,
                                  comparison_id = c("allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                    "allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                    "noTest", "noTest", "noTest",
                                                    "allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest"))
  tree_comp_df$clean_id <- unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){paste0(x[[2]], "_", x[[3]])}))
  tree_comp_df$recombination_test <-  unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){x[[2]]}))
  tree_comp_df$comparison_gene_status <- unlist(lapply(strsplit(tree_comp_df$comparison_tree, "_"), function(x){x[[3]]}))
  tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "ASTRAL")] <- "unfiltered"
  tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "fail")] <- "recombinant"
  # Save csv
  write.csv(tree_comp_df, file = tree_comp_df_file, row.names = FALSE)
}

# Open qCF dataframe
if (file.exists(qcf_df_file)){
  qcf_df <- read.csv(qcf_df_file, stringsAsFactors = FALSE)
} else {
  # Extract qCF values
  qcf_op <- lapply(1:nrow(tree_comp_df), compare.splits.wrapper, df = tree_comp_df)
  qcf_df <- as.data.frame(do.call(rbind, qcf_op))
  qcf_df <- qcf_df[ , c("dataset", "clean_tree", "comparison_tree", "clean_id", "comparison_id", "recombination_test", "comparison_gene_status",
                        "tree", "split_type", "weights", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")]
  # Save csv
  write.csv(qcf_df, file = qcf_df_file, row.names = FALSE)
}



#### Step 4: Plot quartet concordance factors ####
# Remove any splits that are not either "Concordant" or "Conflicting"
qcf_df <- qcf_df[which(qcf_df$split_type %in% c("Concordant", "Conflicting")) , ]
qcf_df$split_type[which(qcf_df$split_type == "Concordant")] <- "Congruent"
# Create a nicely formatted dataset column
qcf_df$dataset_formatted <- factor(qcf_df$dataset,
                                   levels = c("Pease2016", "Vanderpool2020", "1KP", "Whelan2017"),
                                   labels = c("Tomatoes", "Primates", "Plants", "Metazoa"),
                                   ordered = T)
# Add a new columns for faceting/grouping
qcf_df$gene_tree_formatted <- factor(qcf_df$tree,
                                     levels = c("geneconv_pass", "geneconv_fail", "maxchi_pass", "maxchi_fail",
                                                "PHI_pass", "PHI_fail", "allTests_pass", "allTests_fail", 
                                                "noTest"),
                                     labels = c("Clean", "Recombinant", "Clean", "Recombinant",
                                                "Clean", "Recombinant", "Clean", "Recombinant",
                                                "Unfiltered"),
                                     ordered = TRUE)
qcf_df$recombination_test_formatted <- factor(qcf_df$recombination_test,
                                              levels = c("geneconv", "maxchi", "PHI", "allTests"),
                                              labels = c("GENCONV", "MaxChi", "PHI", "All tests"),
                                              ordered = T)

## Separate datasets into separate plots
# Separate into dataframes for the different datasets
shallow_qcf <- qcf_df[which(qcf_df$dataset %in% c("Pease2016", "Vanderpool2020")), ]
plant_qcf   <- qcf_df[which(qcf_df$dataset == "1KP"), ]
metazoa_qcf <- qcf_df[which(qcf_df$dataset == "Whelan2017"), ]
# Remove any rows from the Plants dataset for genes that failed the test, and for GENECONV and allTests
plant_qcf <- plant_qcf[which(plant_qcf$clean_tree %in% c("1KP_maxchi_pass_ASTRAL_qCF.tre", "1KP_PHI_pass_ASTRAL_qCF.tre")), ]
plant_qcf <- plant_qcf[which(plant_qcf$comparison_tree %in% c("1KP_NoTest_ASTRAL_qCF.tre")), ]

## Plot the shallow datasets (Primates and Tomatoes)
shallow_plot <- ggplot(shallow_qcf, aes(x = gene_tree_formatted, y = q1, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "qCF", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size=8)))

## Plot the Metazoan dataset
metazoa_plot <- ggplot(metazoa_qcf, aes(x = gene_tree_formatted, y = q1, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "qCF", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size=8)))

## Plot the Plants dataset
plants_plot <- ggplot(plant_qcf, aes(x = gene_tree_formatted, y = q1, fill = split_type)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Gene filtering") +
  scale_y_continuous(name = "qCF", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "c.") +
  scale_fill_manual(name = "Branch type", values = c("#a6cee3", "#1f78b4"), labels = c("Congruent", "Conflicting")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size=8)))

## Save plots as a patchwork
quilt <- shallow_plot + metazoa_plot + plants_plot + plot_layout(ncol = 1, heights = c(2,1,1))
quilt_pdf <- paste0(plot_dir, "qcf_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 12, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "qcf_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 12, width = 10, units = "in")

## Statistics plots for deep datasets that have sufficient conflicting branches to calculate p-values
# Metazoa
metazoa_stat_plot <- ggplot(metazoa_qcf, aes(x = split_type, y = q1, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "qCF", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "a.") +
  scale_fill_manual(name = "Gene filtering", values = c("#0571b0", "#f7f7f7"), labels = c("Clean", "Unfiltered")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size=8)))
metazoa_p_values <- metazoa_stat_plot  + stat_compare_means(method = "t.test", aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 4, color = "darkgrey")
# Plants
plants_stat_plot <- ggplot(plant_qcf, aes(x = split_type, y = q1, fill = gene_tree_formatted)) +
  geom_boxplot() +
  facet_grid(dataset_formatted~recombination_test_formatted) +
  scale_x_discrete(name = "Branch type") +
  scale_y_continuous(name = "qCF", breaks = seq(0,1,0.2),  labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1), limits = c(0,1)) +
  labs(title = "b.") +
  scale_fill_manual(name = "Gene filtering", values = c("#0571b0", "#f7f7f7"), labels = c("Clean", "Unfiltered")) +
  theme_bw() +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size=8)))
plants_p_values <- plants_stat_plot  + stat_compare_means(method = "t.test", aes(label = paste0("p=", ..p.format..)), label.y = c(0.2), size = 5, color = "darkgrey")
## Save plots as a patchwork
quilt <- metazoa_p_values + plants_p_values + plot_layout(ncol = 1)
quilt_pdf <- paste0(plot_dir, "qcf_p_value_quilt.pdf")
ggsave(filename = quilt_pdf, plot = quilt, height = 8, width = 10, units = "in")
quilt_png <- paste0(plot_dir, "qcf_p_value_quilt.png")
ggsave(filename = quilt_png, plot = quilt, height = 8, width = 10, units = "in")
