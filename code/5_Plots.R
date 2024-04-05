### gene_filtering/code/5_Plots.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2023

## This script creates a variety of plots and figures from the subset trees and species trees for the 4 empirical datasets

## ggtree notes:
# To add node labels (for defining clade labels): geom_text(aes(label=node), hjust=-.3)
# To add clade labels (for defining clades with a vertical bar): geom_cladelab(node=34, label="Clade")



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                 <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir                <- for saving plots and analyses.
# annotations_csv_file    <- location of misc/annotations.csv file from the Leebens-Mack (2019) "Data from 1000 Plants Transcriptomes" data repository
# roots_by_groups         <- set which taxa is outgroup for each dataset: Primates, Tomatoes, Metazoans, and Plants

### Caitlin's paths ###
# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
annotation_csv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/misc/annotations.csv"

# Dataset information
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ape) # functions: read.tree, Ntip, root
library(Cairo) # for plotting png
library(phangorn)
library(ggplot2) # for nice plots
library(ggtree) # for plotting phylogenetic trees and densitress (ggdensitree)
library(ggtext) # for nice tree plots
library(ggpubr) # for textGrob function
library(patchwork) # for collating plots
library(TreeTools) # for CollapseNode function
library(colorBlindness) # For Plants clade labels

# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "species_trees/")

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




#### Step 3: Plotting Primates dataset ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
primate_tree_files <- grep("Primates", tree_files, value = TRUE)
primate_astral_trees <- grep("ASTRAL", primate_tree_files, value = TRUE)
primate_concat_trees <- grep("CONCAT", primate_tree_files, value = TRUE)

## Primate Plot 1: ASTRAL No Test
# Assemble file path and open tree
pt1_file <- grep("NoTest", primate_astral_trees, value = TRUE)
pt1 <- read.tree(pt1_file)
# Root tree (as in Vanderpool 2020 paper)
pt1 <- root(pt1, outgroup = roots_by_group[["Primates"]])
# Change edge.length to 0.5
pt1$edge.length[which(is.nan(pt1$edge.length))] <- 0.5
# Remove underscores from tips
pt1$tip.label <- gsub("_", " ", pt1$tip.label)
# Color code clades
pt1_colors <- primate_colour_palette
pt1_labs <-color.primates.by.clades(pt1, color_palette = pt1_colors)
# Create plot
p <- ggtree(pt1) %<+% pt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0.002, geom = "text", size = 5) + 
  geom_rootedge(rootedge = 0.3, size = 0.5) +
  scale_y_reverse() +  
  scale_x_continuous(breaks = seq(0,15,3)) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) + 
  scale_color_manual(values = pt1_colors) +
  theme(axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 15), legend.text = element_text (size = 12), legend.position = c(0.1,0.2)) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.")))
quilt <-  p + 
  plot_annotation(title = "Primate dataset - ASTRAL tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
# Create plot name
p_name <- paste0(plot_dir, "Primates_ASTRAL_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf")

## Primate Plot 2: ASTRAL variable clade (Cebidae)
# Three possible topologies:
#   1. No Test; P,PHI; P,MaxChi; P,GENECONV; and P,All tests 
#   2. F, MaxChi 
#   3. F,PHI; F,GENECONV; and F,All tests

# Assemble one file path to one tree for each possible topology
pt1_file <- grep("NoTest", primate_astral_trees, value = TRUE)
pt2_file <- grep("maxchi_fail", primate_astral_trees, value = TRUE)
pt3_file <- grep("PHI_fail", primate_astral_trees, value = TRUE)
# Open one tree from each of the possible topologies
pt1 <- read.tree(file = pt1_file)
pt2 <- read.tree(file = pt2_file)
pt3 <- read.tree(file = pt3_file)
# Root trees (as in Vanderpool 2020 paper)
pt1 <- root(pt1, outgroup = roots_by_group[["Primates"]])
pt2 <- root(pt2, outgroup = roots_by_group[["Primates"]])
pt3 <- root(pt3, outgroup = roots_by_group[["Primates"]])
# Drop all tips except the four tips inside the variable clade
pt1 <- keep.tip(pt1, c("Callithrix_jacchus", "Aotus_nancymaae", "Saimiri_boliviensis", "Cebus_capucinus_imitator"))
pt2 <- keep.tip(pt2, c("Callithrix_jacchus", "Aotus_nancymaae", "Saimiri_boliviensis", "Cebus_capucinus_imitator"))
pt3 <- keep.tip(pt3, c("Callithrix_jacchus", "Aotus_nancymaae", "Saimiri_boliviensis", "Cebus_capucinus_imitator"))
# Rename tips to remove underscores
pt1$tip.label <- gsub("_", " ", pt1$tip.label)
pt2$tip.label <- gsub("_", " ", pt2$tip.label)
pt3$tip.label <- gsub("_", " ", pt3$tip.label)
# Change edge.length to 0.1
pt1$edge.length[which(is.nan(pt1$edge.length))] <- 0.1
pt2$edge.length[which(is.nan(pt2$edge.length))] <- 0.1
pt3$edge.length[which(is.nan(pt3$edge.length))] <- 0.1
# Create plot for each topology
p1_0 <- ggtree(pt1, linewidth = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, linewidth = 1) + xlim(-0.1, 1.2)
p2 <- ggtree(pt2, linewidth = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, linewidth = 1) + xlim(-0.1, 1.2)
p3 <- ggtree(pt3, linewidth = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, linewidth = 1) + xlim(-0.1, 1.2)
# Flip pt1 branch so all 3 plots have same order
p1 <- flip(p1_0, 2, 3)
# Create plot name 
p_name <- paste0(plot_dir, "Primates_ASTRAL_Cebidae_variable_clade_a.NoTest_b.F.maxchi_c.F.PHI_plot")
# Assemble plot in patchwork
p <- p1 + p2 + p3 + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 60))
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = p, device = "pdf", height = 10, width = 32, units = "in")

## Primate Plot 3: CONCAT No Test
# Assemble file path and open tree
pt1_file <- grep("NoTest", primate_concat_trees, value = TRUE)
pt1 <- read.tree(pt1_file)
# Root tree (as in Vanderpool 2020 paper)
pt1 <- root(pt1, outgroup = roots_by_group[["Primates"]])
# Remove underscores from tips
pt1$tip.label <- gsub("_", " ", pt1$tip.label)
# Color code clades
pt1_colors <- primate_colour_palette
pt1_labs <-color.primates.by.clades(pt1, color_palette = pt1_colors)
# Create plot
p <- ggtree(pt1) %<+% pt1_labs +
  geom_rootedge(rootedge = 0.005, size = 0.5) +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0.002, geom = "text", size = 5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) + 
  scale_color_manual(values = pt1_colors) +
  theme(axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 15), legend.text = element_text (size = 12), legend.position = c(0.1,0.2)) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.")))
quilt <-  p + 
  plot_annotation(title = "Primate dataset - Concatenated tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
# Create plot name
p_name <- paste0(plot_dir, "Primates_CONCAT_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf")

## Primate Plot 4: CONCAT variable clade (Papionini)
# Three possible topologies:
#   1. No Test; P,PHI; P,MaxChi; P,GENECONV; P,All tests; F,MaxChi; F,GENECONV; and F,All tests 
#   2. F,PHI

# Assemble one file path to one tree for each possible topology
pt1_file <- grep("NoTest", primate_concat_trees, value = TRUE)
pt2_file <- grep("PHI_fail", primate_concat_trees, value = TRUE)
# Open one tree from each of the possible topologies
pt1 <- read.tree(file = pt1_file)
pt2 <- read.tree(file = pt2_file)
# Root trees (as in Vanderpool 2020 paper)
pt1 <- root(pt1, outgroup = roots_by_group[["Primates"]])
pt2 <- root(pt2, outgroup = roots_by_group[["Primates"]])
# Drop all tips except the four tips inside the variable clade
pt1 <- keep.tip(pt1, c("Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"))
pt2 <- keep.tip(pt2, c("Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"))
# Rename tips to remove underscores
pt1$tip.label <- gsub("_", " ", pt1$tip.label)
pt2$tip.label <- gsub("_", " ", pt2$tip.label)
# Create plot for each topology
p1 <- ggtree(pt1, size = 1) + geom_tiplab(offset = 0, geom = "text", size = 8) + 
  geom_rootedge(rootedge = 0.0005, size = 1) + xlim(-0.0005, 0.01)
p2 <- ggtree(pt2, size = 1) + geom_tiplab(offset = 0, geom = "text", size = 8) + 
  geom_rootedge(rootedge = 0.0005, size = 1) + xlim(-0.0005, 0.025)
# Create plot name 
p_name <- paste0(plot_dir, "Primates_CONCAT_Papionini_variable_clade_a.NoTest_b.F.PHI_plot")
# Assemble plot in patchwork
p <- p1 + p2 +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 45))
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = p, device = "pdf", height = 10, width = 15, units = "in")



#### Step 4: Plotting differences in Primate trees for supplementary data ####
## Get the paths for the tree files
cebidae_tree_file <- paste0(maindir, "primate_tree_topologies/Cebidae_three_possible_topologies.txt")
comparison_tree_file <- paste0(maindir, "primate_tree_topologies/ComparisonTrees_three_possible_topologies.txt")
comparison_clade_file <- paste0(maindir, "primate_tree_topologies/ComparisonTrees_Clades_phylo.txt")

# Palettes
palette1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#999999")
palette2 <- c("#D81B60", "#1E88E5", "#E0A800", "#004D40", "#b1b1b1")

## Plot Cebidae trees
cebidae_trees <- read.tree(cebidae_tree_file)
# Root tree (as in Vanderpool 2020 paper)
cebidae_trees <- root(cebidae_trees, outgroup = roots_by_group[["Primates"]])
# Create color scheme
cebidae_df <- color.code.comparison.clades(cebidae_trees, variable = "Cebidae")

# Plot each tree
p1 <- ggtree(cebidae_trees[[1]])  %<+% cebidae_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4], clade_e = palette2[5]))

p2 <- ggtree(cebidae_trees[[2]])  %<+% cebidae_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4], clade_e = palette2[5]))

p3 <- ggtree(cebidae_trees[[3]])  %<+% cebidae_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4], clade_e = palette2[5]))

# Assemble patchwork
quilt <- p1 / p2 / p3 +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
# Create plot name
quilt_name <- paste0(plot_dir, "Primates_comparison_trees_Cebidae_clade")
# Save plot
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 15, width = 7, units = "in")

## Plot comparison trees
comparison_trees <- read.tree(comparison_tree_file)
# Root tree (as in Vanderpool 2020 paper)
comparison_trees <- root(comparison_trees, outgroup = roots_by_group[["Primates"]])
# Create color scheme
comparison_df <- color.code.comparison.clades(comparison_trees, variable = "Comparison")

# Plot each tree
p1 <- ggtree(comparison_trees[[1]])  %<+% comparison_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4]))

p2 <- ggtree(comparison_trees[[2]])  %<+% comparison_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4]))

p3 <- ggtree(comparison_trees[[3]])  %<+% comparison_df + 
  geom_tiplab(aes(label = lab, color = clade), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) +
  theme(axis.text.x = element_text(size = 13)) +
  scale_color_manual(values = c(clade_a = palette2[3], clade_b = palette2[2], clade_c = palette2[1], clade_d = palette2[4]))

# Assemble patchwork
quilt <- p1 / p2 / p3 +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
# Create plot name
quilt_name <- paste0(plot_dir, "Primates_comparison_trees_deep_branch")
# Save plot
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 15, width = 7, units = "in")

## Plot comparison tree clades
# Open clades as phylo objects
comparison_clades <- read.tree(comparison_clade_file)
# Make labels for each clade
c1_labs <- comparison.clade.tip.labels(comparison_clades[[1]])
c2_labs <- comparison.clade.tip.labels(comparison_clades[[2]])
c3_labs <- comparison.clade.tip.labels(comparison_clades[[3]])
c4_labs <- comparison.clade.tip.labels(comparison_clades[[4]])

# Plot each  clade
p1 <- ggpubr::text_grob(c1_labs$clean_taxa_names[[1]], face = "italic", size = 12)

p2 <- ggtree(comparison_clades[2])  %<+% c2_labs + 
  geom_tiplab(aes(label = lab), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  xlim(0, 7.5)

p3 <- ggtree(comparison_clades[3])  %<+% c3_labs + 
  geom_tiplab(aes(label = lab), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  xlim(0, 9)

p4 <- ggtree(comparison_clades[4])  %<+% c4_labs + 
  geom_tiplab(aes(label = lab), offset = 0, geom = "text", size = 4, parse = TRUE, show.legend = FALSE) +
  xlim(0, 30)

# Assemble patchwork 
quilt = (p2 + p1) / (p3 | p4) + plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
# Create plot name
quilt_name <- paste0(plot_dir, "Primates_comparison_trees_deep_branch_individal_clades")
# Save plot
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf")



#### Step 5: Plotting Tomatoes dataset ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
tomatoes_tree_files <- grep("Tomatoes", tree_files, value = TRUE)
tomatoes_astral_trees <- grep("ASTRAL", tomatoes_tree_files, value = TRUE)
tomatoes_concat_trees <- grep("CONCAT", tomatoes_tree_files, value = TRUE)

## ASTRAL Peruvianum topologies
# Need:
#     1. NoTest (identical to F,GENECONV and F,All tests)
#     2. P,PHI
#     3. F,PHI
#     4. P,MaxChi (identical to P,GENECONV)
#     5. F,MaxChi
#     6. P,All tests

# Find file for each tree
tt1_file <- grep("NoTest", tomatoes_astral_trees, value = TRUE)
tt2_file <- grep("PHI_pass", tomatoes_astral_trees, value = TRUE)
tt3_file <- grep("PHI_fail", tomatoes_astral_trees, value = TRUE)
tt4_file <- grep("maxchi_pass", tomatoes_astral_trees, value = TRUE)
tt5_file <- grep("maxchi_fail", tomatoes_astral_trees, value = TRUE)
tt6_file <- grep("allTests_pass", tomatoes_astral_trees, value = TRUE)
# Open each tree
tt1 <- read.tree(tt1_file)
tt1_small <- read.tree(tt1_file)
tt2 <- read.tree(tt2_file)
tt3 <- read.tree(tt3_file)
tt4 <- read.tree(tt4_file)
tt5 <- read.tree(tt5_file)
tt6 <- read.tree(tt6_file)
# Reroot tree 1 at proper outgroup
tt1 <- root(tt1, outgroup = roots_by_group[["Tomatoes"]])
# Reformat tomato clades for trees 1-6
tt1_small <- reformat.congruent.tomato.clades(tt1_small)
tt2 <- reformat.congruent.tomato.clades(tt2)
tt3 <- reformat.congruent.tomato.clades(tt3)
tt4 <- reformat.congruent.tomato.clades(tt4)
tt5 <- reformat.congruent.tomato.clades(tt5)
tt6 <- reformat.congruent.tomato.clades(tt6)
# Rename tip labels to have scientific names (not just numbers)
tt1$tip.label <- rename.tomato.tips(tt1$tip.label)
tt1_small$tip.label <- rename.tomato.tips(tt1_small$tip.label)
tt2$tip.label <- rename.tomato.tips(tt2$tip.label)
tt3$tip.label <- rename.tomato.tips(tt3$tip.label)
tt4$tip.label <- rename.tomato.tips(tt4$tip.label)
tt5$tip.label <- rename.tomato.tips(tt5$tip.label)
tt6$tip.label <- rename.tomato.tips(tt6$tip.label)
# Add small tips of 0.2 to each branch (ASTRAL does not add terminal branch lengths)
tt1$edge.length[which(is.nan(tt1$edge.length))] <- 0.2
tt1_small$edge.length[which(is.nan(tt1_small$edge.length))] <- 0.2
tt2$edge.length[which(is.nan(tt2$edge.length))] <- 0.2
tt3$edge.length[which(is.nan(tt3$edge.length))] <- 0.2
tt4$edge.length[which(is.nan(tt4$edge.length))] <- 0.2
tt5$edge.length[which(is.nan(tt5$edge.length))] <- 0.2
tt6$edge.length[which(is.nan(tt6$edge.length))] <- 0.2
# Color code clades
tt1_labs <- color.code.tomato.clades(tt1, taxa.numbers = FALSE, trimmed = FALSE)
tt1_small_labs <- color.code.tomato.clades(tt1, taxa.numbers = FALSE, trimmed = TRUE)
# Save plot of the unfiltered tree
# To colour code tip labels without richtext:
#     ggtree(tt1_small) %<+% tt1_small_labs + xlim(NA, 6) + geom_tiplab(aes(label=lab), parse=T)
# To colour code tip labels with richtext:
#    ggtree(tt1_small) %<+% tt1_small_labs + geom_richtext(data=td_filter(isTip), aes(label = name), label.color=NA) + hexpand(.3)
# For more details, see: https://yulab-smu.top/treedata-book/faq.html#faq-formatting-label
p <- ggtree(tt1) %<+% tt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 5) + 
  geom_rootedge(rootedge = 0.3, size = 0.5) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  scale_color_manual(values = tomato_colour_palette) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) +
  theme(axis.text.x = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12), legend.position = c(0.1,0.2))
quilt <-  p + 
  plot_annotation(title = "Tomato dataset - ASTRAL tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
p_name <- paste0(plot_dir, "Tomatoes_ASTRAL_NoTest_plot")
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf", width = 8, height = 7, units = "in")

# Create a small plot of each of the six trees
p1 <- ggtree(tt1_small) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
p2 <- ggtree(tt2) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
p3 <- ggtree(tt3) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
p4 <- ggtree(tt4) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
p5 <- ggtree(tt5) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
p6 <- ggtree(tt6) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = tomato_colour_palette)
# Create a patchwork of the trees tt2-tt7
quilt <- (p1 | p2)/(p3 | p4)/(p5 | p6) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Tomatoes_ASTRAL_Peruvianum_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 13, width = 10, units = "in")

## IQ-Tree Peruvianum topologies
# Need:
#     1. NoTest (P,PHI; P,MaxChi; P,GENECONV; F, All tests)
#     2. F,PHI (identical to F,MaxChi)
#     3. F,Geneconv
#     4. P,All tests
tt1_file <- grep("NoTest", tomatoes_concat_trees, value = TRUE)
tt2_file <- grep("PHI_fail", tomatoes_concat_trees, value = TRUE)
tt3_file <- grep("geneconv_fail", tomatoes_concat_trees, value = TRUE)
tt4_file <- grep("allTests_pass", tomatoes_concat_trees, value = TRUE)
# Open each tree
tt1 <- read.tree(tt1_file)
tt1_small <- read.tree(tt1_file)
tt2 <- read.tree(tt2_file)
tt3 <- read.tree(tt3_file)
tt4 <- read.tree(tt4_file)
# Reroot tree 1 at proper outgroup
tt1 <- root(tt1, outgroup = roots_by_group[["Tomatoes"]])
# Reformat tomato clades for trees 1-4
tt1_small <- reformat.congruent.tomato.clades(tt1_small)
tt2 <- reformat.congruent.tomato.clades(tt2)
tt3 <- reformat.congruent.tomato.clades(tt3)
tt4 <- reformat.congruent.tomato.clades(tt4)
# Rename tip labels to have scientific names (not just numbers)
tt1$tip.label <- rename.tomato.tips(tt1$tip.label)
tt1_small$tip.label <- rename.tomato.tips(tt1_small$tip.label)
tt2$tip.label <- rename.tomato.tips(tt2$tip.label)
tt3$tip.label <- rename.tomato.tips(tt3$tip.label)
tt4$tip.label <- rename.tomato.tips(tt4$tip.label)
# Color code clades
tt1_labs <- color.code.tomato.clades(tt1, taxa.numbers = FALSE, trimmed = FALSE)
tt1_small_labs <- color.code.tomato.clades(tt1, taxa.numbers = FALSE, trimmed = TRUE)
# Save plot of the unfiltered tree
p <- ggtree(tt1) %<+% tt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 5) + 
  geom_rootedge(rootedge = 0.001, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  scale_y_reverse() +  
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = tomato_colour_palette) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) +
  theme(axis.text.x = element_text(size = 15), legend.position = c(0.1,0.2), 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12))
quilt <-  p + 
  plot_annotation(title = "Tomato dataset - Concatenated tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
p_name <- paste0(plot_dir, "Tomatoes_CONCAT_NoTest_plot")
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf", height = 7, width = 8, units = "in")
# Create small plot of each of the four trees
p1 <- ggtree(tt1_small, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = tomato_colour_palette)
p2 <- ggtree(tt2, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = tomato_colour_palette)
p3 <-ggtree(tt3, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = tomato_colour_palette)
p4 <- ggtree(tt4, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = tomato_colour_palette)
# Create a patchwork of the trees tt2-tt7
quilt <- ((p1 | p2)/(p3 | p4)) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Tomatoes_CONCAT_Peruvianum_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 13, width = 10, units = "in")



#### Step 6: Cloudogram of tomato trees ####
# Open file containing all tomato gene trees
all_tomato_gene_trees_file <- paste0(maindir, "tomato_cloudogram/Tomatoes_all_gene_trees.txt")
tomato_gts <- read.tree(all_tomato_gene_trees_file)
# Open the coalescent species tree to be the consensus tree and root it
consensus_tree_file <- paste0(maindir, "species_trees/Tomatoes_NoTest_ASTRAL_species.tre")
consensus_tree <- read.tree(consensus_tree_file)
consensus_tree <- root(consensus_tree, roots_by_group[["Tomatoes"]])
consensus_tree$edge.length[which(is.na(consensus_tree$edge.length))] <- 0.1
# Change tip names
consensus_tree$tip.label <- rename.tomato.tips(consensus_tree$tip.label)
for (i in 1:length(tomato_gts)){
  # Extract tree
  i_tree <- tomato_gts[[i]]
  # Get vector of new tip labels
  i_new_tips <- rename.tomato.tips(i_tree$tip.label)
  # Replace old tips with new tip labels
  i_tree$tip.label <- i_new_tips
  # Replace tree
  tomato_gts[[i]] <- i_tree
}
# Get tip color list
tip_color <- rep("Black", length(consensus_tree$tip.label))
special_tips <- c(15, 18)
tip_color[special_tips] <- "Red"
# Plot densitree of all tomato gene trees
plot_file <- paste0(plot_dir, "Tomatoes_all_gene_trees_densiTree_plot")
# Save as pdf
pdf(file = paste0(plot_file, ".pdf"), width = 12, height = 10)
densiTree(tomato_gts, type = "cladogram", alpha = 0.1, consensus = consensus_tree, scaleX = TRUE, col = "steelblue", cex = 1.2, 
          tip.color = "Black", scale.bar = FALSE)
dev.off()

Cairo::CairoPNG(file = paste0(plot_file, "_small.png"), width = 850, height = 950)
densiTree(tomato_gts, type = "cladogram", alpha = 0.1, consensus = consensus_tree, scaleX = TRUE, col = "steelblue", cex = 1, 
          font = 3, tip.color = "Black", scale.bar = FALSE)
dev.off()


#### Step 7: Plotting Metazoan dataset ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_astral_trees <- grep("ASTRAL", metazoan_tree_files, value = TRUE)
metazoan_concat_trees <- grep("CONCAT", metazoan_tree_files, value = TRUE)
metazoan_palette <- metazoan_palette

## Plot Metazoan ASTRAL tree
# Find files
mt1_file <- grep("NoTest", metazoan_astral_trees, value = TRUE)
# Open trees
mt1 <- read.tree(mt1_file)
# Root tree and drop tips from clades (keep tips within Ctenophora only, keep one species per other clade)
mt1 <- reformat.congruent.metazoan.clades(mt1, trim = "FALSE")
# Add small tips of 0.1 to each branch (ASTRAL does not add terminal branch lengths)
mt1$edge.length[which(is.nan(mt1$edge.length))] <- 0.2
# Relabel clades
mt1_labs <- color.code.metazoan.clades(mt1, trimmed = "FALSE")
# Plot and save
p <- ggtree(mt1) %<+% mt1_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 4) + 
  geom_rootedge(rootedge = 0.3, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(0, 40, 0, 0)) + 
  scale_y_reverse() +  
  scale_color_manual(values = metazoan_colour_palette) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) + 
  theme(axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text (size = 12), legend.position = c(0.08,0.12))
quilt <-  p + 
  plot_annotation(title = "Metazoan dataset - ASTRAL tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
p_name <- paste0(plot_dir, "Metazoan_ASTRAL_NoTest_plot")
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf", width = 12, height = 12, units = "in")

## Plot Metazoan CONCAT tree
# Find files
mt1_file <- grep("NoTest", metazoan_concat_trees, value = TRUE)
# Open trees
mt1 <- read.tree(mt1_file)
# Root tree and drop tips from clades (keep tips within Ctenophora only, keep one species per other clade)
mt1 <- reformat.congruent.metazoan.clades(mt1, trim = "FALSE")
# Relabel clades
mt1_labs <- color.code.metazoan.clades(mt1, trimmed = "FALSE")
# Plot and save
p <- ggtree(mt1) %<+% mt1_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 4) + 
  geom_rootedge(rootedge = 0.05, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(0, 80, 0, 0)) + 
  scale_y_reverse() +  
  scale_color_manual(values = metazoan_colour_palette) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) + 
  theme(axis.text.x = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text (size = 12), legend.position = c(0.08,0.12))
quilt <-  p + 
  plot_annotation(title = "Metazoan dataset - Concatenated tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
p_name <- paste0(plot_dir, "Metazoan_CONCAT_NoTest_plot")
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf", width = 13, height = 12, units = "in")

## Plotting differences in relationships between clades for CONCAT trees
# Find files
mt1_file <- grep("NoTest", metazoan_concat_trees, value = TRUE)
mt2_file <- grep("geneconv", metazoan_concat_trees, value = TRUE)
# Open trees
mt1 <- read.tree(mt1_file)
mt2 <- read.tree(mt2_file)
# Root tree and drop tips from clades (keep tips within Ctenophora only, keep one species per other clade)
mt1 <- reformat.congruent.metazoan.clades(mt1, trim = "Trim_all")
mt2 <- reformat.congruent.metazoan.clades(mt2, trim = "Trim_all")
# Relabel clades
mt1_labs <- color.code.metazoan.clades(mt1, trimmed = "Trim_all")
mt2_labs <- color.code.metazoan.clades(mt2, trimmed = "Trim_all")
# Plot using ggtree
p1 <- ggtree(mt1, size = 0.75) %<+% mt1_labs + 
  geom_tiplab(aes(label = short_lab), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 8) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"))
p2 <- ggtree(mt2, size = 0.75) %<+% mt2_labs + 
  geom_tiplab(aes(label = short_lab), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 8) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"))
quilt <- (p1 | p2) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Metazoan_CONCAT_topology_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf")



#### Step 8: Plotting Plants species trees ####
# Extract tree files
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
plant_tree_files <- grep("Plants", tree_files, value = TRUE)
plant_astral_trees <- grep("ASTRAL", plant_tree_files, value = TRUE)
plant_concat_trees <- grep("CONCAT", plant_tree_files, value = TRUE)

# Open the annotation csv file as a dataframe
annotation_df <- read.csv(annotation_csv_file, stringsAsFactors = FALSE)

## Plants Plot 1: ASTRAL No Test
# Assemble file path and open tree
plants_notest_astral_file <- grep("NoTest", plant_astral_trees, value = TRUE)
p_n_a_tree <- read.tree(plants_notest_astral_file)
# Change edge.length to 0.5
p_n_a_tree <- add.terminal.branches(p_n_a_tree, 0.5)
# Color code clades
plant_labs <- color.plants.by.clades(tree = p_n_a_tree, color_palette = plants_color_palette, clade_df = annotation_df)
# Root tree
chromista_tips <- plant_labs$Code[which(plant_labs$Very.Brief.Classification == "Chromista")]
p_n_a_tree <- root(p_n_a_tree, outgroup = chromista_tips)
# Create plot
p <- ggtree(p_n_a_tree) %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 3, alpha = 1) +
  geom_rootedge(rootedge = 1, linewidth = 0.5) +
  scale_color_manual(values = plants_color_palette) +
  scale_y_reverse() +
  theme(axis.text.x = element_text(size = 12, color = "white"),
        legend.title = element_text(size = 16), legend.text = element_text (size = 13), legend.position = c(0.1,0.3)) +
  guides(color = guide_legend(title = "Clade legend"))
# Create quilt
quilt <-  p + 
  plot_annotation(title = "Plants dataset - ASTRAL tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
# Create plot name
p_name <- paste0(plot_dir, "Plants_ASTRAL_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf")

## Plants Plot 2: CONCAT No Test
# Assemble file path and open tree
plants_notest_concat_file <- grep("NoTest", plant_concat_trees, value = TRUE)
p_n_c_tree <- read.tree(plants_notest_concat_file)
# Color code clades
plant_labs <- color.plants.by.clades(tree = p_n_c_tree, color_palette = plants_color_palette, clade_df = annotation_df)
# Root tree
chromista_tips <- plant_labs$Code[which(plant_labs$Very.Brief.Classification == "Chromista")]
p_n_c_tree <- root(p_n_c_tree, outgroup = chromista_tips)
# Create plot
p <- ggtree(p_n_c_tree) %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 3, alpha = 1) +
  geom_rootedge(rootedge = 0.1, linewidth = 0.5) +
  scale_color_manual(values = plants_color_palette) +
  scale_y_reverse() +
  theme(axis.text.x = element_text(size = 12, color = "white"),
        legend.title = element_text(size = 16), legend.text = element_text (size = 13), legend.position = c(0.1,0.3)) +
  guides(color = guide_legend(title = "Clade legend"))
# Create quilt
quilt <-  p + 
  plot_annotation(title = "Plants dataset - Concatenated tree", subtitle = "Unfiltered dataset",
                  theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
                                plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 0.5)) )
# Create plot name
p_name <- paste0(plot_dir, "Plants_CONCAT_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = quilt, device = "pdf")



#### Step 9: Pretty plotting for Plants comparing deep ASTRAL trees ####
## Supplementary Figure 16
# Get list of trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Find trees
notest_tree_file <- grep("ASTRAL", grep("NoTest", grep("Plants", all_trees, value = TRUE), value = TRUE), value = TRUE)
phi_tree_file <- grep("ASTRAL", grep("PHI_pass", grep("Plants", all_trees, value = TRUE), value = TRUE), value = TRUE)
# Open trees and drop taxa
notest_tree <- read.tree(notest_tree_file)
phi_tree <- read.tree(test_tree_file)
# Extract tips
phi_clade_test <- extract.clade(phi_tree, 1306)
phi_clade_test$edge.length[which(is.na(phi_clade_test$edge.length))] <- 1
# Extract clade
notest_clade <- extract.clade(notest_tree, getMRCA(notest_tree, phi_clade_test$tip.label)) 
phi_clade <- extract.clade(phi_tree, getMRCA(phi_tree, notest_clade$tip.label)) 
# Ladderize clades
notest_clade <- ladderize(notest_clade)
phi_clade <- ladderize(phi_clade)
# Create plant labels
plant_labs <- color.plants.by.clades(tree = notest_tree, color_palette = plants_color_palette, clade_df = annotation_df)
# Trim to only labels in the clade
clade_labs <- plant_labs[which(plant_labs$Code %in% phi_clade$tip.label), ]
# Add arbitrary terminal branch length (ASTRAL does not calculate terminal branch length therefore terminal branches have length of NA)
phi_clade$edge.length[which(is.na(phi_clade$edge.length))] <- 1
notest_clade$edge.length[which(is.na(notest_clade$edge.length))] <- 1
# Change color to red for labs that are involved in branch change
clade_labs$Color_PHI <- "black"
clade_labs$Color_PHI[which(clade_labs$Code %in% phi_clade_test$tip.label)] <- "grey"
clade_labs$Color_PHI[which(clade_labs$Code == "YGAT")] <- "red"
# Plot NoTest ASTRAL clade
notest_plot <- ggtree(notest_clade)  %<+% clade_labs + 
  geom_tiplab(aes(label = Species, color = Color_PHI), size = 4) + 
  geom_rootedge(rootedge = 0.5) +
  scale_color_manual(values = c("red" = "red", "black" = "black")) +
  guides(color = "none") +
  labs(title = "ASTRAL Unfiltered") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  xlim(-0.5, 14)
# Plot P_PHI ASTRAL clade
phi_plot <- ggtree(phi_clade)  %<+% clade_labs + 
  geom_tiplab(aes(label = Species, color = Color_PHI), size = 4) + 
  geom_rootedge(rootedge = 0.5) +
  scale_color_manual(values = c("red" = "red", "black" = "black")) +
  guides(color = "none") +
  labs(title = "ASTRAL P_PHI") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(b = 20))) +
  xlim(-0.5, 15)
# Save tree
quilt <- (notest_plot + phi_plot) + plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Plants_SuppFigure_OutlierBranch_plot")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, width = 15)


