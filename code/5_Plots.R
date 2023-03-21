### gene_filtering/5_Plots.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2022

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
library(ggplot2) # for nice plots
library(ggtree) # for plotting phylogenetic trees
library(ggtext) # for nice tree plots
library(patchwork) # for collating plots
# library(phangorn) # for densiTree function - replaced by ggdensitree
library(TreeTools) # for CollapseNode function

# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "species_trees/")



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
# Color code incongruent clade
pt1_labs <-color.code.primate.clades(pt1, concatenated = FALSE)
# Create plot
p <- ggtree(pt1) %<+% pt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0.002, geom = "text", size = 5) + 
  scale_y_reverse() +  
  scale_x_continuous(breaks = seq(0,15,3)) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) + 
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = c(Variable = "gray50", Congruent = "black"))
# Create plot name
p_name <- paste0(plot_dir, "Primates_ASTRAL_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = p, device = "pdf")

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
p1_0 <- ggtree(pt1, size = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, size = 1) + xlim(-0.1, 1.2)
p2 <- ggtree(pt2, size = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, size = 1) + xlim(-0.1, 1.2)
p3 <- ggtree(pt3, size = 1) + geom_tiplab(offset = 0, geom = "text", size = 10) + 
  geom_rootedge(rootedge = 0.1, size = 1) + xlim(-0.1, 1.2)
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
# Color code incongruent clade
pt1_labs <-color.code.primate.clades(pt1, concatenated = TRUE)
# Create plot
p <- ggtree(pt1) %<+% pt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0.002, geom = "text", size = 5.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) + 
  theme(axis.text.x = element_text(size = 15)) +
  scale_color_manual(values = c(Variable = "gray50", Congruent = "black"))
# Create plot name
p_name <- paste0(plot_dir, "Primates_CONCAT_NoTest_plot")
# Save plot
ggsave(filename = paste0(p_name, ".pdf"), plot = p, device = "pdf")

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



#### Step 4: Plotting Tomatoes dataset ####
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
p1 <- ggtree(tt1) %<+% tt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.3, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Outgroup = "black")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) +
  theme(axis.text.x = element_text(size = 15), legend.position = "left", 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12))
p1_name <- paste0(plot_dir, "Tomatoes_ASTRAL_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf", width = 8, height = 7, units = "in")


# Create a small plot of each of the six trees
p1 <- ggtree(tt1_small) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p2 <- ggtree(tt2) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p3 <- ggtree(tt3) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                vClade_Outgroup = "black"))
p4 <- ggtree(tt4) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p5 <- ggtree(tt5) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p6 <- ggtree(tt6) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.2, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 175, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
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
p1 <- ggtree(tt1) %<+% tt1_labs +
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.001, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Outgroup = "black")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) +
  theme(axis.text.x = element_text(size = 15), legend.position = "left", 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12))
p1_name <- paste0(plot_dir, "Tomatoes_CONCAT_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf", height = 7, width = 8, units = "in")
# Create small plot of each of the four trees
p1 <- ggtree(tt1_small, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Outgroup = "black"))
p2 <- ggtree(tt2, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy",
                                Outgroup = "black"))
p3 <-ggtree(tt3, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Outgroup = "black"))
p4 <- ggtree(tt4, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = F, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white"), legend.position = "bottom") +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Outgroup = "black"))
# Create a patchwork of the trees tt2-tt7
quilt <- ((p1 | p2)/(p3 | p4)) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Tomatoes_CONCAT_Peruvianum_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", height = 13, width = 10, units = "in")



#### Step 5: Plotting Metazoan dataset ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_astral_trees <- grep("ASTRAL", metazoan_tree_files, value = TRUE)
metazoan_concat_trees <- grep("CONCAT", metazoan_tree_files, value = TRUE)

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
p1 <- ggtree(mt1) %<+% mt1_labs +
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 4) + 
  geom_rootedge(rootedge = 0.001, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) + 
  theme(axis.text.x = element_text(size = 15), legend.position = "left", 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12))
p1_name <- paste0(plot_dir, "Metazoan_ASTRAL_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf", width = 8, height = 10, units = "in")

## Plotting differences in topology for ASTRAL trees
# Find files
mt1_file <- grep("NoTest", metazoan_astral_trees, value = TRUE)
mt2_file <- grep("PHI", metazoan_astral_trees, value = TRUE)
mt3_file <- grep("maxchi", metazoan_astral_trees, value = TRUE)
mt4_file <- grep("geneconv", metazoan_astral_trees, value = TRUE)
# Open trees
mt1 <- read.tree(mt1_file)
mt2 <- read.tree(mt2_file)
mt3 <- read.tree(mt3_file)
mt4 <- read.tree(mt4_file)
# Root tree and drop tips from clades (keep tips within Ctenophora only, keep one species per other clade)
mt1 <- reformat.congruent.metazoan.clades(mt1, trim = "Keep_Ctenophora")
mt2 <- reformat.congruent.metazoan.clades(mt2, trim = "Keep_Ctenophora")
mt3 <- reformat.congruent.metazoan.clades(mt3, trim = "Keep_Ctenophora")
mt4 <- reformat.congruent.metazoan.clades(mt4, trim = "Keep_Ctenophora")
# Add small tips of 0.2 to each branch (ASTRAL does not add terminal branch lengths)
mt1$edge.length[which(is.nan(mt1$edge.length))] <- 0.2
mt2$edge.length[which(is.nan(mt2$edge.length))] <- 0.2
mt3$edge.length[which(is.nan(mt3$edge.length))] <- 0.2
mt4$edge.length[which(is.nan(mt4$edge.length))] <- 0.2
# Relabel clades
mt_labs <- color.code.metazoan.clades(mt1, trimmed = "Keep_Ctenophora")
# Plot using ggtree
p1 <- ggtree(mt1, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60"))
p2 <- ggtree(mt2, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60"))
p3 <- ggtree(mt3, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60"))
p4 <- ggtree(mt4, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60"))
quilt <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Metazoan_ASTRAL_Ctenophora_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", width = 8, height = 10, units = "in")

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
p1 <- ggtree(mt1) %<+% mt1_labs +
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = T, offset = 0, geom = "text", size = 4) + 
  geom_rootedge(rootedge = 0.001, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Choanoflagellata = "gray60")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp."))) + 
  theme(axis.text.x = element_text(size = 15), legend.position = "left", 
        legend.title = element_text(size = 15), legend.text = element_text (size = 12))
p1_name <- paste0(plot_dir, "Metazoan_CONCAT_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf", width = 8, height = 10, units = "in")

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

## Plotting differences in Ctenophora topology for CONCAT trees
# Find files
mt1_file <- grep("NoTest", metazoan_concat_trees, value = TRUE)
mt2_file <- grep("PHI", metazoan_concat_trees, value = TRUE)
mt3_file <- grep("maxchi", metazoan_concat_trees, value = TRUE)
mt4_file <- grep("geneconv", metazoan_concat_trees, value = TRUE)
# Open trees
mt1 <- read.tree(mt1_file)
mt2 <- read.tree(mt2_file)
mt3 <- read.tree(mt3_file)
mt4 <- read.tree(mt4_file)
# Root tree and drop tips from clades (keep tips within Ctenophora only, keep one species per other clade)
mt1 <- reformat.congruent.metazoan.clades(mt1, trim = "Keep_Ctenophora")
mt2 <- reformat.congruent.metazoan.clades(mt2, trim = "Keep_Ctenophora")
mt3 <- reformat.congruent.metazoan.clades(mt3, trim = "Keep_Ctenophora")
mt4 <- reformat.congruent.metazoan.clades(mt4, trim = "Keep_Ctenophora")
# Relabel clades
mt_labs <- color.code.metazoan.clades(mt1, trimmed = "Keep_Ctenophora")
# Plot using ggtree
p1 <- ggtree(mt1, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "gray60"))
p2 <- ggtree(mt2, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "gray60"))
p3 <- ggtree(mt3, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "gray60"))
p4 <- ggtree(mt4, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "gray60"))
quilt <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Metazoan_CONCAT_Ctenophora_comparison_plots")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", width = 8, height = 10, units = "in")



#### Step 6: Plotting differences in Primate trees for supplementary data ####
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
p1 <- textGrob(c1_labs$clean_taxa_names[[1]], gp = gpar(fontface = "italic", size = 12))

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
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf",)



#### Step 7: Pretty plotting for Plants comparing deep ASTRAL trees ####
## Distinct edges 2:8 - movement of taxa YGAT
## Supplementary Figure 14

# Get list of trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Find trees
notest_tree_file <- grep("ASTRAL", grep("NoTest", grep("Plants", all_trees, value = TRUE), value = TRUE), value = TRUE)
test_tree_file <- grep("ASTRAL", grep("PHI_pass", grep("Plants", all_trees, value = TRUE), value = TRUE), value = TRUE)
# Open trees and drop taxa
notest_tree <- read.tree(notest_tree_file)
test_tree <- read.tree(test_tree_file)
# Get list of taxa involved
keep_tips <- c("YGAT", "VNMY", "RVGH", "ZTLR", "RPPC", "AXPJ", "KCPT", "TXMP", "VXOD", "HNCF", "BVOF", "BHYC", "AEPI", "MYVH", "POZS", 
               "OODC", "XPBC", "HBUQ", "ZBVT", "CKDK", "TVCU","FWCQ", "OSIP", "BNDE", "NFXV", "PXYR", "RHAU", "PAZJ", "VVPY", "XNLP",
               "Manes_v4.1", "NJLF", "LPGY", "COAQ", "EZZT", "SIZE", "ZIWB", "Poptr_v3.0", "INQX", "RZTJ", "LFOG", "TDTF", "GLVK", 
               "IEPQ", "KKDQ")
# Drop unneeded taxa
notest_tree <- keep.tip(notest_tree, keep_tips)
test_tree <- keep.tip(test_tree, keep_tips)
# Add arbitrary terminal branch length (ASTRAL does not calculate terminal branch length therefore terminal branches have length of NA)
notest_tree$edge.length[which(is.na(notest_tree$edge.length))] <- 0.1
test_tree$edge.length[which(is.na(test_tree$edge.length))] <- 0.1
# Use the annotation_csv_file to extract the relevant tips and tip full names/classifications
annotations_df <- read.csv(annotation_csv_file)
lab_df <- annotations_df[which(annotations_df$Code %in% keep_tips),]
lab_df <- lab_df[match(keep_tips, lab_df$Code),]
lab_df$Species[which(lab_df$Species == "Passiflora_sp")] <- "Passiflora sp"
# Use lab_df as labels for the ggtree plots
lab_df$Color <- c(rep("A", 5), rep("B", (length(keep_tips) - 5) ))
lab_df <- dplyr::mutate(lab_df, 
                        lab = glue('italic("{Species}")'))
# Plot tree
p1 <- ggtree(notest_tree, branch.length = "none")  %<+% lab_df + 
  geom_tiplab(aes(label = lab, color = Color), size = 4, parse = T, show.legend = F) + 
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(subset = !isTip, label=label), nudge_x = -0.6, nudge_y = 0.5, color = "Gray60", size = 3) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) +
  theme(axis.text.x = element_text(size = 13, color = "White"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_color_manual(values = c(A = "Red", B = "Black"))

p2 <- ggtree(test_tree, branch.length = "none")  %<+% lab_df + 
  geom_tiplab(aes(label = lab, color = Color), size = 4, parse = T, show.legend = F) + 
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(subset = !isTip, label=label), nudge_x = -0.6, nudge_y = 0.5, color = "Gray60", size = 3) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 170, 6, 6)) +
  theme(axis.text.x = element_text(size = 13, color = "White"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_color_manual(values = c(A = "Red", B = "Black"))
# Save tree
quilt <- (p1 + p2) + plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_name <- paste0(plot_dir, "Plants_SuppFigure_VNMY_YGAT_OutlierBranch_plot")
ggsave(filename = paste0(quilt_name, ".pdf"), plot = quilt, device = "pdf", width = 12, height = 10, units = "in")



#### Step 8: Cloudogram of tomato trees ####
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
pdf(file = paste0(plot_file, ".pdf"), width = 10, height = 10)
densiTree(tomato_gts, type = "cladogram", alpha = 0.1, consensus = consensus_tree, scaleX = TRUE, col = "steelblue", cex = 1.2, 
          tip.color = tip_color, scale.bar = FALSE)
dev.off()
# Save as png
png(file = paste0(plot_file, ".png"), width = 900, height = 800)
densiTree(tomato_gts, type = "cladogram", alpha = 0.1, consensus = consensus_tree, scaleX = TRUE, col = "steelblue", cex = 1.5, 
          tip.color = tip_color, scale.bar = FALSE)
dev.off()



#### Step 9: Plotting key differences between all trees ####
## ggdensitree for Primates dataset ##
# Get all species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract trees for that combination of dataset and tree method
plot_tree_files <- grep("Primates", all_trees, value = TRUE)
astral_trees_files <- grep("ASTRAL", plot_tree_files, value = TRUE)
concat_trees_files <- grep("CONCAT", plot_tree_files, value = TRUE)
# Extract text file for each tree
astral_trees_text <- unlist(lapply(astral_trees_files, readLines))
concat_trees_text <- unlist(lapply(concat_trees_files, readLines))
# Read trees into a multiphylo object
astral_trees <- read.tree(text = astral_trees_text)
concat_trees <- read.tree(text = concat_trees_text)
# Extract the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree_file <- grep("NoTest", astral_trees_files, value = TRUE)
notest_concat_tree_file <- grep("NoTest", concat_trees_files, value = TRUE)
# Open the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree <- read.tree(notest_astral_tree_file)
notest_concat_tree <- read.tree(notest_concat_tree_file)
# Root trees
astral_trees <- lapply(astral_trees, root, roots_by_group[["Primates"]])
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
concat_trees <- lapply(concat_trees, root, roots_by_group[["Primates"]])
class(concat_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- root(notest_astral_tree, roots_by_group[["Primates"]])
notest_concat_tree <- root(notest_concat_tree, roots_by_group[["Primates"]])
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
astral_trees <- lapply(1:length(astral_trees), 
                       function(i){astral_trees[[i]] <- reformat.ASTRAL.tree.for.plotting(astral_trees[[i]], 
                                                                                          add.arbitrary.terminal.branches = TRUE, 
                                                                                          terminal.branch.length = 1, 
                                                                                          strip.nodes = FALSE,
                                                                                          scale.tree.length = FALSE, 
                                                                                          new.tree.length = NA)} )
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- reformat.ASTRAL.tree.for.plotting(notest_astral_tree, 
                                                        add.arbitrary.terminal.branches = TRUE, 
                                                        terminal.branch.length = 1, 
                                                        strip.nodes = FALSE,
                                                        scale.tree.length = FALSE, 
                                                        new.tree.length = NA)
# Get the order for the tips (bottom species first, top species last)
primates_tip_order <- c("Pan troglodytes", "Pan paniscus", "Homo sapiens", "Gorilla gorilla",
                        "Pongo abelii", "Nomascus leucogenys", "Macaca mulatta", "Macaca fascicularis",
                        "Macaca nemestrina", "Theropithecus gelada", "Papio anubis", "Mandrillus leucophaeus",
                        "Cercocebus atys", "Chlorocebus sabaeus", "Rhinopithecus roxellana", "Rhinopithecus bieti",
                        "Piliocolobus tephrosceles", "Colobus angolensis palliatus", "Saimiri boliviensis", "Cebus capucinus imitator",
                        "Aotus nancymaae", "Callithrix jacchus", "Carlito syrichta", "Microcebus murinus", 
                        "Propithecus coquereli", "Otolemur garnettii", "Galeopterus variegatus", "Tupaia chinensis",
                        "Mus musculus")
# Create labels for the tips
tip_labels_df <- color.code.primate.clades(notest_concat_tree, color = FALSE)
# Reorder tip_labels_df to make sure ggtree will be in desired order
labs <- tip_labels_df[match(primates_tip_order, tip_labels_df$taxa),]
# Format taxa names
labs$taxa <- gsub(" ", "_", labs$taxa)
# Plot a nice densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = labs$taxa, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.2) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Plot a nice densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = labs$taxa, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.2) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Primates_Species_tree_comparison_ggdensitree")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 10, height = 8, units = "in")
# Plot a nice annotated densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = labs$taxa, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.2) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 42, label = "Hominidae", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "deepskyblue", barcolor = "deepskyblue")+
  geom_cladelab(node = 12, label = "Hylobatidae", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "steelblue3", barcolor = "steelblue3") +
  geom_cladelab(node = 50, label = "Cercopithecinae", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "darkslategray4", barcolor = "darkslategray4") +
  geom_cladelab(node = 47, label = "Colobinae", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "darkslategray3", barcolor = "darkslategray3") +
  geom_cladelab(node = 32, label = "Cebidae", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "darkseagreen4", barcolor = "darkseagreen4") +
  geom_cladelab(node = 5, label = "Tarsiiformes", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "springgreen4", barcolor = "springgreen4") +
  geom_cladelab(node = 39, label = "Lemuriformes", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "olivedrab4", barcolor = "olivedrab4") +
  geom_cladelab(node = 9, label = "Lorisiformes", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "olivedrab3", barcolor = "olivedrab3") +
  geom_cladelab(node = 6, label = "Dermoptera", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4") +
  geom_cladelab(node = 8, label = "Scandentia", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4") +
  geom_cladelab(node = 7, label = "Rodentia", align = TRUE, offset = 10, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4")
# Plot a nice annotated densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = labs$taxa, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.2) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 52, label = "Hominidae", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "deepskyblue", barcolor = "deepskyblue")+
  geom_cladelab(node = 27, label = "Hylobatidae", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "steelblue3", barcolor = "steelblue3") +
  geom_cladelab(node = 41, label = "Cercopithecinae", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "darkslategray4", barcolor = "darkslategray4") +
  geom_cladelab(node = 48, label = "Colobinae", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "darkslategray3", barcolor = "darkslategray3") +
  geom_cladelab(node = 32, label = "Cebidae", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "darkseagreen4", barcolor = "darkseagreen4") +
  geom_cladelab(node = 3, label = "Tarsiiformes", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "springgreen4", barcolor = "springgreen4") +
  geom_cladelab(node = 38, label = "Lemuriformes", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "olivedrab4", barcolor = "olivedrab4") +
  geom_cladelab(node = 9, label = "Lorisiformes", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "olivedrab3", barcolor = "olivedrab3") +
  geom_cladelab(node = 4, label = "Dermoptera", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4") +
  geom_cladelab(node = 6, label = "Scandentia", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4") +
  geom_cladelab(node = 5, label = "Rodentia", align = TRUE, offset = 11, offset.text = 0.2, fontsize = 4, textcolor = "snow4", barcolor = "snow4")
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Primates_Species_tree_comparison_ggdensitree_annotated")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 14, height = 8, units = "in")


## ggdensitree for Tomatoes dataset ##
# Get all species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract trees for that combination of dataset and tree method
plot_tree_files <- grep("Tomatoes", all_trees, value = TRUE)
astral_trees_files <- grep("ASTRAL", plot_tree_files, value = TRUE)
concat_trees_files <- grep("CONCAT", plot_tree_files, value = TRUE)
# Extract text file for each tree
astral_trees_text <- unlist(lapply(astral_trees_files, readLines))
concat_trees_text <- unlist(lapply(concat_trees_files, readLines))
# Read trees into a multiphylo object
astral_trees <- read.tree(text = astral_trees_text)
concat_trees <- read.tree(text = concat_trees_text)
# Extract the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree_file <- grep("NoTest", astral_trees_files, value = TRUE)
notest_concat_tree_file <- grep("NoTest", concat_trees_files, value = TRUE)
# Open the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree <- read.tree(notest_astral_tree_file)
notest_concat_tree <- read.tree(notest_concat_tree_file)
# Root trees
astral_trees <- lapply(astral_trees, root, roots_by_group[["Tomatoes"]])
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
concat_trees <- lapply(concat_trees, root, roots_by_group[["Tomatoes"]])
class(concat_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- root(notest_astral_tree, roots_by_group[["Tomatoes"]])
notest_concat_tree <- root(notest_concat_tree, roots_by_group[["Tomatoes"]])
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
astral_trees <- lapply(1:length(astral_trees), 
                       function(i){astral_trees[[i]] <- reformat.ASTRAL.tree.for.plotting(astral_trees[[i]], 
                                                                                          add.arbitrary.terminal.branches = TRUE, 
                                                                                          terminal.branch.length = 1, 
                                                                                          strip.nodes = FALSE,
                                                                                          scale.tree.length = FALSE, 
                                                                                          new.tree.length = NA)} )
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- reformat.ASTRAL.tree.for.plotting(notest_astral_tree, 
                                                        add.arbitrary.terminal.branches = TRUE, 
                                                        terminal.branch.length = 1, 
                                                        strip.nodes = FALSE,
                                                        scale.tree.length = FALSE, 
                                                        new.tree.length = NA)
# Get the order for the tips (bottom species first, top species last)
tomato_tip_order <- c("LA3909", "LA0436", "LA0429", "LA3124", "LA3475", "SL2.50", "LA1589", "LA1269",
                      "LA2933", "LA2133", "LA1322", "LA2172", "LA1316", "LA1028", "LA1364", "LA1782",
                      "LA4117", "LA2744", "LA2964", "LA1358", "LA0107", "LA0444", "LA1777", "LA0407",
                      "LA3778", "LA0716", "LA4126", "LA2951", "LA4116")
tomato_species_order <- rename.tomato.tips(tomato_tip_order)
# Create labels for the tips
tip_labels_df <- color.code.tomato.clades(notest_astral_tree, taxa.numbers = FALSE, trimmed = FALSE, color = FALSE)
# Reorder tip_labels_df to make sure ggtree will be in desired order
labs <- tip_labels_df[match(tomato_species_order, tip_labels_df$taxa),]
# Replace taxa names with taxa numbers
labs$taxa <- tomato_tip_order
# Plot a nice densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = tomato_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.5) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Plot a nice densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = tomato_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.5) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Tomatoes_Species_tree_comparison_ggdensitree")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 10, height = 8, units = "in")
# Plot a nice annotated densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = tomato_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.5) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 110, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 57, label = "sect.\nLycopersicoides", align = TRUE, horizontal = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "black", barcolor = "black") +
  geom_cladelab(node = 47, label = "Hirsutum", align = TRUE, horizontal = FALSE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "navy", barcolor = "navy") +
  geom_cladelab(node = 50, label = "Peruvianum", align = TRUE, horizontal = FALSE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "darkgreen", barcolor = "darkgreen") +
  geom_cladelab(node = 40, label = "Arcanum", align = TRUE, horizontal = FALSE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "goldenrod3", barcolor = "goldenrod3") +
  geom_cladelab(node = 36, label = "Esculentum", align = TRUE, horizontal = FALSE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "firebrick3", barcolor = "firebrick3")
# Plot a nice annotated densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = tomato_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.5) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 110, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 57, label = "sect.\nLycopersicoides", align = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "black", barcolor = "black") +
  geom_cladelab(node = 44, label = "Hirsutum", align = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "navy", barcolor = "navy") +
  geom_cladelab(node = 36, label = "Peruvianum", align = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "darkgreen", barcolor = "darkgreen") +
  geom_cladelab(node = 48, label = "Arcanum", align = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "goldenrod3", barcolor = "goldenrod3") +
  geom_cladelab(node = 33, label = "Esculentum", align = TRUE, offset = 6.5, offset.text = 0.2, fontsize = 6, textcolor = "firebrick3", barcolor = "firebrick3")
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Tomatoes_Species_tree_comparison_ggdensitree_annotated")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 16, height = 8, units = "in")


## ggdensitree for Metazoans dataset ##
# Get all species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract trees for that combination of dataset and tree method
plot_tree_files <- grep("Metazoan", all_trees, value = TRUE)
astral_trees_files <- grep("ASTRAL", plot_tree_files, value = TRUE)
concat_trees_files <- grep("CONCAT", plot_tree_files, value = TRUE)
# Extract text file for each tree
astral_trees_text <- unlist(lapply(astral_trees_files, readLines))
concat_trees_text <- unlist(lapply(concat_trees_files, readLines))
# Read trees into a multiphylo object
astral_trees <- read.tree(text = astral_trees_text)
concat_trees <- read.tree(text = concat_trees_text)
# Extract the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree_file <- grep("NoTest", astral_trees_files, value = TRUE)
notest_concat_tree_file <- grep("NoTest", concat_trees_files, value = TRUE)
# Open the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree <- read.tree(notest_astral_tree_file)
notest_concat_tree <- read.tree(notest_concat_tree_file)
# Root trees
astral_trees <- lapply(astral_trees, root, roots_by_group[["Metazoan"]])
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
concat_trees <- lapply(concat_trees, root, roots_by_group[["Metazoan"]])
class(concat_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- root(notest_astral_tree, roots_by_group[["Metazoan"]])
notest_concat_tree <- root(notest_concat_tree, roots_by_group[["Metazoan"]])
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
astral_trees <- lapply(1:length(astral_trees), 
                       function(i){astral_trees[[i]] <- reformat.ASTRAL.tree.for.plotting(astral_trees[[i]], 
                                                                                          add.arbitrary.terminal.branches = TRUE, 
                                                                                          terminal.branch.length = 1, 
                                                                                          strip.nodes = FALSE,
                                                                                          scale.tree.length = FALSE, 
                                                                                          new.tree.length = NA)} )
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- reformat.ASTRAL.tree.for.plotting(notest_astral_tree, 
                                                        add.arbitrary.terminal.branches = TRUE, 
                                                        terminal.branch.length = 1, 
                                                        strip.nodes = FALSE,
                                                        scale.tree.length = FALSE, 
                                                        new.tree.length = NA)
# Get the order for the tips (bottom species first, top species last)
metazoans_tip_order <- c("Monosiga_ovata", "Acanthoeca_sp", "Monosiga_brevicolis", "Salpingoeca_rosetta", "Salpingoeca_pyxidium",
                         "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata", "Beroe_forskalii", "Beroe_sp_Queensland_Australia",
                         "Lobata_sp_Punta_Arenas_Argentina", "Bolinopsis_ashleyi", "Ctenophora_sp_Florida_USA", "Cestum_veneris",
                         "Eurhamphaea_vexilligera", "Mnemiopsis_leidyi", "Bolinopsis_infundibulum", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA", 
                         "Ocyropsis_sp_Bimini_Bahamas", "Lobatolampea_tetragona", "Dryodora_glandiformis", "Cydippida_sp", "Mertensiidae_sp_Washington_USA",
                         "Mertensiidae_sp_Antarctica", "Callianira_Antarctica", "Cydippida_sp_Maryland_USA", "Pleurobrachia_sp_South_Carolina_USA",
                         "Pleurobrachia_bachei", "Pleurobrachia_pileus", "Hormiphora_californica", "Hormiphora_palmata", "Coeloplana_astericola",
                         "Vallicula_sp", "Euplokamis_dunlapae",
                         "Amphimedon_queenslandica", "Petrosia_ficiformis", "Crella_elegans", "Kirkpatrickia_variolosa", "Latrunculia_apicalis",
                         "Mycale_phylophylla", "Pseudospongosorites_suberitoides", "Cliona_varians", "Spongilla_lacustris", "Chondrilla_nucula",
                         "Ircinia_fasciculata", "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Hyalonema_populiferum",
                         "Corticium_candelabrum", "Oscarella_carmela", "Sycon_ciliatum", "Sycon_coactum",
                         "Trichoplax_adhaerens",
                         "Nanomia_bijuga", "Agalma_elegans", "Abylopsis_tetragona", "Craseo_lathetica", "Physalia_physalia", "Hydra_oligactis",
                         "Hydra_vulgaris", "Hydra_viridissima", "Periphyla_periphyla", "Aiptasia_pallida", "Hormathia_digitata", "Bolocera_tuediae",
                         "Nematostella_vectensis", "Acropora_digitifera", "Eunicella_verrucosa",
                         "Capitella_teleta", "Hemithris_psittacea", "Drosophila_melanogaster", "Daphnia_pulex", "Homo_sapiens", "Strongylocentrotus_purpatus")
# Create labels for the tips
tip_labels_df <- color.code.metazoan.clades(notest_concat_tree, trimmed = "FALSE", color = FALSE)
# Reorder tip_labels_df to make sure ggtree will be in desired order
labs <- tip_labels_df[match(metazoans_tip_order, tip_labels_df$taxa),]
# Plot a nice densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = metazoans_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = long_lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Plot a nice densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = metazoans_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = long_lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 180, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Metazoan_Species_tree_comparison_ggdensitree")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 14, height = 14, units = "in")
# Plot a nice annotated densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, tip.order = metazoans_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = short_lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "ASTRAL species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 151, label = "Choanoflagellata", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = 122, label = "Ctenophora", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "navy", barcolor = "navy") +
  geom_cladelab(node = 100, label = "Porifera", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "goldenrod", barcolor = "goldenrod") +
  geom_cladelab(node = 83, label = "Cnidaria", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "firebrick", barcolor = "firebrick") +
  geom_cladelab(node = 93, label = "Bilateria", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "black", barcolor = "black") +
  geom_cladelab(node = 22, label = "Placozoa", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "darkgreen", barcolor = "darkgreen")
# Plot a nice annotated densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, tip.order = metazoans_tip_order, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") %<+% labs +
  geom_tiplab(aes(label = short_lab), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold")) +
  geom_cladelab(node = 151, label = "Choanoflagellata", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = 114, label = "Ctenophora", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "navy", barcolor = "navy") +
  geom_cladelab(node = 90, label = "Porifera", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "goldenrod", barcolor = "goldenrod") +
  geom_cladelab(node = 82, label = "Cnidaria", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "firebrick", barcolor = "firebrick") +
  geom_cladelab(node = 84, label = "Bilateria", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "black", barcolor = "black") +
  geom_cladelab(node = 62, label = "Placozoa", align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "darkgreen", barcolor = "darkgreen")
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Metazoan_Species_tree_comparison_ggdensitree_annotated")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 16.5, height = 12, units = "in")


## ggdensitree for Plants dataset ##
# Open the annotations file for the plants dataset
annotations_df <- read.csv(annotation_csv_file)
# Extract Chromista taxa: these will be the outgroups
outgroup_taxa <- annotations_df[annotations_df$Very.Brief.Classification == "Chromista ",]$Code

# Get all species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract trees for that combination of dataset and tree method
plot_tree_files <- grep("Plants", all_trees, value = TRUE)
astral_trees_files <- grep("ASTRAL", plot_tree_files, value = TRUE)
concat_trees_files <- grep("CONCAT", plot_tree_files, value = TRUE)
# Extract text file for each tree
astral_trees_text <- unlist(lapply(astral_trees_files, readLines))
concat_trees_text <- unlist(lapply(concat_trees_files, readLines))
# Read trees into a multiphylo object
astral_trees <- read.tree(text = astral_trees_text)
concat_trees <- read.tree(text = concat_trees_text)
# Extract the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree_file <- grep("NoTest", astral_trees_files, value = TRUE)
notest_concat_tree_file <- grep("NoTest", concat_trees_files, value = TRUE)
# Open the NoTest tree (the tree estimated from the unfiltered set of loci)
notest_astral_tree <- read.tree(notest_astral_tree_file)
notest_concat_tree <- read.tree(notest_concat_tree_file)
# Trim the outgroup - get rid of any taxa that aren't in the trees
notest_astral_tree_outgroup <- outgroup_taxa[outgroup_taxa %in% notest_astral_tree$tip.label]
notest_concat_tree_outgroup <- outgroup_taxa[outgroup_taxa %in% notest_concat_tree$tip.label]
# Root multiphylo trees using the trimmed outgroups
notest_astral_tree <- root(notest_astral_tree, notest_astral_tree_outgroup)
notest_concat_tree <- root(notest_concat_tree, notest_concat_tree_outgroup)
astral_trees <- lapply(astral_trees, root, notest_astral_tree_outgroup)
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
concat_trees <- lapply(concat_trees, root, notest_concat_tree_outgroup)
class(concat_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
astral_trees <- lapply(1:length(astral_trees), 
                       function(i){astral_trees[[i]] <- reformat.ASTRAL.tree.for.plotting(astral_trees[[i]], 
                                                                                          add.arbitrary.terminal.branches = TRUE, 
                                                                                          terminal.branch.length = 1, 
                                                                                          strip.nodes = FALSE,
                                                                                          scale.tree.length = FALSE, 
                                                                                          new.tree.length = NA)} )
class(astral_trees) <- "multiPhylo" # change object class from "list" into "multiPhylo
notest_astral_tree <- reformat.ASTRAL.tree.for.plotting(notest_astral_tree, 
                                                        add.arbitrary.terminal.branches = TRUE, 
                                                        terminal.branch.length = 1, 
                                                        strip.nodes = FALSE,
                                                        scale.tree.length = FALSE, 
                                                        new.tree.length = NA)

# Create labels for the tips
tip_labels_df <- label.plant.taxa(notest_concat_tree$tip.label, annotations_df)
# Plot a nice densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") +
  labs(title = "ASTRAL species trees")
# Plot a nice densitree of the concatenated species trees
concat_densitree <- ggdensitree(concat_trees, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") +
  labs(title = "Concatenated species trees")
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Plants_Species_tree_comparison_ggdensitree")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 14, height = 14, units = "in")

# Replot but with clade annotations added
# Identify unique clades
unique_clades <- unique(annotations_df$Very.Brief.Classification)
# Extract most recent commmon ancestor (mrca) for each clade
astral_mrca_nodes <- lapply(unique_clades, 
                            function(c){
                              findMRCA(notest_astral_tree, annotations_df[(annotations_df$Code %in% notest_astral_tree$tip.label & annotations_df$Very.Brief.Classification == c), ]$Code)
                            } )
concat_mrca_nodes <- lapply(unique_clades, 
                            function(c){
                              findMRCA(notest_concat_tree, annotations_df[(annotations_df$Code %in% notest_concat_tree$tip.label & annotations_df$Very.Brief.Classification == c), ]$Code)
                            } )
# Format mrca results
astral_mrca_nodes[sapply(astral_mrca_nodes, is.null)] <- NA
astral_mrca_nodes <- unlist(astral_mrca_nodes)
concat_mrca_nodes[sapply(concat_mrca_nodes, is.null)] <- NA
concat_mrca_nodes <- unlist(concat_mrca_nodes)
# Format unique_clades to remove underscores
mrca_labels <- gsub("_", " ", unique_clades)
# Create dataframe
mrca_df <- data.frame(astral_nodes = mrca_nodes, concat_nodes = concat_mrca_nodes, labels = mrca_labels)
# Remove clades without a most common recent ancestor
mrca_df <- mrca_df[(is.na(mrca_df$concat_nodes) == FALSE), ]
# Remove ANAGrade labels
mrca_df[(mrca_df$labels == "ANAGrade"), ]$labels <- "ANA lineages"
# Plot the astral nodes manually
astral_nodes = c("Streptophyte_algae" = NA, "Chlorophyta" = NA, "Chromista " = 1179, 
                 "Gymnos" = 1872, "ANAGrade" = 1865, "Eudicots" = 1700, "Lycophytes" = 2027, 
                 "Monocots" = 1757, "Mosses" = 2081, "Glaucophyta " = 2174, "Ceratophyllales" = 550,
                 "Chloranthales" = 1729, "Monilophytes" = 1954, "Hornworts" = 2050, 
                 "Liverworts " = 2060, "Magnoliids" = 1733, "Rhodophyta " = 2178)
# Plot a nice densitree of the astral species trees
astral_densitree <- ggdensitree(astral_trees, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") +
  labs(title = "ASTRAL species trees") +
  geom_cladelab(node = astral_nodes[3], label = astral_nodes[3], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[4], label = astral_nodes[4], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[5], label = astral_nodes[5], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[6], label = astral_nodes[6], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[7], label = astral_nodes[7], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[8], label = astral_nodes[8], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[9], label = astral_nodes[9], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[10], label = astral_nodes[10], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[11], label = astral_nodes[11], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[12], label = astral_nodes[12], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[13], label = astral_nodes[13], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[14], label = astral_nodes[14], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[15], label = astral_nodes[15], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[16], label = astral_nodes[16], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = astral_nodes[17], label = astral_nodes[17], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))

# Output test
test_name <- paste0(plot_dir, "test_nodes")
# Assemble the figure
ggsave(filename = paste0(test_name, ".pdf"), plot = astral_densitree, device = "pdf")

# Plot a nice densitree of the concatenated species trees
concat_densitree <- ggdensitree(astral_trees, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue") +
  labs(title = "Concatenated species trees") +
  geom_cladelab(node = mrca_df$concat_nodes[1], label = mrca_df$labels[1], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[2], label = mrca_df$labels[2], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[3], label = mrca_df$labels[3], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[4], label = mrca_df$labels[4], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[5], label = mrca_df$labels[5], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[6], label = mrca_df$labels[6], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[7], label = mrca_df$labels[7], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[8], label = mrca_df$labels[8], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[9], label = mrca_df$labels[9], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[10], label = mrca_df$labels[10], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[11], label = mrca_df$labels[11], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[12], label = mrca_df$labels[12], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[13], label = mrca_df$labels[13], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[14], label = mrca_df$labels[14], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[15], label = mrca_df$labels[15], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  geom_cladelab(node = mrca_df$concat_nodes[16], label = mrca_df$labels[16], align = TRUE, offset = 8.5, offset.text = 0.2, fontsize = 4, textcolor = "gray50", barcolor = "gray50") +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  labs(title = "Concatenated species trees") +
  theme(axis.text.x = element_text(color = "white"), axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        plot.title = element_text(hjust = 0, size = 14, face = "bold"))
# Construct file name for this densitree plot
densitree_name <- paste0(plot_dir, "Plants_Species_tree_comparison_ggdensitree_annotated")
# Assemble the figure
quilt <- (astral_densitree + concat_densitree) + 
  plot_annotation(tag_levels = "a", tag_suffix = ".") & theme(plot.tag = element_text(size = 20))
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 14, height = 8, units = "in")


# Test to get nodes
test_tree <- ggdensitree(astral_trees, align.tips = TRUE, branch.length = "none", alpha = 0.5, color = "steelblue")  %<+% tip_labels_df +
  labs(title = "ASTRAL species trees") +
  geom_tiplab(aes(label = clade, color = clade), parse = FALSE, show.legend = FALSE, offset = 0.2, geom = "text", size = 4) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 80, 6, 6)) +
  theme(axis.text.x = element_text(size = 13, color = "White"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()) +
  scale_color_manual(values = c("Streptophyte algae" = "White", "Chlorophyta" = "white", "Chromista " = "red",
                                "Gymnos" = "white", "ANAGrade" = "red", "Eudicots" = "white", "Lycophytes" = "white", "Monocots" = "white",
                                "Mosses" = "white", "Glaucophyta " = "white", "Ceratophyllales" = "white", "Chloranthales" = "white", 
                                "Monilophytes" = "white", "Hornworts" = "white", "Liverworts " = "white", "Magnoliids" = "white", 
                                "Rhodophyta " = "white")) +
  geom_text(aes(label=node), hjust=-.3, fill = "lightgreen", label.size = 0.5)
test_name <- paste0(plot_dir, "test_nodes")
# Assemble the figure
ggsave(filename = paste0(test_name, ".pdf"), plot = test_tree, device = "pdf", width = 20, height = 49, units = "in")


is.monophyletic(notest_astral_tree, tip_labels_df[tip_labels_df$clade == "Euglenozoa ",]$code)
tree_df <- fortify(notest_astral_tree)


