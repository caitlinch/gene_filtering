### gene_filtering/5_Plots.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2022

## This script:
# 1. Creates a variety of plots and figures



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                 <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir                <- for saving plots and analyses.
# datasets                <- set name(s) for the dataset(s)
# roots                   <- set which taxa is outgroup for each dataset

### Caitlin's paths ###
# Folders and filepaths
# # For work computer:
# maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
# plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
# For laptop:
maindir <- "/Users/caitlin/Repositories/gene_filtering/"
plot_dir <- "/Users/caitlin/Documents/PhD/Ch01_EmpiricalTreelikeness/plots/"

# Dataset information
datasets <- c("Vanderpool2020", "Pease2016", "Whelan2017", "1KP")
roots <- list("1KP" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                        "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                        "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
              "Whelan2017" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
              "Vanderpool2020" = c("Mus_musculus"), 
              "Pease2016" = c("LA4116", "LA2951", "LA4126"))
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ape) # functions: read.tree, Ntip, root
library(ggplot2) # for nice plots
library(ggtree) # for plotting phylogenetic trees
library(ggtext) # for nice tree plots
library(patchwork) # for collating plots
# Source functions
source(paste0(maindir,"code/func_plots.R"))
# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "trees/")



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
pt1 <- root(pt1, outgroup = roots[["Vanderpool2020"]])
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
pt1 <- root(pt1, outgroup = roots[["Vanderpool2020"]])
pt2 <- root(pt2, outgroup = roots[["Vanderpool2020"]])
pt3 <- root(pt3, outgroup = roots[["Vanderpool2020"]])
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
pt1 <- root(pt1, outgroup = roots[["Vanderpool2020"]])
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
pt1 <- root(pt1, outgroup = roots[["Vanderpool2020"]])
pt2 <- root(pt2, outgroup = roots[["Vanderpool2020"]])
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
tt1 <- root(tt1, outgroup = roots[["Pease2016"]])
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
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.3, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 150, 6, 6)) + 
  theme(axis.text.x = element_text(size = 15)) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p1_name <- paste0(plot_dir, "Tomatoes_ASTRAL_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf")
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
tt1 <- root(tt1, outgroup = roots[["Pease2016"]])
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
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.001, size = 0.5) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 12)) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p1_name <- paste0(plot_dir, "Tomatoes_CONCAT_NoTest_plot")
ggsave(filename = paste0(p1_name, ".pdf"), plot = p1, device = "pdf")
# Create small plot of each of the four trees
p1 <- ggtree(tt1_small, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p2 <- ggtree(tt2, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy",
                                Clade_Outgroup = "black"))
p3 <-ggtree(tt3, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
p4 <- ggtree(tt4, size = 0.75) %<+% tt1_small_labs + 
  geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 6) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Esculentum = "firebrick3", Arcanum = "goldenrod3", 
                                Peruvianum = "darkgreen", Hirsutum = "navy", 
                                Clade_Outgroup = "black"))
# Create a patchwork of the trees tt2-tt7
quilt <- (p1 | p2)/(p3 | p4) +
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


## Plotting differences in topology for CONCAT trees
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
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "black"))
p2 <- ggtree(mt2, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "black"))
p3 <- ggtree(mt3, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "black"))
p4 <- ggtree(mt4, size = 0.4) %<+% mt_labs + 
  geom_tiplab(aes(label = short_lab, color = clade), parse=T, show.legend = FALSE, offset = 0, geom = "text", size = 3) + 
  geom_rootedge(rootedge = 0.005, size = 0.5) + 
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 160, 6, 6)) + 
  theme(axis.text.x = element_text(size = 0), axis.line.x = element_line(colour = "white"), 
        axis.ticks.x = element_line(colour = "white")) +
  scale_color_manual(values = c(Bilateria = "black", Cnidaria = "firebrick3", Placozoa = "darkgreen", 
                                Porifera = "goldenrod3", Ctenophora = "navy", Clade_Outgroup = "black"))
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


## Plotting differences in topology for CONCAT trees
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

## Plotting differences in topology for CONCAT trees
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


