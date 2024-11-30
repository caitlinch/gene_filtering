### gene_filtering/code/5_Densitrees.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2023

## This script creates a densitree plots for the 4 empirical datasets

##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                 <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir                <- for saving plots and analyses.
# annotations_csv_file    <- location of misc/annotations.csv file from the Leebens-Mack (2019) "Data from 1000 Plants Transcriptomes" data repository
# roots_by_groups         <- set which taxa is outgroup for each dataset: Primates, Tomatoes, Metazoans, and Plants

### Caitlin's paths ###
# Folders and filepaths
location = "macbook"
if (location == "work"){
  maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
  plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
  annotation_csv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/misc/annotations.csv"
} else if (location == "macbook"){
  maindir <- "/Users/caitlin/Repositories/gene_filtering/"
  plot_dir <- "/Users/caitlin/Documents/PhD/Ch01_EmpiricalTreelikeness/thesis_revision_plots/"
  annotation_csv_file <- "/Users/caitlin/Documents/PhD/Ch01_EmpiricalTreelikeness/annotations.csv"
}

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
library(ggtree) # for plotting phylogenetic trees and densitress (ggdensitree)
library(patchwork)
library(colorBlindness) # For Plants clade labels
#library(ggtext) # for nice tree plots
#library(patchwork) # for collating plots
#library(TreeTools) # for CollapseNode function

# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Save original graphical parameters
reset_graph_params <- par()

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
                             "Porifera" = "#E69F00", "Placozoa" = "#000000", "Outgroup" = "#999999",
                             "Choanoflagellata" = "#999999")
plants_color_palette <- c(SteppedSequential5Steps[c(1,3,5,6,8,10,11,13,15,16,18,20,21,23,25)], "black", "grey40", "grey70")
plant_classifications <-  c("Chromista", "Rhodophyta", "Glaucophyta", "Chlorophyta", "Streptophyte algae",
                            "Hornworts", "Liverworts", "Mosses", "Lycophytes", "Monilophytes", "Gymnos",
                            "ANAGrade", "Monocots", "Magnoliids", "Chloranthales", "Eudicots", "Ceratophyllales")
names(plants_color_palette) <- plant_classifications

# Read in plant annotations
annotation_df <- read.csv(annotation_csv_file)



#### Step 3: Primates dataset ####
## Open Primate dataset trees
# Get all files for species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract files for this dataset
plot_tree_files <- grep("Primates", all_trees, value = TRUE)
# Extract trees for that combination of dataset and tree method
astral_trees_files <- grep("NoTest", grep("ASTRAL", plot_tree_files, value = TRUE), value = T, invert = T)
concat_trees_files <- grep("NoTest", grep("CONCAT", plot_tree_files, value = TRUE), value = T, invert = T)

# Separate into pass and fail trees
astral_pass_tree_files <- grep("pass", astral_trees_files, value = T)
astral_fail_tree_files <- grep("fail", astral_trees_files, value = T)
concat_pass_tree_files <- grep("pass", concat_trees_files, value = T)
concat_fail_tree_files <- grep("fail", concat_trees_files, value = T)
# Extract text file for each tree
astral_pass_trees_text <- unlist(lapply(astral_pass_tree_files, readLines))
astral_fail_trees_text <- unlist(lapply(astral_fail_tree_files, readLines))
concat_pass_trees_text <- unlist(lapply(concat_pass_tree_files, readLines))
concat_fail_trees_text <- unlist(lapply(concat_fail_tree_files, readLines))

## Open trees
# Read trees into a multiphylo object
a_p_trees <- read.tree(text = astral_pass_trees_text)
a_f_trees <- read.tree(text = astral_fail_trees_text)
c_p_trees <- read.tree(text = concat_pass_trees_text)
c_f_trees <- read.tree(text = concat_fail_trees_text)

## Root trees
a_p_trees <- lapply(a_p_trees, root, roots_by_group[["Primates"]])
a_f_trees <- lapply(a_f_trees, root, roots_by_group[["Primates"]])
c_p_trees <- lapply(c_p_trees, root, roots_by_group[["Primates"]])
c_f_trees <- lapply(c_f_trees, root, roots_by_group[["Primates"]])
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(a_f_trees) <- "multiPhylo"
class(c_p_trees) <- "multiPhylo"
class(c_f_trees) <- "multiPhylo"

## Add terminal branch lengths for ASTRAL trees
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
a_p_trees <- lapply(1:length(a_p_trees), function(i){add.terminal.branches(a_p_trees[[i]], 1)})
a_f_trees <- lapply(1:length(a_f_trees), function(i){add.terminal.branches(a_f_trees[[i]], 1)})
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(a_f_trees) <- "multiPhylo"

## Order tips for plots
# Get the order for the tips (bottom species first, top species last)
primates_tip_order <- c("Mus musculus", "Tupaia chinensis", "Galeopterus variegatus",
                        "Otolemur garnettii", "Propithecus coquereli", "Microcebus murinus",
                        "Carlito syrichta",
                        "Callithrix jacchus","Aotus nancymaae", "Cebus capucinus imitator", "Saimiri boliviensis",
                        "Nomascus leucogenys", "Pongo abelii", "Gorilla gorilla",  "Homo sapiens", "Pan paniscus",  "Pan troglodytes",
                        "Colobus angolensis palliatus", "Piliocolobus tephrosceles", "Rhinopithecus bieti", "Rhinopithecus roxellana",
                        "Chlorocebus sabaeus", "Cercocebus atys", "Mandrillus leucophaeus", "Papio anubis", "Theropithecus gelada", "Macaca nemestrina", "Macaca fascicularis", "Macaca mulatta")

## Open ASTRAL unfiltered tree
unfiltered_astral_tree_file <- grep("ASTRAL", grep("NoTest", plot_tree_files, value = T), value = T)
unfiltered_astral_tree <- read.tree(unfiltered_astral_tree_file)
unfiltered_astral_tree <- root(unfiltered_astral_tree, roots_by_group[["Primates"]])
unfiltered_astral_tree$edge.length[which(is.nan(unfiltered_astral_tree$edge.length))] <- 1

## Create labels for densitree plots
# Create labels for the tips
primate_labels <- color.primates.by.clades(notest_concat_tree, color_palette = primate_colour_palette)
# Reorder primate_labels dataframe to make sure ggtree will be in desired order
primate_labels <- primate_labels[match(primates_tip_order, primate_labels$taxa),]
# Format taxa names
primate_labels$taxa <- gsub(" ", "_", primate_labels$taxa)

## Plot main figure ASTRAL pass
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% primate_labels +
  scale_y_reverse() +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.3) +
  scale_color_manual(values = primate_colour_palette) +
  labs(title = "Primates: ASTRAL Unfiltered") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(size = 16), legend.text = element_text (size = 12),
        legend.position = "left",
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 5))) +
  xlim(0,17)
astral_pass_plot <- ggdensitree(a_p_trees, tip.order = primate_labels$taxa, align.tips = TRUE, alpha = 0.8, color = "black") %<+% primate_labels +
  scale_y_reverse() +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.3) +
  scale_color_manual(values = primate_colour_palette, guide="none") +
  labs(title = "Primates: ASTRAL Pass") +
  xlim(-15.4, 3) +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
# Save plot
quilt <- unfiltered_astral_plot + astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Primates_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 14.5)
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", width = 14.5)
# Save long plot
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% primate_labels +
  scale_y_reverse() +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.3) +
  scale_color_manual(values = primate_colour_palette) +
  labs(title = "Primates: ASTRAL Unfiltered") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(size = 16), legend.text = element_text (size = 12),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 5))) +
  xlim(0,17)
quilt <- unfiltered_astral_plot / astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Primates_ggdensitree_long")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 12, width = 8)
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 12, width = 8)

#### Plot all densitrees
# Plot: ASTRAL, pass
a_p_densitree <- ggdensitree(a_p_trees, tip.order = primate_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% primate_labels +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4) +
  scale_color_manual(values = primate_colour_palette) +
  xlim(-15.4, 3) +
  labs(title = "Pass tests - ASTRAL trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(size = 18), legend.text = element_text (size = 16), legend.position = c(0.12,0.25),
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 6)))
# Plot: ASTRAL, fail
a_f_densitree <- ggdensitree(a_f_trees, tip.order = primate_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% primate_labels +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = FALSE, offset = 0.2, geom = "text", size = 4) +
  scale_color_manual(values = primate_colour_palette) +
  xlim(-15.4, 3) +
  labs(title = "Fail tests - ASTRAL trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
# Plot: CONCAT, pass
c_p_densitree <- ggdensitree(c_p_trees, tip.order = primate_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% primate_labels +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = FALSE, offset = 0.002, geom = "text", size = 4) +
  scale_color_manual(values = primate_colour_palette) +
  xlim(-0.189, 0.034) +
  labs(title = "Pass tests - Concatenated trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
# Plot: CONCAT, fail
c_f_densitree <- ggdensitree(c_f_trees, tip.order = primate_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% primate_labels +
  geom_tiplab(aes(label = short_lab, color = clade), parse = TRUE, show.legend = FALSE, offset = 0.002, geom = "text", size = 4) +
  scale_color_manual(values = primate_colour_palette) +
  xlim(-0.26, 0.038) +
  labs(title = "Fail tests - Concatenated trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
# Save plots
quilt <- (a_p_densitree | a_f_densitree) / (c_p_densitree | c_f_densitree) +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "Primates_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 12, width = 16, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 12, width = 16, units = "in")





#### Step 4: Tomatoes dataset ####
## Open Tomatoes dataset trees
# Get all files for species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract files for this dataset
plot_tree_files <- grep("Tomato", all_trees, value = TRUE)
# Extract trees for that combination of dataset and tree method
astral_trees_files <- grep("NoTest", grep("ASTRAL", plot_tree_files, value = TRUE), value = T, invert = T)
concat_trees_files <- grep("NoTest", grep("CONCAT", plot_tree_files, value = TRUE), value = T, invert = T)

# Separate into pass and fail trees
astral_pass_tree_files <- grep("pass", astral_trees_files, value = T)
astral_fail_tree_files <- grep("fail", astral_trees_files, value = T)
concat_pass_tree_files <- grep("pass", concat_trees_files, value = T)
concat_fail_tree_files <- grep("fail", concat_trees_files, value = T)
# Extract text file for each tree
astral_pass_trees_text <- unlist(lapply(astral_pass_tree_files, readLines))
astral_fail_trees_text <- unlist(lapply(astral_fail_tree_files, readLines))
concat_pass_trees_text <- unlist(lapply(concat_pass_tree_files, readLines))
concat_fail_trees_text <- unlist(lapply(concat_fail_tree_files, readLines))

## Open trees
# Read trees into a multiphylo object
a_p_trees <- read.tree(text = astral_pass_trees_text)
a_f_trees <- read.tree(text = astral_fail_trees_text)
c_p_trees <- read.tree(text = concat_pass_trees_text)
c_f_trees <- read.tree(text = concat_fail_trees_text)

## Root trees
a_p_trees <- lapply(a_p_trees, root, roots_by_group[["Tomatoes"]])
a_f_trees <- lapply(a_f_trees, root, roots_by_group[["Tomatoes"]])
c_p_trees <- lapply(c_p_trees, root, roots_by_group[["Tomatoes"]])
c_f_trees <- lapply(c_f_trees, root, roots_by_group[["Tomatoes"]])
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(a_f_trees) <- "multiPhylo"
class(c_p_trees) <- "multiPhylo"
class(c_f_trees) <- "multiPhylo"

## Add terminal branch lengths for ASTRAL trees
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
a_p_trees <- lapply(1:length(a_p_trees), function(i){add.terminal.branches(a_p_trees[[i]], 1)})
a_f_trees <- lapply(1:length(a_f_trees), function(i){add.terminal.branches(a_f_trees[[i]], 1)})
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(a_f_trees) <- "multiPhylo"

## Order tips for plots
# Get the order for the tips (bottom species first, top species last)
tomato_tip_order <- c("LA3909", "LA0436", "LA0429", "LA3124", "LA3475", "SL2.50", "LA1589", "LA1269",
                      "LA2933", "LA2133", "LA1322", "LA2172", "LA1316", "LA1028", "LA0444", "LA0107",
                      "LA1358", "LA2964", "LA2744", "LA4117", "LA1782", "LA1364", "LA1777", "LA0407",
                      "LA3778", "LA0716", "LA4116", "LA4126", "LA2951")
tomato_species_order <- rename.tomato.tips(tomato_tip_order)

## Create labels for densitree plots
tomato_labels <- color.code.tomato.clades(notest_astral_tree, taxa.numbers = FALSE, trimmed = FALSE, color = FALSE)
# Reorder tip_labels_df to make sure ggtree will be in desired order
tomato_labels <- tomato_labels[match(tomato_species_order, tomato_labels$taxa),]
# Replace taxa names with taxa numbers
tomato_labels$taxa <- tomato_tip_order

## Open ASTRAL unfiltered tree
unfiltered_astral_tree_file <- grep("ASTRAL", grep("NoTest", plot_tree_files, value = T), value = T)
unfiltered_astral_tree <- read.tree(unfiltered_astral_tree_file)
unfiltered_astral_tree <- root(unfiltered_astral_tree, roots_by_group[["Tomatoes"]])
unfiltered_astral_tree$edge.length[which(is.nan(unfiltered_astral_tree$edge.length))] <- 1

## Plot ASTRAL pass main figure
# Plot: ASTRAL, unfiltered
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 3.8) +
  scale_y_reverse() +
  scale_color_manual(values = tomato_colour_palette) +
  labs(title = "Tomatoes: ASTRAL Unfiltered") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text (size = 12),
        legend.position = "right",
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 6))) +
  xlim(0, 24)
# Plot: ASTRAL, pass
astral_pass_plot <- ggdensitree(a_p_trees, tip.order = tomato_labels$taxa, align.tips = TRUE, alpha = 0.8, color = "black") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 3.8) +
  scale_color_manual(values = tomato_colour_palette, guide = "none") +
  xlim(-16.6, 4.8) +
  labs(title = "Tomatoes: ASTRAL Pass") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
# Assemble the plot using patchwork
quilt <- unfiltered_astral_plot + astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
# Save the plot
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Tomatoes_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", width = 15)
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", width = 15)
# Save a long version too
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 3.8) +
  scale_y_reverse() +
  scale_color_manual(values = tomato_colour_palette) +
  labs(title = "Tomatoes: ASTRAL Unfiltered") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text (size = 12),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 6))) +
  xlim(0, 24)
quilt <- unfiltered_astral_plot / astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Tomatoes_ggdensitree_long")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 11, width = 9, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 11, width = 9, units = "in")

## Plot densitress
# Plot: ASTRAL, pass
a_p_densitree <- ggdensitree(a_p_trees, tip.order = tomato_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = TRUE, offset = 0.2, geom = "text", size = 4.5) +
  scale_color_manual(values = tomato_colour_palette) +
  xlim(-16.6, 4.8) +
  labs(title = "Pass tests - ASTRAL trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.title = element_text(size = 18), legend.text = element_text (size = 16), legend.position.inside = c(0.08,0.25),
        legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 6)))
# Plot: ASTRAL, fail
a_f_densitree <- ggdensitree(a_f_trees, tip.order = tomato_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = FALSE, offset = 0.2, geom = "text", size = 4.5) +
  scale_color_manual(values = tomato_colour_palette) +
  xlim(-19.58, 5.5) +
  labs(title = "Fail tests - ASTRAL trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))
# Plot: CONCAT, pass
c_p_densitree <- ggdensitree(c_p_trees, tip.order = tomato_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, offset = 0.0002, show.legend = FALSE, geom = "text", size = 4.5) +
  scale_color_manual(values = tomato_colour_palette) +
  xlim(-0.0287, 0.0078) +
  labs(title = "Pass tests - Concatenated trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))
# Plot: CONCAT, fail
c_f_densitree <- ggdensitree(c_f_trees, tip.order = tomato_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% tomato_labels +
  geom_tiplab(aes(label = lab, color = clade), parse = TRUE, show.legend = FALSE, offset = 0.0002, geom = "text", size = 4.5) +
  scale_color_manual(values = tomato_colour_palette) +
  xlim(-0.0272, 0.0078) +
  labs(title = "Fail tests - Concatenated trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))
# Assemble the plot using patchwork
quilt <- (a_p_densitree | a_f_densitree) / (c_p_densitree | c_f_densitree) +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
# Save the plot
densitree_name <- paste0(plot_dir, "Tomatoes_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 12, width = 18, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 12, width = 18, units = "in")



#### Step 5: Metazoans dataset ####
## Open Metazoan dataset trees
# Get all files for species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract files for this dataset
plot_tree_files <- grep("Metazoan", all_trees, value = TRUE)
# Extract trees for that combination of dataset and tree method
astral_trees_files <- grep("NoTest", grep("ASTRAL", plot_tree_files, value = TRUE), value = T, invert = T)
concat_trees_files <- grep("NoTest", grep("CONCAT", plot_tree_files, value = TRUE), value = T, invert = T)

# Separate into pass and fail trees
astral_pass_tree_files <- grep("pass", astral_trees_files, value = T)
concat_pass_tree_files <- grep("pass", concat_trees_files, value = T)
# Extract text file for each tree
astral_pass_trees_text <- unlist(lapply(astral_pass_tree_files, readLines))
concat_pass_trees_text <- unlist(lapply(concat_pass_tree_files, readLines))

## Open trees
# Read trees into a multiphylo object
a_p_trees <- read.tree(text = astral_pass_trees_text)
c_p_trees <- read.tree(text = concat_pass_trees_text)

## Root trees
a_p_trees <- lapply(a_p_trees, root, roots_by_group[["Metazoan"]])
c_p_trees <- lapply(c_p_trees, root, roots_by_group[["Metazoan"]])
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(c_p_trees) <- "multiPhylo"

## Add terminal branch lengths for ASTRAL trees
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
a_p_trees <- lapply(1:length(a_p_trees), function(i){add.terminal.branches(a_p_trees[[i]], 1)})
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"

## Order tips for plots
# Get the order for the tips (bottom species first, top species last)
metazoans_tip_order <- c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis",
                         "Euplokamis_dunlapae",  "Coeloplana_astericola", "Vallicula_sp", "Hormiphora_californica", "Hormiphora_palmata",
                         "Pleurobrachia_pileus", "Pleurobrachia_bachei",  "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA",
                         "Callianira_Antarctica", "Mertensiidae_sp_Antarctica", "Mertensiidae_sp_Washington_USA",  "Cydippida_sp", "Dryodora_glandiformis",
                         "Lobatolampea_tetragona",  "Beroe_abyssicola", "Beroe_sp_Antarctica", "Beroe_ovata",  "Beroe_sp_Queensland_Australia", "Beroe_forskalii",
                         "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA", "Bolinopsis_infundibulum", "Mnemiopsis_leidyi",
                         "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris", "Ctenophora_sp_Florida_USA",
                         "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum", "Aphrocallistes_vastus",
                         "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica", "Petrosia_ficiformis",
                         "Spongilla_lacustris", "Cliona_varians", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans",
                         "Kirkpatrickia_variolosa",
                         "Trichoplax_adhaerens",
                         "Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster", "Daphnia_pulex",
                         "Eunicella_verrucosa", "Acropora_digitifera", "Nematostella_vectensis", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata",
                         "Periphyla_periphyla", "Hydra_viridissima", "Hydra_vulgaris", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona", "Craseo_lathetica",
                         "Nanomia_bijuga", "Agalma_elegans")
metazoans_tip_order <- rev(metazoans_tip_order)

## Create labels for densitree plots
metazoan_labels <- color.code.metazoan.clades(notest_concat_tree, trimmed = "FALSE", color = FALSE)
# Reorder tip_labels_df to make sure ggtree will be in desired order
metazoan_labels <- metazoan_labels[match(metazoans_tip_order, metazoan_labels$taxa),]
# Format taxa names
metazoan_labels$taxa <- gsub(" ", "_", metazoan_labels$taxa)
metazoan_labels$short_lab_noformat <- shorten.short.names(metazoan_labels$short_lab_noformat)

## Open ASTRAL unfiltered tree
unfiltered_astral_tree_file <- grep("ASTRAL", grep("NoTest", plot_tree_files, value = T), value = T)
unfiltered_astral_tree <- read.tree(unfiltered_astral_tree_file)
unfiltered_astral_tree <- root(unfiltered_astral_tree, roots_by_group[["Metazoan"]])
unfiltered_astral_tree$edge.length[which(is.nan(unfiltered_astral_tree$edge.length))] <- 1

## Plot main figure
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% metazoan_labels +
  geom_tiplab(aes(label = short_lab_noformat, color = clade), parse = FALSE, show.legend = TRUE, offset = 0.2, geom = "text", size = 3, fontface = 3) +
  scale_y_reverse() +
  scale_color_manual(values = metazoan_colour_palette) +
  labs(title = "Metazoa: ASTRAL Unfiltered") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12), legend.text = element_text (size = 10),
        legend.position = c("left")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 4))) +
  xlim(0,14)
astral_pass_plot <- ggdensitree(a_p_trees, tip.order = metazoan_labels$taxa, align.tips = TRUE, alpha = 0.8, color = "black", layout = "slanted") %<+% metazoan_labels +
  geom_tiplab(aes(label = short_lab_noformat, color = clade), align = TRUE, parse = FALSE, show.legend = TRUE, geom = "text", size = 3, fontface = 3) +
  scale_color_manual(values = metazoan_colour_palette, guide = "none") +
  labs(title = "Metazoa: ASTRAL Pass") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  xlim(-12,3.8)
# Save the plot
quilt <- unfiltered_astral_plot + astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Metazoan_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 12, width = 14, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 12, width = 14, units = "in")
# Save long plot
quilt <- unfiltered_astral_plot / astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Metazoan_ggdensitree_long")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 20, width = 10, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 20, width = 10, units = "in")


## Plot densitress
# Plot: ASTRAL, pass
a_p_densitree <- ggdensitree(a_p_trees, tip.order = metazoan_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% metazoan_labels +
  geom_tiplab(aes(label = short_lab_noformat, color = clade), parse = FALSE, show.legend = TRUE, offset = 0.2, geom = "text", size = 3.5, fontface = 3) +
  scale_color_manual(values = metazoan_colour_palette) +
  xlim(-12, 3.2) +
  labs(title = "Pass tests - ASTRAL trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.title = element_text(size = 12), legend.text = element_text (size = 10), legend.position.inside = c(0.05,0.05),
        legend.key.size = unit(0.9, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 4)))
# Plot: CONCAT, pass
c_p_densitree <- ggdensitree(c_p_trees, tip.order = metazoan_labels$taxa, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% metazoan_labels +
  geom_tiplab(aes(label = short_lab_noformat, color = clade), parse = FALSE, offset = 0.0002, show.legend = FALSE, geom = "text", size = 3.5, fontface = 3) +
  scale_color_manual(values = metazoan_colour_palette) +
  xlim(-1.25, 0.35) +
  labs(title = "Pass tests - Concatenated trees") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))
# Assemble the plot using patchwork
quilt <- (a_p_densitree | c_p_densitree) +  plot_annotation(tag_levels = 'a', tag_suffix = ".") &  theme(plot.tag = element_text(size = 30))
# Save the plot
densitree_name <- paste0(plot_dir, "Metazoans_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 11, width = 16, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png")



#### Step 6: Plants dataset ####
## Open Plant dataset trees
# Get all files for species trees
all_trees <- paste0(maindir, "species_trees/", list.files(paste0(maindir, "species_trees/")))
# Extract files for this dataset
plot_tree_files <- grep("Plant", all_trees, value = TRUE)
# Extract trees for that combination of dataset and tree method
astral_trees_files <- grep("NoTest", grep("ASTRAL", plot_tree_files, value = T), value = T, invert = T)
concat_trees_files <- grep("NoTest", grep("CONCAT", plot_tree_files, value = T), value = T, invert = T)

# Separate into pass and fail trees
astral_pass_tree_files <- grep("pass", astral_trees_files, value = T)
concat_pass_tree_files <- grep("pass", concat_trees_files, value = T)
# Extract text file for each tree
astral_pass_trees_text <- unlist(lapply(astral_pass_tree_files, readLines))
concat_pass_trees_text <- unlist(lapply(concat_pass_tree_files, readLines))

## Open trees
# Read pass test trees into a multiphylo object
a_p_trees <- read.tree(text = astral_pass_trees_text)
c_p_trees <- read.tree(text = concat_pass_trees_text)

## Create labels for densitree plots
# Color code clades
plant_labs <- color.plants.by.clades(tree = a_p_trees[[1]], color_palette = plants_color_palette, clade_df = annotation_df)

## Root trees
# Extract Chromista taxa: these will be the outgroups
chromista_tips <- plant_labs$Code[grep("Chromista", plant_labs$Very.Brief.Classification)]
a_p_trees <- lapply(a_p_trees, root, chromista_tips)
c_p_trees <- lapply(c_p_trees, root, chromista_tips)
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"
class(c_p_trees) <- "multiPhylo"

## Add terminal branch lengths for ASTRAL trees
# If tree estimation method is ASTRAL, add an arbitrary terminal branch length
a_p_trees <- lapply(1:length(a_p_trees), function(i){add.terminal.branches(a_p_trees[[i]], 1)})
# Convert object class from "list" into "multiPhylo
class(a_p_trees) <- "multiPhylo"

## Open ASTRAL unfiltered tree
unfiltered_astral_tree_file <- grep("ASTRAL", grep("NoTest", plot_tree_files, value = T), value = T)
unfiltered_astral_tree <- read.tree(unfiltered_astral_tree_file)
unfiltered_astral_tree <- root(unfiltered_astral_tree, chromista_tips)
unfiltered_astral_tree$edge.length[which(is.nan(unfiltered_astral_tree$edge.length))] <- 1

## Plot main figure
unfiltered_astral_plot <- ggtree(unfiltered_astral_tree, alpha = 1, color = "black") %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 2, alpha = 0.75) +
  scale_color_manual(values = plants_color_palette) +
  scale_y_reverse() +
  labs(title = "Pass tests - ASTRAL trees\n") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text (size = 14),
        legend.position.inside = c(0.15,0.23)) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 4)))
astral_pass_plot <- ggdensitree(a_p_trees, align.tips = TRUE, alpha = 0.8, color = "black") %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 2, alpha = 0.75) +
  scale_color_manual(values = plants_color_palette, guide = "none") +
  scale_y_reverse() +
  labs(title = "Pass tests - ASTRAL trees\n") +
  theme(axis.ticks.x = element_line(color = "white"),
        axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# Save the plot
quilt <- unfiltered_astral_plot + astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Plants_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 10, width = 10, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 10, width = 10, units = "in")
# Save long plot
quilt <- unfiltered_astral_plot / astral_pass_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ".") &
  theme(plot.tag = element_text(size = 30))
densitree_name <- paste0(plot_dir, "mainfig_ASTRAL_Plants_ggdensitree_long")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf", height = 20, width = 10, units = "in")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png", height = 20, width = 10, units = "in")

## Plot densitress
# Plot: ASTRAL, pass
a_p_densitree <- ggdensitree(a_p_trees, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 2, alpha = 0.75) +
  scale_color_manual(values = plants_color_palette) +
  xlim(-39.2, 0) +
  scale_y_reverse() +
  labs(title = "Pass tests - ASTRAL trees\n") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.title = element_text(size = 18), legend.text = element_text (size = 16), legend.position.inside = c(0.15,0.3),
        legend.key.size = unit(0.9, "lines")) +
  guides(color = guide_legend(title = "Clade legend", override.aes=list(label = "Sp.", size = 4)))
# Plot: CONCAT, pass
c_p_densitree <- ggdensitree(c_p_trees, align.tips = TRUE, alpha = 0.5, color = "steelblue") %<+% plant_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 2, alpha = 0.75, show.legend = FALSE) +
  scale_color_manual(values = plants_color_palette) +
  xlim(-3.18, 0) +
  scale_y_reverse() +
  labs(title = "Pass tests - Concatenated trees\n") +
  theme(axis.ticks.x = element_line(color = "white"), axis.line.x = element_line(color = "white"),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

## Assemble the plot using patchwork
quilt <- (a_p_densitree | c_p_densitree) +  plot_annotation(tag_levels = 'a', tag_suffix = ".") &  theme(plot.tag = element_text(size = 30))

## Save the plot
densitree_name <- paste0(plot_dir, "Plants_ggdensitree")
ggsave(filename = paste0(densitree_name, ".pdf"), plot = quilt, device = "pdf")
ggsave(filename = paste0(densitree_name, ".png"), plot = quilt, device = "png")



