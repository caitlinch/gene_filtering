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
library(ggtree) # for plotting phylogenetic trees and densitress (ggdensitree)
library(ggtext) # for nice tree plots
library(patchwork) # for collating plots
library(TreeTools) # for CollapseNode function

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



#### Step 3: Primates dataset ####
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



#### Step 4: Tomatoes dataset ####
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



#### Step 5: Metazoans dataset ####
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



#### Step 5: Plants dataset ####
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