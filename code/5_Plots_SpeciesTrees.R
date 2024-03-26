### gene_filtering/code/5_Plots_SpeciesTrees.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2024

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

# Laptop
maindir <- "~/Repositories/gene_filtering/"
plot_dir <- "/Users/caitlin/Documents/PhD/Ch01_EmpiricalTreelikeness/plots/"
annotation_csv_file <- "/Users/caitlin/Documents/PhD/Ch01_EmpiricalTreelikeness/annotations.csv"

### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ape) # functions: read.tree, Ntip, root
library(phangorn)
library(ggplot2) # for nice plots
library(ggtree) # for plotting phylogenetic trees and densitress (ggdensitree)
library(ggtext) # for nice tree plots
library(patchwork) # for collating plots
library(TreeTools) # for CollapseNode function
library(colorBlindness) # For Plants clade labels

# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "species_trees/")

# Add roots for each dataset
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))

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



#### Step 3: Plot all ASTRAL Primates trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
primate_tree_files <- grep("Primates", tree_files, value = TRUE)
primate_astral_trees <- grep("ASTRAL", primate_tree_files, value = TRUE)
plot_titles <- paste0("Primates - ", c("All tests, Fail", "All tests, Pass", 
                                       "GENECONV, Fail", "GENECONV, Pass",
                                       "MaxChi, Fail", "MaxChi, Pass", 
                                       "No Tests (Unfiltered)", 
                                       "PHI, Fail", "PHI, Pass"))
for (i in 1:length(primate_astral_trees)){
  # Select tree
  f_tree_path <- primate_astral_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
  # Change terminal edge.length to 0.5
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
  # Remove underscores from tips
  f_tree$tip.label <- gsub("_", " ", f_tree$tip.label)
  # Color code clades
  f_colors <- primate_colour_palette
  f_labs <-color.primates.by.clades(f_tree, color_palette = f_colors)
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0.002, geom = "text", size = 5) + 
    scale_y_reverse() +  
    scale_x_continuous(breaks = seq(0,15,3)) +
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(6, 180, 6, 6)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "ASTRAL tree") +
    theme(axis.text.x = element_text(size = 12),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}



#### Step 4: Plot all CONCAT Primates trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
primate_tree_files <- grep("Primates", tree_files, value = TRUE)
primate_concat_trees <- grep("CONCAT", primate_tree_files, value = TRUE)
plot_titles <- paste0("Primates - ", c("All tests, Fail", "All tests, Pass", 
                                       "GENECONV, Fail", "GENECONV, Pass",
                                       "MaxChi, Fail", "MaxChi, Pass", 
                                       "No Tests (Unfiltered)", 
                                       "PHI, Fail", "PHI, Pass"))
for (i in 1:length(primate_concat_trees)){
  # Select tree
  f_tree_path <- primate_concat_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
  # Change terminal edge.length to 0.5
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
  # Remove underscores from tips
  f_tree$tip.label <- gsub("_", " ", f_tree$tip.label)
  # Color code clades
  f_colors <- primate_colour_palette
  f_labs <-color.primates.by.clades(f_tree, color_palette = f_colors)
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0.002, geom = "text", size = 5) + 
    scale_x_continuous(breaks = seq(0,0.25,0.05)) +
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(6, 120, 6, 6)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.contree", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}



#### Step 5: Plot all ASTRAL Tomatoes trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
tomato_tree_files <- grep("Tomatoes", tree_files, value = TRUE)
tomato_astral_trees <- grep("ASTRAL", tomato_tree_files, value = TRUE)
plot_titles <- paste0("Tomatoes - ", c("All tests, Fail", "All tests, Pass", 
                                       "GENECONV, Fail", "GENECONV, Pass",
                                       "MaxChi, Fail", "MaxChi, Pass", 
                                       "No Tests (Unfiltered)", 
                                       "PHI, Fail", "PHI, Pass"))
for (i in 1:length(tomato_astral_trees)){
  # Select tree
  f_tree_path <- tomato_astral_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Tomatoes"]], resolve.root = TRUE)
  # Rename tip labels to have scientific names (not just numbers)
  f_tree$tip.label <- rename.tomato.tips(f_tree$tip.label)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 1.5
  # Color code clades
  f_colors <- tomato_colour_palette
  f_labs <-color.code.tomato.clades(f_tree, taxa.numbers = FALSE, trimmed = FALSE)
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(6, 180, 6, 6)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "ASTRAL tree") +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}


#### Step 6: Plot all CONCAT Tomatoes trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
tomato_tree_files <- grep("Tomatoes", tree_files, value = TRUE)
tomato_concat_trees <- grep("CONCAT", tomato_tree_files, value = TRUE)
plot_titles <- paste0("Tomatoes - ", c("All tests, Fail", "All tests, Pass", 
                                       "GENECONV, Fail", "GENECONV, Pass",
                                       "MaxChi, Fail", "MaxChi, Pass", 
                                       "No Tests (Unfiltered)", 
                                       "PHI, Fail", "PHI, Pass"))
for (i in 1:length(tomato_concat_trees)){
  # Select tree
  f_tree_path <- tomato_concat_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Tomatoes"]], resolve.root = TRUE)
  # Rename tip labels to have scientific names (not just numbers)
  f_tree$tip.label <- rename.tomato.tips(f_tree$tip.label)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 1.5
  # Color code clades
  f_colors <- tomato_colour_palette
  f_labs <-color.code.tomato.clades(f_tree, taxa.numbers = FALSE, trimmed = FALSE)
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(6, 180, 6, 6)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.contree", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}



#### Step 7: Plot all ASTRAL Metazoan trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_astral_trees <- grep("ASTRAL", metazoa_tree_files, value = TRUE)
plot_titles <- paste0("Metazoa - ", c("GENECONV, Pass", 
                                      "MaxChi, Pass", 
                                      "No Tests (Unfiltered)", 
                                      "PHI, Pass"))
for (i in 1:length(metazoan_astral_trees)){
  # Select tree
  f_tree_path <- metazoan_astral_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Metazoan"]], resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 1
  # Color code clades
  f_colors <- metazoan_colour_palette
  f_labs <-color.code.metazoan.clades(f_tree, trimmed = "FALSE")
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(0, 10, 0, 0)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 9)
}



#### Step 8: Plot all CONCAT Metazoan trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_concat_trees <- grep("CONCAT", metazoa_tree_files, value = TRUE)
plot_titles <- paste0("Metazoa - ", c("GENECONV, Pass", 
                                      "MaxChi, Pass", 
                                      "No Tests (Unfiltered)", 
                                      "PHI, Pass"))
for (i in 1:length(metazoan_concat_trees)){
  # Select tree
  f_tree_path <- metazoan_concat_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Metazoan"]], resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 1
  # Color code clades
  f_colors <- metazoan_colour_palette
  f_labs <-color.code.metazoan.clades(f_tree, trimmed = "FALSE")
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(0, 10, 0, 0)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.contree", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 9)
}



#### Step 7: Plot all ASTRAL Plants trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
plants_tree_files <- grep("Plants", tree_files, value = TRUE)
plants_astral_trees <- grep("ASTRAL", plants_tree_files, value = TRUE)
plot_titles <- paste0("Plants - ", c("MaxChi, Pass", 
                                      "No Tests (Unfiltered)", 
                                      "PHI, Pass"))
annotation_df <- read.csv(annotation_csv_file)
for (i in 1:length(plants_astral_trees)){
  # Select tree
  f_tree_path <- plants_astral_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Color code clades
  f_colors <- plants_color_palette
  f_labs <- color.plants.by.clades(tree = f_tree, color_palette = plants_color_palette, clade_df = annotation_df)
  # Root tree by outgroup
  f_root <- f_labs$Code[which(f_labs$Very.Brief.Classification == "Chromista")]
  f_tree <- root(f_tree, outgroup = f_root, resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tippoint(aes(color = Very.Brief.Classification), size = 3, alpha = 1) +
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(0, 5, 0, 0)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 10)
}



#### Step 7: Plot all CONCAT Plants trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
plants_tree_files <- grep("Plants", tree_files, value = TRUE)
plants_concat_trees <- grep("CONCAT", plants_tree_files, value = TRUE)
plot_titles <- paste0("Plants - ", c("MaxChi, Pass", 
                                     "No Tests (Unfiltered)", 
                                     "PHI, Pass"))
annotation_df <- read.csv(annotation_csv_file)
for (i in 1:length(plants_concat_trees)){
  # Select tree
  f_tree_path <- plants_concat_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Color code clades
  f_colors <- plants_color_palette
  f_labs <- color.plants.by.clades(tree = f_tree, color_palette = plants_color_palette, clade_df = annotation_df)
  # Root tree by outgroup
  f_root <- f_labs$Code[which(f_labs$Very.Brief.Classification == "Chromista")]
  f_tree <- root(f_tree, outgroup = f_root, resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
  # Create plot
  f_plot <- ggtree(f_tree) %<+% f_labs +
    geom_tippoint(aes(color = Very.Brief.Classification), size = 3, alpha = 1) +
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    theme_tree2(plot.margin=margin(0, 5, 0, 0)) + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[i], subtitle = "CONCAT tree") +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",
          legend.position.inside = c(0.1, 0.2),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 10)
}
