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
# annotations_csv_file    <- location of misc/annotations.csv file from the Leebens-Mack (2019) "Data from 1000 Plants Transcriptomes" data repository - used to assign taxa names and clades
# roots_by_groups         <- set which taxa is outgroup for each dataset: Primates, Tomatoes, Metazoans, and Plants

### Caitlin's paths ###
# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
annotation_csv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/misc/annotations.csv"
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
old.par <- par(mar = c(0, 0, 0, 0))
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}
# Plot all plots in two figures
for (i in 1:2){
  if (i == 1){
    primate_ordered_trees <- c(grep("geneconv_pass", primate_astral_trees, value = T), 
                              grep("geneconv_fail", primate_astral_trees, value = T),
                              grep("maxchi_pass", primate_astral_trees, value = T),
                              grep("maxchi_fail", primate_astral_trees, value = T))
    plot_titles <- c("GENECONV, pass", "GENECONV, fail", "MaxChi, Pass", "MaxChi, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(20, 16, 17, 16)
  } else if (i == 2){
    primate_ordered_trees <- c(grep("PHI_pass", primate_astral_trees, value = T), 
                              grep("PHI_fail", primate_astral_trees, value = T),
                              grep("allTests_pass", primate_astral_trees, value = T),
                              grep("allTests_fail", primate_astral_trees, value = T))
    plot_titles <- c("PHI, pass", "PHI, fail", "All tests, Pass", "All tests, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(20, 9, 18, 16)
  }
  for (j in 1:length(primate_ordered_trees)){
    # Identify tree path
    f_tree_path <- primate_ordered_trees[j]
    # Read tree
    f_tree <- read.tree(f_tree_path)
    # Root tree
    f_tree <- root(f_tree, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
    # Change terminal edge.length to 0.5
    f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
    # Remove underscores from tips
    f_tree$tip.label <- gsub("_", " ", f_tree$tip.label)
    # Color code clades
    f_colors <- primate_colour_palette
    f_labs <-color.primates.by.clades(f_tree, color_palette = f_colors)
    # Add to tree list
    tree_list[[j]] <- f_tree
  }
  # Plot the four panels
  plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[1]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5))) +
    xlim(0, xlimits[1])
  plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[2]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[2])
  plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[3]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[3])
  plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[4]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[4])
  quilt1 <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2, nrow = 2) + 
    plot_annotation(title = "Primates - ASTRAL", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
  quilt1_path <- paste0(plot_dir, "Primates_ASTRAL_panels_", i, ".pdf")
  ggsave(filename = quilt1_path, plot = quilt1, width = 14, height = 10)
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.treefile", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}
# Plot all plots in two figures
for (i in 1:2){
  if (i == 1){
    primate_ordered_trees <- c(grep("geneconv_pass", primate_concat_trees, value = T), 
                               grep("geneconv_fail", primate_concat_trees, value = T),
                               grep("maxchi_pass", primate_concat_trees, value = T),
                               grep("maxchi_fail", primate_concat_trees, value = T))
    plot_titles <- c("GENECONV, pass", "GENECONV, fail", "MaxChi, Pass", "MaxChi, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(0.22, 0.22, 0.21, 0.22)
  } else if (i == 2){
    primate_ordered_trees <- c(grep("PHI_pass", primate_concat_trees, value = T), 
                               grep("PHI_fail", primate_concat_trees, value = T),
                               grep("allTests_pass", primate_concat_trees, value = T),
                               grep("allTests_fail", primate_concat_trees, value = T))
    plot_titles <- c("PHI, pass", "PHI, fail", "All tests, Pass", "All tests, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(0.22, 0.28, 0.21, 0.22)
  }
  for (j in 1:length(primate_ordered_trees)){
    # Identify tree path
    f_tree_path <- primate_ordered_trees[j]
    # Read tree
    f_tree <- read.tree(f_tree_path)
    # Root tree
    f_tree <- root(f_tree, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
    # Change terminal edge.length to 0.5
    f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
    # Remove underscores from tips
    f_tree$tip.label <- gsub("_", " ", f_tree$tip.label)
    # Color code clades
    f_colors <- primate_colour_palette
    f_labs <-color.primates.by.clades(f_tree, color_palette = f_colors)
    # Add to tree list
    tree_list[[j]] <- f_tree
  }
  # Plot the four panels
  plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[1]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5))) +
    xlim(0, xlimits[1])
  plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[2]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[2])
  plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[3]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[3])
  plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[4]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[4])
  quilt1 <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2, nrow = 2) + 
    plot_annotation(title = "Primates - CONCAT", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
  quilt1_path <- paste0(plot_dir, "Primates_CONCAT_panels_", i, ".pdf")
  ggsave(filename = quilt1_path, plot = quilt1, width = 13.5, height = 10)
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}
# Plot all plots in two figures
for (i in 1:2){
  if (i == 1){
    tomato_ordered_trees <- c(grep("geneconv_pass", tomato_astral_trees, value = T), 
                              grep("geneconv_fail", tomato_astral_trees, value = T),
                              grep("maxchi_pass", tomato_astral_trees, value = T),
                              grep("maxchi_fail", tomato_astral_trees, value = T))
    plot_titles <- c("GENECONV, pass", "GENECONV, fail", "MaxChi, Pass", "MaxChi, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(20, 22, 19, 23)
  } else if (i == 2){
    tomato_ordered_trees <- c(grep("PHI_pass", tomato_astral_trees, value = T), 
                              grep("PHI_fail", tomato_astral_trees, value = T),
                              grep("allTests_pass", tomato_astral_trees, value = T),
                              grep("allTests_fail", tomato_astral_trees, value = T))
    plot_titles <- c("PHI, pass", "PHI, fail", "All tests, Pass", "All tests, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(20, 22, 18, 21)
  }
  for (j in 1:length(tomato_ordered_trees)){
    f_tree_path <- tomato_ordered_trees[j]
    # Read tree
    f_tree <- read.tree(f_tree_path)
    # Root tree by outgroup
    f_tree <- root(f_tree, outgroup = roots_by_group[["Tomatoes"]], resolve.root = TRUE)
    # Change terminal edge.length
    f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.5
    # Rename tip labels to have scientific names (not just numbers)
    f_tree$tip.label <- rename.tomato.tips(f_tree$tip.label)
    # Color code clades
    f_colors <- tomato_colour_palette
    f_labs <-color.code.tomato.clades(f_tree, taxa.numbers = FALSE, trimmed = FALSE)
    # Add to tree list
    tree_list[[j]] <- f_tree
  }
  # Plot the four panels
  plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[1]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))  +
    xlim(0, xlimits[1])
  plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[2]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[2])
  plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[3]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[3])
  plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 3) + 
    scale_y_reverse() +  
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[4]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[4])
  quilt1 <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2, nrow = 2) + 
    plot_annotation(title = "Tomatoes - ASTRAL", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
  quilt1_path <- paste0(plot_dir, "Tomatoes_ASTRAL_panels_", i, ".pdf")
  ggsave(filename = quilt1_path, plot = quilt1, width = 16, height = 10)
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.")))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.treefile", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf")
}
# Plot all plots in two figures
for (i in 1:2){
  if (i == 1){
    tomato_ordered_trees <- c(grep("geneconv_pass", tomato_concat_trees, value = T), 
                              grep("geneconv_fail", tomato_concat_trees, value = T),
                              grep("maxchi_pass", tomato_concat_trees, value = T),
                              grep("maxchi_fail", tomato_concat_trees, value = T))
    plot_titles <- c("GENECONV, pass", "GENECONV, fail", "MaxChi, Pass", "MaxChi, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(0.034, 0.032, 0.032, 0.032)
  } else if (i == 2){
    tomato_ordered_trees <- c(grep("PHI_pass", tomato_concat_trees, value = T), 
                              grep("PHI_fail", tomato_concat_trees, value = T),
                              grep("allTests_pass", tomato_concat_trees, value = T),
                              grep("allTests_fail", tomato_concat_trees, value = T))
    plot_titles <- c("PHI, pass", "PHI, fail", "All tests, Pass", "All tests, Fail")
    tree_list <- list("One" = NA,
                      "Two" = NA,
                      "Three" = NA,
                      "Four" = NA)
    xlimits <- c(0.034, 0.032, 0.032, 0.033)
  }
  for (j in 1:length(tomato_ordered_trees)){
    f_tree_path <- tomato_ordered_trees[j]
    # Read tree
    f_tree <- read.tree(f_tree_path)
    # Root tree by outgroup
    f_tree <- root(f_tree, outgroup = roots_by_group[["Tomatoes"]], resolve.root = TRUE)
    # Rename tip labels to have scientific names (not just numbers)
    f_tree$tip.label <- rename.tomato.tips(f_tree$tip.label)
    # Color code clades
    f_colors <- tomato_colour_palette
    f_labs <-color.code.tomato.clades(f_tree, taxa.numbers = FALSE, trimmed = FALSE)
    # Add to tree list
    tree_list[[j]] <- f_tree
  }
  # Plot the four panels
  plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
    geom_tiplab(aes(label=lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors) +
    labs(title = plot_titles[1]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
          legend.title = element_text(size = 15), 
          legend.text = element_text (size = 12), 
          legend.position	= "left",) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5))) +
    xlim(0, xlimits[1])
  plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[2]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[2])
  plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    coord_cartesian(clip = 'off') + 
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[3]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[3])
  plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
    geom_tiplab(aes(label = lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
    scale_y_reverse() +  
    scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
    labs(title = plot_titles[4]) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5)) +
    xlim(0, xlimits[4])
  quilt1 <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2, nrow = 2) + 
    plot_annotation(title = "Tomatoes - CONCAT", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
  quilt1_path <- paste0(plot_dir, "Tomatoes_CONCAT_panels_", i, ".pdf")
  ggsave(filename = quilt1_path, plot = quilt1, width = 12, height = 10)
}



#### Step 7: Plot all ASTRAL Metazoan trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_astral_trees <- grep("ASTRAL", metazoan_tree_files, value = TRUE)
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 9)
}
# Plot all plots in two figures
metazoan_ordered_trees <- c(grep("NoTest", metazoan_astral_trees, value = T), 
                            grep("geneconv_pass", metazoan_astral_trees, value = T),
                            grep("maxchi_pass", metazoan_astral_trees, value = T),
                            grep("PHI_pass", metazoan_astral_trees, value = T))
plot_titles <- c("No Tests\n(Unfiltered)", "GENECONV, pass", "MaxChi, Pass", "PHI, Pass")
tree_list <- list("NoTest" = NA,
                  "geneconv_pass" = NA,
                  "MaxChi_pass" = NA,
                  "PHI_pass" = NA)
for (i in 1:length(metazoan_ordered_trees)){
  # Select tree
  f_tree_path <- metazoan_ordered_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Metazoan"]], resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 0.4
  f_tree$edge.length <- f_tree$edge.length/2
  # Color code clades
  f_colors <- metazoan_colour_palette
  f_labs <-color.code.metazoan.clades(f_tree, trimmed = "FALSE")
  # Add tree to list
  tree_list[[i]] <- f_tree
}
plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[1]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "bottom",) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[2]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5))
plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[3]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "bottom",) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[4]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5))
quilt1 <- plot1 + plot2  + plot_layout(ncol = 1, nrow = 2) + 
  plot_annotation(title = "Metazoa - ASTRAL", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt1_path <- paste0(plot_dir, "Metazoa_ASTRAL_panels_1.pdf")
ggsave(filename = quilt1_path, plot = quilt1, height = 13, width = 12)
quilt2 <- plot3 + plot4 + plot_layout(ncol = 1, nrow = 2) + 
  plot_annotation(title = "Metazoa - ASTRAL", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt2_path <- paste0(plot_dir, "Metazoa_ASTRAL_panels_2.pdf")
ggsave(filename = quilt2_path, plot = quilt2, height = 13, width = 12)



#### Step 8: Plot all CONCAT Metazoan trees ####
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
metazoan_tree_files <- grep("Metazoan", tree_files, value = TRUE)
metazoan_concat_trees <- grep("CONCAT", metazoan_tree_files, value = TRUE)
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
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.treefile", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 9)
}
# Plot all plots in two figures
metazoan_ordered_trees <- c(grep("NoTest", metazoan_concat_trees, value = T), 
                            grep("geneconv_pass", metazoan_concat_trees, value = T),
                            grep("maxchi_pass", metazoan_concat_trees, value = T),
                            grep("PHI_pass", metazoan_concat_trees, value = T))
plot_titles <- c("No Tests\n(Unfiltered)", "GENECONV, pass", "MaxChi, Pass", "PHI, Pass")
tree_list <- list("NoTest" = NA,
                  "geneconv_pass" = NA,
                  "MaxChi_pass" = NA,
                  "PHI_pass" = NA)
for (i in 1:length(metazoan_ordered_trees)){
  # Select tree
  f_tree_path <- metazoan_ordered_trees[i]
  # Read tree
  f_tree <- read.tree(f_tree_path)
  # Root tree by outgroup
  f_tree <- root(f_tree, outgroup = roots_by_group[["Metazoan"]], resolve.root = TRUE)
  # Change terminal edge.length
  f_tree$edge.length[which(is.nan(f_tree$edge.length))] <- 1
  # Color code clades
  f_colors <- metazoan_colour_palette
  f_labs <-color.code.metazoan.clades(f_tree, trimmed = "FALSE")
  # Add tree to list
  tree_list[[i]] <- f_tree
}
plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[1]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "bottom",) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[2]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5))
plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[3]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "bottom",) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot4 <- ggtree(tree_list[[4]]) %<+% f_labs +
  geom_tiplab(aes(label = long_lab, color = clade), parse=T, show.legend = TRUE, offset = 0, geom = "text", size = 2.5) + 
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[4]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 0.5))
quilt1 <- plot1 + plot2  + plot_layout(ncol = 1, nrow = 2) + 
  plot_annotation(title = "Metazoa - CONCAT", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt1_path <- paste0(plot_dir, "Metazoa_CONCAT_panels_1.pdf")
ggsave(filename = quilt1_path, plot = quilt1, height = 15, width = 13)
quilt2 <- plot3 + plot4 + plot_layout(ncol = 1, nrow = 2) + 
  plot_annotation(title = "Metazoa - CONCAT", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt2_path <- paste0(plot_dir, "Metazoa_CONCAT_panels_2.pdf")
ggsave(filename = quilt2_path, plot = quilt2, height = 15, width = 13)



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
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.tre", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 10)
}

# Plot all plots in one figure
plant_ordered_trees <- c(grep("NoTest", plants_astral_trees, value = T), 
                         grep("maxchi_pass", plants_astral_trees, value = T),
                         grep("PHI_pass", plants_astral_trees, value = T))
plot_titles <- c("No Tests\n(Unfiltered)", "MaxChi, Pass", "PHI, Pass")
tree_list <- list("NoTest" = NA,
                  "MaxChi_pass" = NA,
                  "PHI_pass" = NA)
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
  # Add tree to list
  tree_list[[i]] <- f_tree
}
plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1.5, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[1]) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "left",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1.5, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[2]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5))
plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1.5, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[3]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5))
quilt <- plot1 + plot2 + plot3 + plot_layout(ncol = 3, nrow = 1) + 
  plot_annotation(title = "Plants - ASTRAL", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt_path <- paste0(plot_dir, "Plants_ASTRAL_all_trees.pdf")
ggsave(filename = quilt_path, plot = quilt, width = 10, height = 8)



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
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.subtitle = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 0.5) ) +
    guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
  # Create plot name
  f_name <- paste0(plot_dir, gsub("\\.raxml\\.bestTree", "", basename(f_tree_path)))
  # Save plot
  ggsave(filename = paste0(f_name, ".pdf"), plot = f_plot, device = "pdf", height = 10)
}

# Plot all plots in one figure
plant_ordered_trees <- c(grep("NoTest", plants_concat_trees, value = T), 
                         grep("maxchi_pass", plants_concat_trees, value = T),
                         grep("PHI_pass", plants_concat_trees, value = T))
plot_titles <- c("No Tests\n(Unfiltered)", "MaxChi, Pass", "PHI, Pass")
tree_list <- list("NoTest" = NA,
                  "MaxChi_pass" = NA,
                  "PHI_pass" = NA)
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
  # Add tree to list
  tree_list[[i]] <- f_tree
}
plot1 <- ggtree(tree_list[[1]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors) +
  labs(title = plot_titles[1]) +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 15), 
        legend.text = element_text (size = 12), 
        legend.position	= "left",
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(label = "Sp.", size = 4.5)))
plot2 <- ggtree(tree_list[[2]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[2]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 0.5))
plot3 <- ggtree(tree_list[[3]]) %<+% f_labs +
  geom_tippoint(aes(color = Very.Brief.Classification), size = 1, alpha = 0.7) +
  scale_y_reverse() +  
  coord_cartesian(clip = 'off') + 
  scale_color_manual(name = "Clade legend", values = f_colors, guide = "none") +
  labs(title = plot_titles[3]) +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5, vjust = 0.5))
quilt <- plot1 + plot2 + plot3 + plot_layout(ncol = 3, nrow = 1) + 
  plot_annotation(title = "Plants - CONCAT", theme = theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold")))
quilt_path <- paste0(plot_dir, "Plants_CONCAT_all_trees.pdf")
ggsave(filename = quilt_path, plot = quilt, width = 10, height = 8)



#### Step 8: Plot differences in unfiltered Primate trees ####
# Open trees
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
notest_files <- grep("NoTest", tree_files, value = T)
p_concat <- read.tree(grep("Primates", grep("CONCAT", notest_files, value = T), value = T))
p_astral <- read.tree(grep("Primates", grep("ASTRAL", notest_files, value = T), value = T))
p_concat <- root(p_concat, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
p_astral <- root(p_astral, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
# Extract clade
concat_clade <- keep.tip(p_concat, c("Aotus_nancymaae", "Callithrix_jacchus", "Cebus_capucinus_imitator", "Saimiri_boliviensis"))
astral_clade <- keep.tip(p_astral, c("Aotus_nancymaae", "Callithrix_jacchus", "Cebus_capucinus_imitator", "Saimiri_boliviensis"))
concat_clade$tip.label <- gsub("_", " ", concat_clade$tip.label)
astral_clade$tip.label <- gsub("_", " ", astral_clade$tip.label)
concat_clade <- ladderize(concat_clade)
astral_clade <- ladderize(astral_clade)
astral_clade$edge.length[which(is.nan(astral_clade$edge.length))] <- 1
# Plot clade
concat_plot <- ggtree(concat_clade) +
  geom_rootedge(0.005) +
  geom_tiplab(color = "black", size = 5, fontface = "italic") +
  scale_y_reverse() +
  xlim(-0.005, 0.03) +
  labs(title = "CONCAT") +
  theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5, face = "bold", margin = margin(b = 20)))
astral_plot <- ggtree(astral_clade) +
  geom_rootedge(0.2) +
  geom_tiplab(color = "black", size = 5, fontface = "italic") +
  scale_y_reverse() +
  xlim(-0.22, 2.5) +
  labs(title = "ASTRAL") +
  theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5, face = "bold", margin = margin(b = 20)))
# Save plot
quilt <- astral_plot + concat_plot + plot_layout(ncol = 2, nrow = 1) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
quilt_path <- paste0(plot_dir, "Primates_NoTest_tree_differences.pdf")
ggsave(filename = quilt_path, plot = quilt, width = 13, height = 6)


