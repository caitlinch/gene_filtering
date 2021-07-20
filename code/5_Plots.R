### empirical_treelikeness/4_Plots.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Caitlin Cherryh 2021


##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
print("Set filepaths")
# maindir             <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# plots_dir           <- for saving plots and analyses. This file should contain a folder for each dataset (where the folder name and corresponding dataset name are identical)
# species_tree_folder <- folder containing the species trees estimated in ASTRAL and IQ-Tree
# treelikeness_file   <- file containing results of recombination detection tests
# datasets            <- set name(s) for the dataset(s)

# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"
species_tree_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/ZZ_trimmed_loci_old/"
treelikeness_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/02_1KP_Pease2016_Strassert2021_Vanderpool2020_collated_RecombinationDetection_TrimmedLoci.csv"

# Datasets
datasets <- c("Vanderpool2020", "Pease2016", "Strassert2021", "1KP")
datasets = c("Vanderpool2020", "Pease2016")
dataset = "Pease2016"
roots <- c("Pease2016" = "LA4116", "Vanderpool2020" = "Mus_musculus")
n_tips <- c("Pease2016" = 29, "Vanderpool2020" = 29)
taxa_order <- list("Pease2016" = c("LA4116", "LA4126", "LA2951", "LA3778", "LA0716", "LA1777", "LA0407", "LA4117", "LA1782", "LA2964", "LA0444", "LA0107", "LA1358", "LA2744",
                                   "LA1364", "LA1316", "LA1028", "LA2172", "LA2133", "LA1322", "LA1589", "LA1269", "LA2933", "SL2.50", "LA3475", "LA3124", "LA0429", "LA3909", "LA0436"))

# Determine which plots to create
datasets_for_densitree <- c("Pease2016", "Vanderpool2020")
datasets_for_boxplots <- c("Pease2016", "Vanderpool2020")
datasets_to_compare_tree_topologies <- c("Pease2016", "Vanderpool2020")



#### Step 2: Open files and packages ####
# Open packages
print("Open packages ")
library(ggplot2)
library(cowplot)
library(reshape2) # functions: melt, recast
library(ape) # functions: read.tree, Ntip, root
library(phangorn) # functions: treedist, densiTree
library(phytools) # functions: nodeHeights (to rescale tree height)

#library(ips)
#library(treespace) # phylogenetic tree exploration
#library(adegraphics) # improved graphical functionalities from ade4 (multivariate data analysis)
#library(adegenet) # toolkit for exploring genomic and genetic data
#library(rgl) # for interactive 3D plots

# Source functions
source(paste0(maindir,"code/func_plots.R"))



#### Step 3: Preparation for plotting ####
# Open the treelikeness_df
treelikeness_df <- read.csv(treelikeness_file, stringsAsFactors = FALSE)



#### Step 4: Constructing cloudograms from the filtered loci sets ####
# Iterate through the datasets 
for (dataset in datasets[datasets %in% datasets_for_densitree]){
  # Find the folder for this dataset
  dataset_trees_folder <- paste0(species_tree_folder, dataset, "/species_trees/")
  dataset_files <- list.files(dataset_trees_folder)
  dataset_gene_trees <- paste0(dataset_trees_folder, grep(".txt", dataset_files, value = TRUE))
  
  # Want to compare the pass/fail trees for each test
  tests <- c("PHI", "maxchi", "geneconv", "allTests")
  test = "PHI"
  for (test in tests){
    # Find gene tree sets that passed/failed the test
    pass_test_file <- grep("pass", grep(test, dataset_gene_trees, value = TRUE), value = TRUE)
    fail_test_file <- grep("fail", grep(test, dataset_gene_trees, value = TRUE), value = TRUE)
    
    # Open gene tree sets as multiphylo objects
    pass_trees <- read.tree(pass_test_file)
    fail_trees <- read.tree(fail_test_file)
    
    # Remove trees that do not include every taxa
    # Identify number of tips on tree
    pass_tree_tips <- unlist(lapply(1:length(pass_trees), function(i,pass_trees){Ntip(pass_trees[[i]])}, pass_trees = pass_trees))
    fail_tree_tips <- unlist(lapply(1:length(fail_trees), function(i,fail_trees){Ntip(fail_trees[[i]])}, fail_trees = fail_trees))
    # Find trees with missing tips
    pass_trees_to_remove <- which(pass_tree_tips != n_tips[[dataset]])
    fail_trees_to_remove <- which(fail_tree_tips != n_tips[[dataset]])
    # Remove trees with missing tips
    pass_trees_to_keep <- setdiff(1:length(pass_trees), pass_trees_to_remove)
    pass_trees <- pass_trees[pass_trees_to_keep]
    fail_trees_to_keep <- setdiff(1:length(fail_trees), fail_trees_to_remove)
    fail_trees <- fail_trees[fail_trees_to_keep]
    
    # Root all trees
    pass_trees <- root(pass_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    fail_trees <- root(fail_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    
    # Scale trees to have the same length
    pass_trees <- rescale.multiphylo(pass_trees, 1)
    fail_trees <- rescale.multiphylo(fail_trees, 1)
    
    # Plot trees that passed as a densiTree
    p1 <- densiTree(pass_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
                    direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.01, )
    p2 <- densiTree(fail_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
                    direction = "leftwards", scaleX = TRUE)
    
    # Plot trees as densitrees
    p3 <- ggdensitree(pass_trees[1:10], colour = "steelblue", tip.order = taxa_order[[dataset]], alpha = 0.3, align.tips = TRUE) + geom_tiplab(colour = "black")
    p4 <- ggdensitree(fail_trees[1:10], colour = "steelblue", tip.order = taxa_order[[dataset]], alpha = 0.3)
    # reverse x-axis for p2 and set offset to make the trees face each other
    d4 <- p4$data
    d4$x <- max(d4$x) - d4$x + 1
    p4$data <- d4
    # plot trees next to each other
    trees_plot <- p3 + p4
    multiplot(p3, p4, ncol = 2)
  }
}
