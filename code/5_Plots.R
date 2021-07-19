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

# Software packages



#### Step 2: Open files and packages ####
# Open packages
print("Open packages ")
library(ggplot2)
library(reshape2)
library(ape)
library(phangorn) # treedist, densiTree
library(ggtree)

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
for (dataset in datasets){
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
    
    # Root all trees
    pass_trees <- root(pass_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    fail_trees <- root(fail_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    
    # Plot trees that passed as a densiTree
    densiTree(pass_trees, col = "blue", type = "cladogram")
    densiTree(fail_trees, col = "blue", type = "cladogram")
    ggdensitree(pass_trees)
    
  }
  
  }
