### gene_filtering/5_Plots.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Caitlin Cherryh 2022

## This script requires:
#     - DensiTree (Bouckaert 2010) (https://www.cs.auckland.ac.nz/~remco/DensiTree/)

## This script:
# 1. Creates a variety of plots and figures



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                 <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir                <- for saving plots and analyses. This file should contain a folder for each dataset (where the folder name and corresponding dataset name are identical)
# species_tree_folder     <- folder containing the species trees estimated in ASTRAL and IQ-Tree
# treelikeness_file       <- file containing results of recombination detection tests
# datasets                <- set name(s) for the dataset(s)
# roots                   <- set which taxa is outgroup for each dataset
# n_tips                  <- number of tips in each tree
# taxa_order              <- order of taxa for plotting (if desired)
# densitree_path          <- path to DensiTree software executable for making DensiTree figures
# datasets_for_densitree  <- which datasets you want to create figures for using DensiTree

### Caitlin's paths ###
# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/"
species_tree_folder <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
treelikeness_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/02_AllDatasets_collated_RecombinationDetection_TrimmedLoci.csv"

# Dataset information
datasets <- c("Vanderpool2020", "Pease2016", "Whelan2017", "1KP")
roots <- list(c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
              c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
              c("Mus_musculus"), 
              c("LA4116", "LA2951", "LA4126"))
n_tips <- c("Pease2016" = 29, "Vanderpool2020" = 29)
taxa_order <- list("Pease2016" = c("LA4116", "LA4126", "LA2951", "LA3778", "LA0716", "LA1777", "LA0407", "LA4117", "LA1782", "LA2964", 
                                   "LA0444", "LA0107", "LA1358", "LA2744", "LA1364", "LA1316", "LA1028", "LA2172", "LA2133", "LA1322",
                                   "LA1589", "LA1269", "LA2933", "SL2.50", "LA3475", "LA3124", "LA0429", "LA3909", "LA0436"),
                   "Vanderpool2020" = c("Mus_musculus", "Tupaia_chinensis", "Galeopterus_variegatus", "Otolemur_garnettii",
                                        "Propithecus_coquereli", "Microcebus_murinus", "Carlito_syrichta", "Callithrix_jacchus",
                                        "Aotus_nancymaae", "Cebus_capucinus_imitator", "Saimiri_boliviensis", "Colobus_angolensis_palliatus",
                                        "Piliocolobus_tephrosceles", "Rhinopithecus_bieti", "Rhinopithecus_roxellana", "Chlorocebus_sabaeus",
                                        "Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada", "Macaca_nemestrina",
                                        "Macaca_fascicularis", "Macaca_mulatta", "Nomascus_leucogenys", "Pongo_abelii", "Gorilla_gorilla",
                                        "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"))
# DensiTree paths/variables
densitree_path <- "/Users/caitlincherryh/Documents/Executables/Densitree/DensiTree.v2.2.7.jar"
datasets_for_densitree <- c("Pease2016", "Vanderpool2020")
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ggplot2)
library(ape) # functions: read.tree, Ntip, root
library(phangorn) # functions: treedist, densiTree
library(phytools) # functions: nodeHeights (to rescale tree height)
# Source functions
source(paste0(maindir,"code/func_plots.R"))



#### Step 3: Preparation for plotting ####
# Open the treelikeness_df
treelikeness_df <- read.csv(treelikeness_file, stringsAsFactors = FALSE)



#### Step 4: Constructing cloudograms from the filtered loci sets ####
# Make a folder to store the cloudograms in 
cloudogram_folder <- paste0(plot_dir, "cloudograms/")
if (dir.exists(cloudogram_folder) == FALSE){dir.create(cloudogram_folder)}

# Iterate through the datasets 
for (dataset in datasets[datasets %in% datasets_for_densitree]){
  # Find the folder for this dataset
  dataset_trees_folder <- paste0(species_tree_folder, dataset, "/species_trees/")
  dataset_files <- list.files(dataset_trees_folder)
  dataset_gene_trees <- paste0(dataset_trees_folder, grep(".txt", dataset_files, value = TRUE))
  
  # Plot a densitree of the NoTest gene trees for this dataset
  # Get all gene trees from dataset and open as a multiphylo object
  none_file <- grep("NoTest", dataset_gene_trees, value = TRUE)
  none_trees <- read.tree(none_file)
  # Remove any trees that are missing taxa
  none_tree_tips <- unlist(lapply(1:length(none_trees), function(i,trees){Ntip(trees[[i]])}, trees = none_trees))
  none_trees_to_remove <- which(none_tree_tips != n_tips[[dataset]])
  if (length(none_trees_to_remove) > 0){
    none_trees_to_keep <- setdiff(1:length(none_trees), none_trees_to_remove)
    none_trees <- none_trees[none_trees_to_keep]
  }
  # Reroot the trees and standardise the length
  none_trees <- root(none_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
  none_trees <- rescale.multiphylo(none_trees, 1)
  # Plot a densitree for the NoTest gene trees
  png(filename = paste0(cloudogram_folder, dataset, "_NoTest_Trees_cloudogram.png"), width = 500, height = 500)
  densiTree(none_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
            direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
  dev.off()
  cairo_pdf(filename = paste0(cloudogram_folder, dataset, "_NoTest_Trees_cloudogram.pdf"))
  densiTree(none_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
            direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
  dev.off()
  
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
    pass_tree_tips <- unlist(lapply(1:length(pass_trees), function(i,trees){Ntip(trees[[i]])}, trees = pass_trees))
    fail_tree_tips <- unlist(lapply(1:length(fail_trees), function(i,trees){Ntip(trees[[i]])}, trees = fail_trees))
    # Find trees with missing tips
    pass_trees_to_remove <- which(pass_tree_tips != n_tips[[dataset]])
    fail_trees_to_remove <- which(fail_tree_tips != n_tips[[dataset]])
    # Remove trees with missing tips
    if (length(pass_trees_to_remove) > 0){
      pass_trees_to_keep <- setdiff(1:length(pass_trees), pass_trees_to_remove)
      pass_trees <- pass_trees[pass_trees_to_keep]
    }
    if (length(fail_trees_to_remove) > 0){
      fail_trees_to_keep <- setdiff(1:length(fail_trees), fail_trees_to_remove)
      fail_trees <- fail_trees[fail_trees_to_keep]
    }
    
    # Root all trees
    pass_trees <- root(pass_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    fail_trees <- root(fail_trees, outgroup = roots[[dataset]], resolve.root = TRUE)
    
    # Scale trees to have the same length
    pass_trees <- rescale.multiphylo(pass_trees, 1)
    fail_trees <- rescale.multiphylo(fail_trees, 1)
    
    # Plot gene trees that passed test as a densiTree
    png(filename = paste0(cloudogram_folder, dataset, "_", test, "_passTrees_cloudogram.png"), width = 500, height = 500)
    densiTree(pass_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
              direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
    dev.off()
    cairo_pdf(filename = paste0(cloudogram_folder, dataset, "_", test, "_passTrees_cloudogram.pdf"))
    densiTree(pass_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
              direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
    dev.off()
    # Plot gene trees that failed test as a densiTree
    png(filename = paste0(cloudogram_folder, dataset, "_", test, "_failTrees_cloudogram.png"), width = 500, height = 500)
    densiTree(fail_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
              direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
    dev.off()
    cairo_pdf(filename = paste0(cloudogram_folder, dataset, "_", test, "_failTrees_cloudogram.pdf"))
    densiTree(fail_trees, col = "steelblue", type = "cladogram", alpha = 0.03, scale.bar = FALSE, consensus = taxa_order[[dataset]], 
              direction = "rightwards", scaleX = TRUE, adj = 1, label.offset = 0.02, cex = 1)
    dev.off()
  }
}
