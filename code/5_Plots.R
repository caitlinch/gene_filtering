### gene_filtering/5_Plots.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2022

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

# Dataset information
datasets <- c("Vanderpool2020", "Pease2016", "Whelan2017", "1KP")
roots <- list("1KP" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                        "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                        "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
              "Whelan2017" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
              "Vanderpool2020" = c("Mus_musculus"), 
              "Pease2016" = c("LA4116", "LA2951", "LA4126"))
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
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ggplot2)
library(ggtree) 
library(ape) # functions: read.tree, Ntip, root
library(phangorn) # functions: treedist, densiTree
library(phytools) # functions: nodeHeights (to rescale tree height)
# Source functions
source(paste0(maindir,"code/func_plots.R"))



#### Step 3: Plotting Primates dataset ####
primate_tree_folder <- paste0(species_tree_folder, "Vanderpool2020/species_trees/")
primate_tree_files <- paste0(primate_tree_folder, list.files(primate_tree_folder, recursive = TRUE))
primate_tree_files <- grep("old_geneconv|Old_geneconv|Old_Geneconv|00_|zz_", primate_tree_files, invert = TRUE, value = TRUE)
primate_astral_trees <- grep("species\\.tre", primate_tree_files, value = TRUE)
primate_concat_trees <- grep("nex.contree", primate_tree_files, value = TRUE)

## Primate Plot 1: ASTRAL No Test
# Open Primates ASTRAL No Test tree
primate_astral_no_test_tree_file <- grep("NoTest", primate_astral_trees, value = TRUE)
primate_astral_no_test_tree <- read.tree(file = primate_astral_no_test_tree_file)
# Reroot tree
primate_astral_no_test_tree <- root(primate_astral_no_test_tree, outgroup = roots[["Vanderpool2020"]])


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
ggtree(pt1) + geom_tiplab(offset = 0, geom = "text", size = 6) + geom_rootedge(rootedge = 0.1) + xlim(-0.1, 1.2)
ggtree(pt2) + geom_tiplab(offset = 0, geom = "text", size = 6) + geom_rootedge(rootedge = 0.1) + xlim(-0.1, 1.2)
ggtree(pt3) + geom_tiplab(offset = 0, geom = "text", size = 6) + geom_rootedge(rootedge = 0.1) + xlim(-0.1, 1.2)


## Primate Plot 3: CONCAT No Test

## Primate Plot 4: CONCAT variable clade (Papionini)


