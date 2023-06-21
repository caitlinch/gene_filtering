### gene_filtering/code/5_Tanglegrams.R
## R program to plot tanglegrams from the species trees
# Caitlin Cherryh 2023

## This script creates tanglegrams from the species trees for the 4 empirical datasets

## ggtree notes:
# To add node labels (for defining clade labels): geom_text(aes(label=node), hjust=-.3)
# To add clade labels (for defining clades with a vertical bar): geom_cladelab(node=34, label="Clade")



##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir                 <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# op_dir                <- for saving tanglegrams and analyses.
# annotations_csv_file    <- location of misc/annotations.csv file from the Leebens-Mack (2019) "Data from 1000 Plants Transcriptomes" data repository
# roots_by_groups         <- set which taxa is outgroup for each dataset: Primates, Tomatoes, Metazoans, and Plants

### Caitlin's paths ###
# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
op_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/tanglegrams/"
annotation_csv_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/misc/annotations.csv"

# Dataset information
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))
### End of Caitlin's paths ###

## Code snippets (might be useful)
# # To extract tips in correct order after ladderising
# is_tip <- tt_pass$edge[,2] <= length(tt_pass$tip.label)
# ordered_tips <- tt_pass$edge[is_tip, 2]
# tip_labels_from_ladder <- tree$tip.label[ordered_tips]



#### Step 2: Open files and packages ####
# Open packages: Phylogenetics
library(ape) # functions: read.tree, Ntip, root
# library(TreeTools) # for CollapseNode function

# Open packages: Plotting
#library(ggplot2) # for nice plots
library(ggtree) # for plotting phylogenetic trees
#library(ggtext) # for nice tree plots
#library(patchwork) # for collating plots
# library(phangorn) # for densiTree function - replaced by ggdensitree
# library(TreeTools) # for CollapseNode function

# Open packages: Tanglegrams
library(dendextend)
# library(phylogram)


# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "species_trees/")



#### Step 3: Plotting Primates dataset ####
# Find Primate tree files
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = TRUE))
primate_tree_files <- grep("Primates", tree_files, value = TRUE)
primate_astral_trees <- grep("ASTRAL", primate_tree_files, value = TRUE)
primate_concat_trees <- grep("CONCAT", primate_tree_files, value = TRUE)
# Sort trees into separate objects for pass/fail
paf_files <- grep("fail", primate_astral_trees, value = T)
pap_files <- grep("pass", primate_astral_trees, value = T)
pcf_files <- grep("fail", primate_concat_trees, value = T)
pcp_files <- grep("pass", primate_concat_trees, value = T)
# Assemble phylo objects
paf_trees <- lapply(paf_files, read.tree)
class(paf_trees) <- "multiPhylo"
pap_trees <- lapply(pap_files, read.tree)
class(pap_trees) <- "multiPhylo"
pcf_trees <- lapply(pcf_files, read.tree)
class(pcf_trees) <- "multiPhylo"
pcp_trees <- lapply(pcp_files, read.tree)
class(pcp_trees) <- "multiPhylo"
# Add arbitrary terminal branch lengths for ASTRAL trees
for (i in 1:length(paf_trees)){
  temp_tree <- paf_trees[[i]]
  terminal_branch_inds <- which(is.na(temp_tree$edge.length))
  temp_tree$edge.length[terminal_branch_inds] <- 1
  paf_trees[[i]] <- temp_tree
}
for (i in 1:length(pap_trees)){
  temp_tree <- pap_trees[[i]]
  terminal_branch_inds <- which(is.na(temp_tree$edge.length))
  temp_tree$edge.length[terminal_branch_inds] <- 1
  pap_trees[[i]] <- temp_tree
}
# Make basic plots for the pass and fail trees

tanglegram(force.ultrametric(pap_trees[[1]], method = "extend"), force.ultrametric(paf_trees[[1]], method = "extend"))

# Extract trees
tt_pass <- pap_trees[[1]]
tt_fail <- paf_trees[[1]]
# Root trees
tt_pass <- root(tt_pass, outgroup = roots_by_group[["Primates"]], resolve.root = TRUE)
tt_fail <- root(tt_fail, outgroup= roots_by_group[["Primates"]], resolve.root = TRUE)
# Make trees ultrametric
tt_pass <- force.ultrametric(tt_pass, method = "extend")
tt_fail <- force.ultrametric(tt_fail, method = "extend")
# Ladderise trees
tt_pass <- ladderize(tt_pass, right = TRUE)
tt_fail <- ladderize(tt_fail, right = TRUE)
# Make tanglegram
tanglegram(tt_pass, tt_fail)


#### Step 4: Plotting Tomatoes dataset ####



#### Step 5: Plotting Metazoans dataset ####



#### Step 6: Plotting Plants dataset ####



