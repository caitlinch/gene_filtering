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
library(RColorBrewer) # For nice colour palettes
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

# Make the color palette for tanglegrams
col_pal <- brewer.pal(8, "Dark2")



#### Step 3: Plotting Metazoans dataset ####
## Find Primate tree files
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = T))
m_tree_files <- grep("Metazoan", tree_files, value = T, ignore.case = T)
m_astral_trees <- grep("ASTRAL", m_tree_files, value = T, ignore.case = T)
m_concat_trees <- grep("CONCAT", m_tree_files, value = T, ignore.case = T)
## Sort file paths for the species trees and the trees that pass the tests
ma_sp_file <- grep("NoTest", m_astral_trees, value = T, ignore.case = T)
map_files <- grep("pass", m_astral_trees, value = T, ignore.case = T)
mc_sp_file <- grep("NoTest", m_concat_trees, value = T, ignore.case = T)
mcp_files <- grep("pass", m_concat_trees, value = T, ignore.case = T)
## Assemble phylo objects
# Open species trees
ma_sp_tree <- read.tree(ma_sp_file)
mc_sp_tree <- read.tree(mc_sp_file)
# Open test trees
map_trees <- lapply(map_files, read.tree)
class(map_trees) <- "multiPhylo"
mcp_trees <- lapply(mcp_files, read.tree)
class(mcp_trees) <- "multiPhylo"
## Add arbitrary terminal branch lengths for ASTRAL trees
ma_sp_tree$edge.length[which(is.na(ma_sp_tree$edge.length))] <- 1
for (i in 1:length(map_trees)){
  temp_tree <- map_trees[[i]]
  terminal_branch_inds <- which(is.na(temp_tree$edge.length))
  temp_tree$edge.length[terminal_branch_inds] <- 1
  map_trees[[i]] <- temp_tree
}
## Prepare parameters for plotting
met_params <- list(map_1 = list(cols = c(rep(col_pal[[1]], 26), rep(col_pal[[2]], 3), rep(col_pal[[3]], 8), rep(col_pal[[4]], 7), rep(col_pal[[5]], 10),
                                         rep(col_pal[[6]], 2), rep(col_pal[[7]], 2), rep(col_pal[[8]], 5), rep(col_pal[[1]], 3), rep(col_pal[[2]], 2), 
                                         rep(col_pal[[3]], 2), rep(col_pal[[4]], 1), rep(col_pal[[5]], 2), rep("grey50", 3)),
                                left_title = "Species tree",
                                right_title = "Pass",
                                main_title = "Metazoan dataset",
                                sub_title = "ASTRAL trees - GENECONV",
                                left_tree = ma_sp_tree,
                                right_tree = map_trees[[1]],
                                output_file_path = paste0(op_dir, "Metazoan_ASTRAL_GENECONV_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ),
                   map_2 = list(cols = c(rep(col_pal[[1]], 18), rep(col_pal[[2]], 2), rep(col_pal[[3]], 3), rep(col_pal[[4]], 3), rep("grey50", 3)),
                                left_title = "Species tree",
                                right_title = "Pass",
                                sub_title = "ASTRAL trees - MaxChi",
                                left_tree = ma_sp_tree,
                                right_tree = map_trees[[2]],
                                output_file_path = paste0(op_dir, "Metazoan_ASTRAL_MaxChi_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ),
                   map_3 = list(cols = c(rep(col_pal[[1]], 18), rep(col_pal[[2]], 2), rep(col_pal[[3]], 3), rep(col_pal[[4]], 3), rep("grey50", 3)),
                                left_title = "Species tree",
                                right_title = "Pass",
                                main_title = "Metazoan dataset",
                                sub_title = "ASTRAL trees - PHI",
                                left_tree = ma_sp_tree,
                                right_tree = map_trees[[3]],
                                output_file_path = paste0(op_dir, "Metazoan_ASTRAL_PHI_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ),
                   mcp_1 = list(cols = "grey50",
                                left_title = "Species tree",
                                right_title = "Pass",
                                main_title = "Metazoan dataset",
                                sub_title = "Concatenated trees - GENECONV",
                                left_tree = ma_sp_tree,
                                right_tree = mcp_trees[[1]],
                                output_file_path = paste0(op_dir, "Metazoan_CONCAT_GENECONV_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ),
                   mcp_2 = list(cols = "grey50",
                                left_title = "Species tree",
                                right_title = "Pass",
                                main_title = "Metazoan dataset",
                                sub_title = "Concatenated trees - MaxChi",
                                left_tree = ma_sp_tree,
                                right_tree = mcp_trees[[2]],
                                output_file_path = paste0(op_dir, "Metazoan_CONCAT_MaxChi_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ),
                   mcp_3 = list(cols = c(rep(col_pal[[1]], 18), rep(col_pal[[2]], 2), rep(col_pal[[3]], 3), rep(col_pal[[4]], 3), rep(col_pal[[5]], 3)),
                                left_title = "Species tree",
                                right_title = "Pass",
                                main_title = "Metazoan dataset",
                                sub_title = "Concatenated trees - PHI",
                                left_tree = ma_sp_tree,
                                right_tree = mcp_trees[[3]],
                                output_file_path = paste0(op_dir, "Metazoan_CONCAT_PHI_tanglegram"),
                                outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") ) )
## Format trees
tt_species <- met_params[["map_1"]]$left_tree
tt_pass <- met_params[["map_1"]]$right_tree
# Root trees
tt_species <- root(tt_species, outgroup = met_params[["map_1"]]$outgroup, resolve.root = TRUE)
tt_pass <- root(tt_pass, outgroup = met_params[["map_1"]]$outgroup, resolve.root = TRUE)
# Make trees ultrametric
tt_species <- force.ultrametric(tt_species, method = "extend")
tt_pass <- force.ultrametric(tt_pass, method = "extend")
# Ladderise trees
tt_species <- ladderize(tt_species, right = TRUE)
tt_pass <- ladderize(tt_pass, right = TRUE)
# Reorder trees
tt_species <- reorder(tt_species, "cladewise")
tt_pass <- reorder(tt_pass, "cladewise")
## Plot trees
# Make tanglegram
tanglegram(tt_species, tt_pass, edge.lwd = 1.5, axes = F,
           main = met_params[["map_1"]]$sub_title,
           main_left = met_params[["map_1"]]$left_title, main_right = met_params[["map_1"]]$right_title, 
           rank_branches = T, match_order_by_labels = T, common_subtrees_color_lines = T, type = "r")

tanglegram(tt_species, tt_pass, color_lines = primate_params[["map_1"]]$cols, edge.lwd = 1.5, axes = F,
           main = met_params[["map_1"]]$sub_title,
           main_left = met_params[["map_1"]]$left_title, main_right = met_params[["map_1"]]$right_title, 
           rank_branches = T, match_order_by_labels = T, common_subtrees_color_lines = T, type = "r")



#### Step 4: Plotting Plants dataset ####





