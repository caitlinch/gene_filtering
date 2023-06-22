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
# roots_by_groups         <- set which taxa is outgroup for each dataset: Primates, Tomatoes, Metazoans, and Plants

### Caitlin's paths ###
# Folders and filepaths
maindir <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
op_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/tanglegrams/"

# Dataset information
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages: Phylogenetics
library(ape) # functions: read.tree, Ntip, root

# Open packages: Tanglegrams
library(dendextend)

# Open packages: Plotting
# library(ggtree) # for plotting phylogenetic trees
library(viridisLite) # for color palettes


# Source functions
source(paste0(maindir,"code/func_plots.R"))

# Assemble folder for species trees
species_tree_folder <- paste0(maindir, "species_trees/")



#### Step 3: Plotting Metazoans dataset ####
## Find Metazoan tree files
tree_files <- paste0(species_tree_folder, list.files(species_tree_folder, recursive = T))
metazoan_tree_files <- grep("Metazoan", tree_files, value = T, ignore.case = T)
## Plot all tanglegrams
metazoan_plotting_parameters <- get.Metazoan.plotting.parameters(metazoan_tree_files = metazoan_tree_files)
all_list_vars <- names(metazoan_plotting_parameters)
lapply(all_list_vars, plot.Metazoan.tanglegram, metazoan_plotting_parameters)



#### Step 4: Plotting Plants dataset ####





