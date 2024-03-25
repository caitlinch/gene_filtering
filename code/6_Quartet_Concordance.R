### gene_filtering/code/6_Quartet_Concordance.R
## R program to plot and explore results of the gene filtering project
# Caitlin Cherryh 2024

## This script processes and plots ASTRAL trees with quartet concordance factors


##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
# maindir     <- "gene_filtering" repository location (github.com/caitlinch/gene_filtering)
# plot_dir    <- for saving plots and analyses.
# qcf_dir     <- location of ASTRAL trees with qCFs

### Caitlin's paths ###
# Folders and filepaths
maindir   <- "/Users/caitlincherryh/Documents/Repositories/gene_filtering/"
plot_dir  <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/06_results/plots/"
qcf_dir   <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_qCF/"
### End of Caitlin's paths ###



#### Step 2: Open files and packages ####
# Open packages
library(ape) # read.tree, Ntip, root
library(phangorn) # as.splits

# Source functions
source(paste0(maindir, "code/func_comparison.R"))

# Dataset information
dataset_names <- c("Plants" = "1KP", "Tomatoes" = "Pease2016", "Metazoa" = "Whelan2017", "Primates" = "Vanderpool2020")
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoa" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))



#### Step 3: Convert trees with branch labels into a convenient dataframe ####
# Extract all tree files
all_files <- list.files(qcf_dir)
tree_files <- paste0(qcf_dir, grep(".tre", all_files, value = T))
# Iterate through one dataset at a time
d <- "Pease2016"


# Identify the original dataset tree
original_tree_path <- grep("NoTest", grep(d, tree_files, value = T), value = T)

# Identify the cleaned trees
cleaned_tree_files <- grep("pass", grep(d, tree_files, value = T), value = T)
recomb_tree_files  <- grep("fail", grep(d, tree_files, value = T), value = T)

# For each of the cleaned trees
i = 1
i_tree_file <- cleaned_tree_files[i]

# Identify which splits occur in both trees



