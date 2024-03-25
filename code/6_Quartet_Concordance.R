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
library(ape) # functions: read.tree, Ntip, root

# Source functions


# Dataset information
roots_by_group <- list("Plants" = c("BAKF", "ROZZ", "MJMQ", "IRZA", "IAYV", "BAJW", "APTP", "LXRN", "NMAK", "RFAD", "LLEN", "RAPY", "OGZM",
                                    "QDTV", "FIDQ", "EBWI", "JQFK", "BOGT", "VKVG", "DBYD", "FSQE", "LIRF", "QLMZ", "JCXF", "ASZK", "ULXR",
                                    "VRGZ", "LDRY", "VYER", "FIKG", "RWXW", "FOMH", "YRMA", "HFIK", "JGGD"), 
                       "Metazoan" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), 
                       "Primates" = c("Mus_musculus"), 
                       "Tomatoes" = c("LA4116", "LA2951", "LA4126"))



#### Step 3: Convert trees with branch labels into a convenient dataframe ####









