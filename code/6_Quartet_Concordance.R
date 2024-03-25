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
tree_files <- grep(".tre", all_files, value = T)
# Create dataframe for testing every pair of trees
clean_trees       <- grep("pass", tree_files, value = T)
recomb_trees      <- grep("fail", tree_files, value = T)
alignment_trees   <- grep("NoTest", tree_files, value = T)
tree_comp_df      <- data.frame(dataset = c(rep("1KP", 8),
                                            rep("Pease2016", 8),
                                            rep("Vanderpool2020",  8),
                                            rep("Whelan2017", 3) ),
                                clean_tree = c("1KP_allTests_pass_ASTRAL_qCF.tre", "1KP_allTests_pass_ASTRAL_qCF.tre",
                                               "1KP_geneconv_pass_ASTRAL_qCF.tre", "1KP_geneconv_pass_ASTRAL_qCF.tre",
                                               "1KP_maxchi_pass_ASTRAL_qCF.tre", "1KP_maxchi_pass_ASTRAL_qCF.tre",
                                               "1KP_PHI_pass_ASTRAL_qCF.tre" , "1KP_PHI_pass_ASTRAL_qCF.tre",
                                               "Pease2016_allTests_pass_ASTRAL_qCF", "Pease2016_allTests_pass_ASTRAL_qCF",
                                               "Pease2016_geneconv_pass_ASTRAL_qCF", "Pease2016_geneconv_pass_ASTRAL_qCF",
                                               "Pease2016_maxchi_pass_ASTRAL_qCF", "Pease2016_maxchi_pass_ASTRAL_qCF",
                                               "Pease2016_PHI_pass_ASTRAL_qCF", "Pease2016_PHI_pass_ASTRAL_qCF",
                                               "Vanderpool2020_allTests_pass_ASTRAL_qCF", "Vanderpool2020_allTests_pass_ASTRAL_qCF",
                                               "Vanderpool2020_geneconv_pass_ASTRAL_qCF", "Vanderpool2020_geneconv_pass_ASTRAL_qCF",
                                               "Vanderpool2020_maxchi_pass_ASTRAL_qCF", "Vanderpool2020_maxchi_pass_ASTRAL_qCF",
                                               "Vanderpool2020_PHI_pass_ASTRAL_qCF", "Vanderpool2020_PHI_pass_ASTRAL_qCF",
                                               "Whelan2017_geneconv_pass_ASTRAL_qCF.tre", "Whelan2017_maxchi_pass_ASTRAL_qCF.tre",
                                               "Whelan2017_PHI_pass_ASTRAL_qCF.tre" ),
                                comparison_tree = c("1KP_allTests_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                    "1KP_geneconv_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                    "1KP_maxchi_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                    "1KP_PHI_fail_ASTRAL_qCF.tre", "1KP_NoTest_ASTRAL_qCF.tre",
                                                    "Pease2016_allTests_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                    "Pease2016_geneconv_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                    "Pease2016_maxchi_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                    "Pease2016_PHI_fail_ASTRAL_qCF.tre",  "Pease2016_NoTest_ASTRAL_qCF.tre",
                                                    "Vanderpool2020_allTests_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                    "Vanderpool2020_geneconv_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                    "Vanderpool2020_maxchi_fail_ASTRAL_qCF.tre", "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                    "Vanderpool2020_PHI_fail_ASTRAL_qCF.tre" , "Vanderpool2020_NoTest_ASTRAL_qCF.tre",
                                                    "Whelan2017_NoTest_ASTRAL_qCF.tre", "Whelan2017_NoTest_ASTRAL_qCF.tre",
                                                    "Whelan2017_NoTest_ASTRAL_qCF.tre"),
                                tree_directory = qcf_dir,
                                comparison_id = c("allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                  "allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                  "allTests_fail", "noTest", "geneconv_fail", "noTest", "maxchi_fail", "noTest", "PHI_fail", "noTest",
                                                  "noTest", "noTest", "noTest"))
tree_comp_df$clean_id <- unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){paste0(x[[2]], "_", x[[3]])}))
tree_comp_df$recombination_test <-  unlist(lapply(strsplit(tree_comp_df$clean_tree, "_"), function(x){x[[2]]}))
tree_comp_df$comparison_gene_status <- unlist(lapply(strsplit(tree_comp_df$comparison_tree, "_"), function(x){x[[3]]}))
tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "ASTRAL")] <- "unfiltered"
tree_comp_df$comparison_gene_status[which(tree_comp_df$comparison_gene_status == "fail")] <- "recombinant"

# Identify splits in each tree
qcf_op <- lapply(1:nrow(tree_comp_df), compare.splits.wrapper, df = tree_comp_df)





