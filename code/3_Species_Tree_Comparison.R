### empirical_treelikeness/code/3_Species_Tree_Comparison.R
## R program to identify the best fitting tree for a candidate dataset
## Additional software packages required:
# 
# Caitlin Cherryh 2021


##### Step 1: Set file paths and run variables #####
run = "local"

if (run == "local"){
  # Datasets of interest
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  tree_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  output_dir <-c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_dataAnalysis/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
}



##### Step 2: Open packages and functions #####
# Open packages
library(ape)



##### Step 3: Compare ASTRAL trees #####
# For each test, I have a fail tree, a pass tree. I also have an all tree - one that includes every loci from the dataset
# I want to take my three trees and modify the format to work for QGOF: add arbitrary terminal branches (1) and remove posterior probability values
# Then I need to convert the list of gene trees to a list of gene trees that work for QGOF
# Working in Julia:
#   - Convert the list of gene trees into a set of quartet CFs
#   - Apply the QGOF test and save the output values
#   - Repeat this for the other two trees

# Starting for one example:

# Extract all the ASTRAL trees for one dataset
dataset = "Vanderpool2020"
dataset_folder <- paste0(output_dir, dataset, "/")
if (dir.exists(dataset_folder) == FALSE){dir.create(dataset_folder)}
# Extract the names of the ASTRAL species trees and the .txt files containing the list of gene trees
species_tree_folder <- paste0(tree_dir, dataset, "/", "species_trees", "/")
all_species_trees_files <- list.files(species_tree_folder)
all_astral_trees <- grep("\\.tre", all_species_trees_files, value = TRUE)
all_astral_gene_trees <- gsub("_species.tre", ".txt", all_astral_trees)
all_astral_files <- c(all_astral_trees, all_astral_gene_trees)

# For each test, create a new folder and copy the relevant files into that folder
tests_to_run <- c("allTests", "PHI", "maxchi", "geneconv")
t = "PHI"
for (t in tests_to_run){
  # Collect files for this test
  files <- c(grep(t, all_astral_trees, value = TRUE), grep("NoTest", all_astral_trees, value = TRUE), 
             grep("pass", grep(t, all_astral_gene_trees, value = TRUE), value = TRUE))
  # Name the new folder for running this test
  new_folder <- paste0(dataset_folder, "quarnetGoFtest_", t, "/")
  if (dir.exists(new_folder) == FALSE){dir.create(new_folder)}
  # Move files into this folder
  for (i in files){
    file.copy(from = paste0(species_tree_folder, i), to = paste0(new_folder, i))
  }
  
  # Give files their new full filepath
  files <- paste0(new_folder, files)
  # Rewrite ASTRAL trees to match format for quarnetGoFtest
  lapply(files[1:3], reformat.ASTRAL.tree.for.Julia)
  
}



