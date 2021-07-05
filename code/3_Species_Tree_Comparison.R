### empirical_treelikeness/code/3_Species_Tree_Comparison.R
## R program to identify the best fitting tree for a candidate dataset
## Additional software packages required:
# 
# Caitlin Cherryh 2021


##### Step 1: Set file paths and run variables #####
run = "local"

if (run == local){
  input_names <- c("1KP", "Strassert2021","Vanderpool2020", "Pease2016")
  
  # File and directory locations
  alignment_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_1KP/alignments/alignments-FAA-masked_genes/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Strassert2021/02_trimAL_Divvier_filtered_genes_only/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/",
                     "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Pease2016/all_window_alignments/")
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  output_dir <- c("/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/")
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
}
