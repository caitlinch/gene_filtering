### empirical_treelikeness/4_Plots.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
library(ggplot2) # data visualisation and better plotting
library(adegraphics) # improved graphical functionalities from ade4 (multivariate data analysis)
library(adegenet) # toolkit for exploring genomic and genetic data
library(rgl) # for interactive 3D plots
library(reshape2)
library(treespace) # phylogenetic tree exploration
library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
# these methods all assume an unrooted tree so trees can be used as is for this analysis



##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
print("set filepaths")
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# output_dir        <- for collated output and results from treelikeness analysis. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# plots_dir         <- for saving plots and analyses. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# datasets          <- set name(s) for the dataset(s)
# plots             <- output plots - TRUE or FALSE

treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_results/"
datasets <- c("Vanderpool2020")
plots <- FALSE



#### Step 3: Open files
source(paste0(maindir,"code/func_plots.R"))

# Construct the folder names for the output dirs (where the treelikeness data is stored)
output_dirs <- paste0(output_dir,datasets,"/")
names(output_dirs) <- datasets
# Create a set of output folders
plot_dirs <- paste0(plot_dir,datasets,"_plots/")
names(plot_dirs) <- datasets
for (d in plot_dirs){
  if (dir.exists(d) == FALSE){
    dir.create(d)
  }
}


