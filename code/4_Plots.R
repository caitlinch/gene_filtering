### empirical_treelikeness/4_Plots.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs
# Caitlin Cherryh 2021


##### Step 1: Set the file paths for input and output files, and necessary functions/directories #####
print("set filepaths")
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# plots_dir         <- for saving plots and analyses. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# alignment_dir     <- location of sequence alignments for each dataset (in same order as datasets)
# datasets          <- set name(s) for the dataset(s)

treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
treelikeness_file <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/RecombinationDetection_Vanderpool2020_collated_results_complete_trimmedLoci_trimmedTaxa.csv"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_results/"
alignment_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/01_Data_Vanderpool2020/1730_Alignments_FINAL/"
datasets <- c("Vanderpool2020")


#### Step 2: Open files and packages ####
## Open packages
library(ggplot2) # data visualisation and better plotting
library(reshape2)
library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
                  # these methods all assume an unrooted tree so trees can be used as is for this analysis
library(ips)
#library(treespace) # phylogenetic tree exploration
#library(adegraphics) # improved graphical functionalities from ade4 (multivariate data analysis)
#library(adegenet) # toolkit for exploring genomic and genetic data
#library(rgl) # for interactive 3D plots
## Open functions
source(paste0(maindir,"code/func_plots.R"))


#### Step 3: Preparation for plotting ####
# Give vectors with information about the datasets names 
names(alignment_dir) <- datasets

# Open the treelikeness_df
treelikeness_df <- read.csv(treelikeness_file, stringsAsFactors = FALSE)


#### Step 4: Compare test results with the number of informative sites ####
# Extract the number of parsimony informative sites for each locus
site_info_list <- lapply(treelikeness_df$loci_name, get.DNA.site.info, alignment_dir)
site_info_df <- as.data.frame(do.call(rbind, site_info_list))
# Add new columns to treelikeness_df
treelikeness_df$n_seg_sites <- site_info_df$n_seg_sites
treelikeness_df$n_pis <- site_info_df$n_pis

# Make a dataframe of the test statistic values 
treelikeness_df$bins <- cut(treelikeness_df$n_pis, breaks = c(0,50,100,150,200,250,300,400,500,2500))
windows <- levels(treelikeness_df$bins)
vars <- c("X3SEQ_p_value", "PHI_permutation_p_value", "max_chi_squared_p_value", "NSS_p_value", "geneconv")
proportion_significant_values <- c()
n_loci_in_window <- c()
window_record <- c()
window_start <- c()
var_record <- c()
for (v in vars){
  # Iterate through each variable
  print(v)
  # Iterate through each window
  for (i in 1:length(levels(treelikeness_df$bins))){
    # Record the variable and the window
    var_record <- c(var_record, v)
    window_record <- c(window_record, levels(treelikeness_df$bins)[i])
    window_start_temp <- gsub("\\(", "", strsplit(levels(treelikeness_df$bins)[i], ",")[[1]][1])
    window_start <- c(window_start, window_start_temp)
    # Subset the dataframe by this window
    w <- levels(treelikeness_df$bins)[i]
    w_df <- treelikeness_df[(treelikeness_df$bins == w),]
    print(paste0(w, " : ", nrow(w_df)))
    if (nrow(w_df) == 0){
      temp_prop <- NA
      temp_n_loci_in_window <- 0
    } else {
      if (v == "geneconv"){
        temp_prop <- length(which((w_df$geneconv_inner_fragment_simulated_p_value <= 0.05) & (w_df$geneconv_outer_fragment_simulated_p_value <= 0.05)))/(nrow(w_df))
        temp_n_loci_in_window <- nrow(w_df)
      } else {
        temp_prop <- (length(which(w_df[[v]] <= 0.05)))/(nrow(w_df))
        temp_n_loci_in_window <- nrow(w_df)
      }
    }
    proportion_significant_values <- c(proportion_significant_values, temp_prop)
    n_loci_in_window <- c(n_loci_in_window, temp_n_loci_in_window)
  }
}
pis_df <- data.frame(windows = window_record, 
                       window_start = as.numeric(window_start), 
                       proportion_significant_p_values = proportion_significant_values,
                       n_loci = n_loci_in_window,
                       var = var_record,
                       plot_labels = paste0(window_record, " \n n = ",n_loci_in_window ))

# Make plots of test statistics against number of parsimony informative sites
facet_labels <- c("3SEQ", "PHI", "MaxChi", "NSS", "GeneConv")
names(facet_labels) <- c("X3SEQ_p_value", "PHI_permutation_p_value", "max_chi_squared_p_value", "NSS_p_value", "geneconv")
formula = y~x

ggplot(data = pis_df, aes(x = window_start, y = proportion_significant_p_values)) +
  geom_point(na.rm = TRUE, cex = 3) +
  geom_smooth(method = "lm", formula = formula) +
  ggpmisc::stat_poly_eq(formula = formula, parse = TRUE) +
  facet_grid(var~., labeller = labeller(var = facet_labels)) +
  scale_x_continuous(name = "Number of parsimony informative sites", breaks = pis_df$window_start[1:9], labels = pis_df$plot_labels[1:9]) +
  scale_y_continuous(name = "Proportion of statistically significant p-values (per window)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



