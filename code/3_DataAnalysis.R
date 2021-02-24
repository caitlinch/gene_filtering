### empirical_treelikeness/3_DataAnalysis.R
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
# datasets          <- set name(s) for the dataset(s)

treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
datasets <- c("Vanderpool2020")



#### Step 3: Open files
source(paste0(maindir,"code/func_plots.R"))

# Construct the folder names for the output dirs (where the treelikeness data is stored)
output_dirs <- paste0(output_dir,datasets,"/")
names(output_dirs) <- datasets


##### Step 3: Data exploration ####
# Open data
for (dataset in datasets){
  ### Open files and get dataframes ready for plotting
  # Make sure folder for plots and results exists
  
  # Open output from treelikeness analysis
  # Dataframes in this file contain treelikeness data as follows: 
  #   * t29_df <- all loci with 29 taxa, wide format
  #   * t29_long_df <- all loci with 29 taxa, long format
  #   * all_loci_df <- all loci, wide format
  #   * all_loci_long_df <- all loci, long format
  op_files <- list.files(output_dir)
  results_files <- grep("collated_results", op_files, value = TRUE)
  ds_files <- paste0(output_dir, grep(dataset, results_files, value = TRUE))
  ds_result_file <- grep("melted",ds_files, invert = TRUE, value = TRUE)
  trimmed_result_file <- grep("trimmedLoci",ds_result_file, value = TRUE)
  t29_df <- read.csv(trimmed_result_file, stringsAsFactors = FALSE) 
  all_loci_file <- grep("trimmedLoci",ds_result_file, invert = TRUE, value = TRUE)
  all_loci_df <- read.csv(all_loci_file, stringsAsFactors = FALSE)
  all_loci_long_op_name <- gsub(".csv","_melted.csv", all_loci_file)
  all_loci_long_df <- read.csv(all_loci_long_op_name, stringsAsFactors = FALSE)
  
  # Calculate derived variables (if they haven't been already calculated)
  # Check this by seeing if columns for these variables already exist
  if (length(which(names(t29_df) == "proportion_constant_sites")) == 0){
    t29_df$proportion_constant_sites <- round(t29_df$num_constant_sites / t29_df$n_sites, 3)
    t29_df$proportion_invariant_sites <- round(t29_df$num_invariant_sites / t29_df$n_sites, 3)
    t29_df$proportion_informative_sites <- round(t29_df$num_parsimony_informative_sites / t29_df$n_sites, 3)
    # Save t29_df with derived variables
    write.csv(t29_df, file = trimmed_result_file, row.names = FALSE) 
  }
  
  # Make name for melted dataset 
  long_trimmed_op_name <- gsub(".csv","_melted.csv", trimmed_result_file)
  # Open long data if it already exists
  if (file.exists(long_trimmed_op_name)){
    t29_long_df <- read.csv(long_trimmed_op_name)
  } else {
    # Put data in long format for ggplot
    id_vars <- c("dataset", "loci", "sequence_type", "n_taxa", "n_sites", "num_constant_sites", "proportion_constant_sites", "num_invariant_sites", "proportion_invariant_sites",
                 "num_parsimony_informative_sites", "proportion_informative_sites", "num_site_patterns", "total_tree_length", "sum_of_internal_branch_lengths",
                 "proportion_internal_branches", "substitution_model", "AC_rate", "AG_rate", "AT_rate", "CG_rate", "CT_rate", "GT_rate", "A_freq",
                 "C_freq", "G_freq", "T_freq", "GC_content_mean", "GC_content_variance", "GC_content_sd")
    measure_vars <- c("X3SEQ_prop_recombinant_sequences", "X3SEQ_p_value", "tree_proportion",
                      "tree_proportion_p_value", "sCF_mean","sCF_median")
    t29_long_df <- melt(t29_df, id = id_vars, measure.vars = measure_vars)
    write.csv(t29_long_df, file = long_trimmed_op_name, row.names = FALSE)
  }
}


### Classify trees into 4 groups: 3seq only, tp only, both, neither
## plot and compare treespace - do trees group together?
## plot and compare groves (change colour to be based on group!) - do trees group together?
# Start by adding columns to the dataframe that classify each value as either TREELIKE or NON-TREELIKE
t29_df <- classify.treelikeness.statistics(t29_df, 0.7)

# Now plot the treelikeness test statistics against each other and colour by group
pretty_colours <- RColorBrewer::brewer.pal(5,"YlGnBu")[2:5]
dark_colours <- RColorBrewer::brewer.pal(4,"Dark2")
ggplot(data = t29_df, aes(x = tree_proportion_p_value, y = X3SEQ_p_value, color = sorted_p_value)) + geom_point() + theme_bw() + 
  scale_color_manual(labels = c("Neither (n = 779)","3seq only  (n = 587)","Tree proportion only  (n = 121)","Both (n = 223)"), values = dark_colours) + 
  guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test p-value results", subtitle = "Vanderpool 2020") +
  scale_x_continuous(name = "Tree proportion p-value") + scale_y_continuous(name = "3seq p-value")
# Now plot the tree proportion test statistics against the 3seq p-value and colour by group
ggplot(data = t29_df, aes(x = tree_proportion, y = X3SEQ_p_value, color = sorted)) + geom_point() + theme_bw() + 
  scale_color_manual(labels = c("Neither (n = 378)","3seq only  (n = 380)","Tree proportion only  (n = 542)","Both (n = 430)"), values = dark_colours) + 
  guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test statistic results", subtitle = "Vanderpool 2020") +
  scale_x_continuous(name = "Tree proportion") + scale_y_continuous(name = "3seq p-value")

### Examine treespace
# Break down trees to those with all species
t29_df <- t29_df[t29_df$n_taxa == 29,]
#t29_df <- t29_df[1:20,] # small dataframe for testing - comment/uncomment as needed
# Read in all trees at once
mp <- read.tree(text = t29_df$tree)
# use treespace <- have to pick a method that works on unrooted trees. Do this by picking one that assumes trees are unrooted
# "Warning message: In is.euclid(distmat) : Zero distance(s)" may appear when doing "RF" method <- this means you have duplicate rows in your distance matrix
# Basically means two trees have the exact same values so the function thinks that you have a duplicate row
res <- treespace(mp, nf=3, method = "wRF")
# Use ggplot to plot the treespace using colours to sort loci as above
# plotGroves(res$pco, lab.show=FALSE) # plotting treespace using plotGroves function = simple but hard to add colours for factors
# Extract principal component vectors and add the sorted column so you know how to group the taxa
pc_df <- res$pco$tab
pc_df$sorted <- t29_df$sorted
pc_df$sorted_p_value <- t29_df$sorted_p_value
# sorted_p_value plot
ggplot(pc_df,aes(x = A1, y = A2, color = sorted_p_value)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
  scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace: A2 ~ A1", subtitle = "Grouped by 3seq and tree proportion p-values") + 
  scale_color_brewer(palette = "Dark2")


# Now find the groves
groves <- findGroves(res, nclust=5)
# Plot using plotGroves to see how the groves are organised
plotGroves(groves)
# Create a data frame to use for plotting
groves_pc_df <- groves$treespace$pco$li
groves_pc_df$groups <- groves$groups
# plot
ggplot(groves_pc_df,aes(x = A1, y = A2, color = groups)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Groves")) +
  scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace with groves: A2 ~ A1") + 
  scale_color_brewer(palette = "Paired")

# Plot and save a 3D plot
colours <- fac2col(groves$groups, col.pal=colorRampPalette(RColorBrewer::brewer.pal(5,"Paired")))
plot3d(groves$treespace$pco$li[,1],
       groves$treespace$pco$li[,2],
       groves$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")




