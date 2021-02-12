### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs

# Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
# BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments



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
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_results/"
datasets <- c("Vanderpool2020")
plots <- FALSE



#### Step 3: Open files
# Create a set of output folders
output_dirs <- paste0(output_dir,datasets,"/")
names(output_dirs) <- datasets
plot_dirs <- paste0(plot_dir,datasets,"/")
names(plot_dirs) <- datasets

##### Step 3: Data exploration ####
# Open data
for (dataset in datasets){
  ### Open files and get dataframes ready for plotting
  # Make sure folder for plots and results exists
  if (!dir.exists(plot_dirs[[dataset]])){
    dir.create(plot_dirs[[dataset]])
  }
  
  # Open output from treelikeness analysis
  op_files <- list.files(output_dir)
  results_files <- grep("collated_results", op_files, value = TRUE)
  ds_files <- paste0(output_dir, grep(dataset, results_files, value = TRUE))
  ds_result_file <- grep("melted",ds_files, invert = TRUE, value = TRUE)
  ds_result_file <- grep("trimmedLoci",ds_result_file, invert = TRUE, value = TRUE)
  op_df <- read.csv(ds_result_file, stringsAsFactors = FALSE)
  
  # Calculate derived variables (if they haven't been already calculated)
  # Check this by seeing if columns for these variables already exist
  if (length(which(names(dataset) == "proportion_constant_sites")) > 0){
    op_df$proportion_constant_sites <- round(op_df$num_constant_sites / op_df$n_sites, 3)
    op_df$proportion_invariant_sites <- round(op_df$num_invariant_sites / op_df$n_sites, 3)
    op_df$proportion_informative_sites <- round(op_df$num_parsimony_informative_sites / op_df$n_sites, 3)
    # Save op_df with derived variables
    write.csv(op_df, file = ds_result_file, row.names = FALSE) 
  }
  
  # Make name for melted dataset 
  long_op_name <- gsub(".csv","_melted.csv", ds_result_file)
  # Open long data if it already exists
  if (file.exists(long_op_name)){
    long_df <- read.csv(long_op_name)
  } else {
    # Put data in long format for ggplot
    id_vars <- c("dataset", "loci", "sequence_type", "n_taxa", "n_sites", "num_constant_sites", "proportion_constant_sites", "num_invariant_sites", "proportion_invariant_sites",
                 "num_parsimony_informative_sites", "proportion_informative_sites", "num_site_patterns", "total_tree_length", "sum_of_internal_branch_lengths",
                 "proportion_internal_branches", "substitution_model", "AC_rate", "AG_rate", "AT_rate", "CG_rate", "CT_rate", "GT_rate", "A_freq",
                 "C_freq", "G_freq", "T_freq", "GC_content_mean", "GC_content_variance", "GC_content_sd")
    measure_vars <- c("X3SEQ_prop_recombinant_sequences", "X3SEQ_p_value", "tree_proportion",
                      "tree_proportion_p_value", "sCF_mean","sCF_median")
    long_df <- melt(op_df, id = id_vars, measure.vars = measure_vars)
    write.csv(long_df, file = long_op_name, row.names = FALSE)
  }
  
  # Select which variables to plot
  plotting_vars <- c("X3SEQ_prop_recombinant_sequences", "X3SEQ_p_value", "tree_proportion", "tree_proportion_p_value", "sCF_mean","sCF_median")
  
  ### Plot and save treelikeness histograms
  if (plots == TRUE){
    p <- ggplot(long_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free") + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treelikeness_histogram_freeScales.png") , plot = p)
    p <- ggplot(long_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free_x") + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treelikeness_histogram_freeX.png") , plot = p)
    ggplot(long_df, aes(x = n_sites)) + geom_histogram() + theme_bw()
    
    ## #Plot and save information about the alignments as histograms
    info_df <- melt(op_df, id = c("dataset","loci"), 
                    measure.vars = c("n_taxa", "n_sites", "proportion_constant_sites", "proportion_informative_sites", "num_site_patterns", "total_tree_length",
                                     "proportion_internal_branches", "GC_content_mean", "GC_content_sd")
    )
    info_df[info_df$variable == "proportion_internal_branches",4] <- round(info_df[info_df$variable == "proportion_internal_branches",4]/100,3)
    p <- ggplot(info_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free") + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_alignmentInfo_histogram_freeScales.png") , plot = p)
    p <- ggplot(info_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free_x") + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_alignmentInfo_histogram_freeX.png") , plot = p)
    
    ### Plot some variables against each other to investigate relationships
    # look for a correlation between the 3seq p-value and the proportion of recombinant sequences
    p <- ggplot(data = op_df, aes(x = X3SEQ_prop_recombinant_sequences, y = X3SEQ_p_value)) + geom_point() + theme_bw() + 
      scale_x_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05)) + scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_value_comparison.png") , plot = p)
    
    # look for a correlation between the tree proportion and 3seq p-values
    p <- ggplot(data = op_df, aes(x = tree_proportion_p_value, y = X3SEQ_p_value, colour = tree_proportion)) + geom_point() + theme_bw() + 
      scale_x_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05)) + scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_pvalue_comparison.png") , plot = p)
    
    # look for a correlation between the tree proportion test statistic and p-value
    p <- ggplot(data = op_df, aes(x = tree_proportion, y = tree_proportion_p_value)) + geom_point() + theme_bw() + 
      scale_x_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05)) + scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_value_comparison.png") , plot = p)
    
    # look for a correlation between the tree proportion test and the proportion of internal branches
    p <- ggplot(data = op_df, aes(x = proportion_internal_branches, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_proportionInternalBranches_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = proportion_internal_branches, y = tree_proportion_p_value)) + geom_point() + theme_bw() + 
      scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_proportionInternalBranches_comparison.png") , plot = p)
    
    # look for a correlation between the 3seq and the proportion of internal branches
    p <- ggplot(data = op_df, aes(x = proportion_internal_branches, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_proportionInternalBranches_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = proportion_internal_branches, y = X3SEQ_p_value)) + geom_point() + theme_bw() +
      scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_proportionInternalBranches_comparison.png") , plot = p)
    
    # look for a correlation between the tree proportion test and the proportion of informative sites branches
    p <- ggplot(data = op_df, aes(x = proportion_informative_sites, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_proportionInformativeSites_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = proportion_informative_sites, y = tree_proportion_p_value)) + geom_point() + theme_bw() +
      scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_proportionInformativeSites_comparison.png") , plot = p)
    
    # look for a correlation between the 3seq test and the proportion of informative sites branches
    p <- ggplot(data = op_df, aes(x = proportion_informative_sites, y = X3SEQ_prop_recombinant_sequences, )) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_proportionInformativeSites_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = proportion_informative_sites, y = X3SEQ_p_value)) + geom_point() + theme_bw() +
      scale_y_continuous(breaks = seq(0,1,0.25), minor_breaks = seq(0,1,0.05))
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_proportionInformativeSites_comparison.png") , plot = p)
    
    # look for a correlation between number of sites and tree proportion/3seq
    x <- c(rep("n_sites",4), rep("n_taxa", 4), rep("total_tree_length", 4), rep("num_site_patterns",4), rep("GC_content_mean", 4))
    y <- rep(c("tree_proportion","tree_proportion_p_value","X3SEQ_prop_recombinant_sequences","X3SEQ_p_value"),5)
    x_names <- c(rep("nSites",4), rep("nTaxa", 4), rep("treeLength",4), rep("numSitePatterns",4), rep("GC_content",4))
    y_names <- rep(c("tree_proportion","tree_proportion_p_value","3seq_prop_recombinant_sequences","3seq_p_value"),5)
    for (i in 1:3){
      p <- ggplot(data = op_df, aes(x = .data[[x[i]]], y = .data[[y[i]]])) + geom_point() + theme_bw()
      p_name <- paste0(plot_dirs[[dataset]], dataset, "_", y_names[i], "_", x_names[i], "_comparison.png")
      p
      #ggsave(filename = p_name , plot = p)
    }
    
    p <- ggplot(data = op_df, aes(x = n_sites, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_nSites_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_sites, y = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_nSites_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_sites, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_nSites_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_sites, y = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_nSites_comparison.png") , plot = p)
    
    # look for a correlation between number of taxa and tree proportion/3seq
    p <- ggplot(data = op_df, aes(x = n_taxa, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_nTaxa_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_taxa, y = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_nTaxa_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_taxa, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_nTaxa_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = n_taxa, y = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_nTaxa_comparison.png") , plot = p)
    
    # look for a correlation between tree length and tree proportion/3seq
    p <- ggplot(data = op_df, aes(x = total_tree_length, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_treeLength_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = total_tree_length, y = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_treeLength_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = total_tree_length, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_treeLength_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = total_tree_length, y = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_treeLength_comparison.png") , plot = p)
    
    # look for a correlation between number of site patterns and tree proportion/3seq
    p <- ggplot(data = op_df, aes(x = num_site_patterns, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_NumSitePatterns_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = num_site_patterns, y = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_NumSitePatterns_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = num_site_patterns, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_NumSitePatterns_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = num_site_patterns, y = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_NumSitePatterns_comparison.png") , plot = p)
    
    # look for a correlation between GC content and tree proportion/3seq
    p <- ggplot(data = op_df, aes(x = GC_content_mean, y = tree_proportion, color = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_GC_content_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = GC_content_mean, y = tree_proportion_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_pValue_GC_content_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = GC_content_mean, y = X3SEQ_prop_recombinant_sequences, color = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_GC_content_comparison.png") , plot = p)
    p <- ggplot(data = op_df, aes(x = GC_content_mean, y = X3SEQ_p_value)) + geom_point() + theme_bw()
    ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_pValue_GC_content_comparison.png") , plot = p)
  }
  
  ### Classify trees into 4 groups: 3seq only, tp only, both, neither
  # plot and compare treespace - do trees group together?
  # plot and compare groves (change colour to be based on group!) - do trees group together?
  # Classify loci using p-values as cut off <- significant p-value = treelike
  op_df$X3SEQ_treelike <- as.numeric(op_df$X3SEQ_p_value)
  op_df$X3SEQ_treelike[op_df$X3SEQ_p_value <= 0.05] <- "TREELIKE"
  op_df$X3SEQ_treelike[op_df$X3SEQ_p_value > 0.05] <- "NON-TREELIKE"
  op_df$X3SEQ_treelike <- factor(op_df$X3SEQ_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  op_df$tree_proportion_p_value_treelike <- as.numeric(op_df$tree_proportion_p_value)
  op_df$tree_proportion_p_value_treelike[op_df$tree_proportion_p_value <= 0.05] <- "TREELIKE"
  op_df$tree_proportion_p_value_treelike[op_df$tree_proportion_p_value > 0.05] <- "NON-TREELIKE"
  op_df$tree_proportion_p_value_treelike <- factor(op_df$tree_proportion_p_value_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  
  # For Vanderpool data:
  # When cut-off is 0.7, 972/1730 are treelike
  # When cut-off is 0.75, 580/1730 are treelike
  # Median tree proportion is 0.7116612
  op_df$tree_proportion_treelike <- as.numeric(op_df$tree_proportion)
  op_df$tree_proportion_treelike[op_df$tree_proportion > 0.70] <- "TREELIKE"
  op_df$tree_proportion_treelike[op_df$tree_proportion <= 0.70] <- "NON-TREELIKE"
  op_df$tree_proportion_treelike <- factor(op_df$tree_proportion_treelike, levels = c("NON-TREELIKE","TREELIKE"), ordered = TRUE)
  
  # Sort loci by p-values for both tree proportion and 3seq
  op_df$sorted_p_value <- op_df$loci
  op_df$sorted_p_value[op_df$tree_proportion_p_value_treelike == "TREELIKE" & op_df$X3SEQ_treelike == "TREELIKE"] <- "Both"
  op_df$sorted_p_value[op_df$tree_proportion_p_value_treelike == "NON-TREELIKE" & op_df$X3SEQ_treelike == "TREELIKE"] <- "3seq only"
  op_df$sorted_p_value[op_df$tree_proportion_p_value_treelike == "TREELIKE" & op_df$X3SEQ_treelike == "NON-TREELIKE"] <- "Tree proportion only"
  op_df$sorted_p_value[op_df$tree_proportion_p_value_treelike == "NON-TREELIKE" & op_df$X3SEQ_treelike == "NON-TREELIKE"] <- "Neither"
  op_df$sorted_p_value <- factor(op_df$sorted_p_value, levels = c("Neither","3seq only","Tree proportion only","Both"), ordered = TRUE)
  op_df$sorted_p_value_col <- fac2col(op_df$sorted_p_value, col.pal = grDevices::colorRampPalette(RColorBrewer::brewer.pal(10,"PRGn")))
  
  # Sort loci by p-values for 3seq and test statistic value for tree proportion
  op_df$sorted <- op_df$loci
  op_df$sorted[op_df$tree_proportion_treelike == "TREELIKE" & op_df$X3SEQ_treelike == "TREELIKE"] <- "Both"
  op_df$sorted[op_df$tree_proportion_treelike == "NON-TREELIKE" & op_df$X3SEQ_treelike == "TREELIKE"] <- "3seq only"
  op_df$sorted[op_df$tree_proportion_treelike == "TREELIKE" & op_df$X3SEQ_treelike == "NON-TREELIKE"] <- "Tree proportion only"
  op_df$sorted[op_df$tree_proportion_treelike == "NON-TREELIKE" & op_df$X3SEQ_treelike == "NON-TREELIKE"] <- "Neither"
  op_df$sorted <- factor(op_df$sorted, levels = c("Neither","3seq only","Tree proportion only","Both"), ordered = TRUE)
  op_df$sorted_col <- fac2col(op_df$sorted, col.pal = grDevices::colorRampPalette(RColorBrewer::brewer.pal(10,"PRGn")))
  
  # Now plot the treelikeness test statistics against each other and colour by group
  pretty_colours <- RColorBrewer::brewer.pal(5,"YlGnBu")[2:5]
  dark_colours <- RColorBrewer::brewer.pal(4,"Dark2")
  p <- ggplot(data = op_df, aes(x = tree_proportion_p_value, y = X3SEQ_p_value, color = sorted_p_value)) + geom_point() + theme_bw() + 
    scale_color_manual(labels = c("Neither (n = 779)","3seq only  (n = 587)","Tree proportion only  (n = 121)","Both (n = 223)"), values = dark_colours) + 
    guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test p-value results", subtitle = "Vanderpool 2020") +
    scale_x_continuous(name = "Tree proportion p-value") + scale_y_continuous(name = "3seq p-value")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_sortingLoci_p-values.png") , plot = p)
  # Now plot the tree proportion test statistics against the 3seq p-value and colour by group
  p <- ggplot(data = op_df, aes(x = tree_proportion, y = X3SEQ_p_value, color = sorted)) + geom_point() + theme_bw() + 
    scale_color_manual(labels = c("Neither (n = 378)","3seq only  (n = 380)","Tree proportion only  (n = 542)","Both (n = 430)"), values = dark_colours) + 
    guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test statistic results", subtitle = "Vanderpool 2020") +
    scale_x_continuous(name = "Tree proportion") + scale_y_continuous(name = "3seq p-value")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_sortingLoci.png") , plot = p)
  
  ### Examine treespace
  # Break down trees to those with all species
  all_taxa_df <- op_df[op_df$n_taxa == 29,]
  #all_taxa_df <- all_taxa_df[1:20,] # small dataframe for testing - comment/uncomment as needed
  # Read in all trees at once
  mp <- read.tree(text = all_taxa_df$tree)
  # use treespace <- have to pick a method that works on unrooted trees. Do this by picking one that assumes trees are unrooted
  # "Warning message: In is.euclid(distmat) : Zero distance(s)" may appear when doing "RF" method <- this means you have duplicate rows in your distance matrix
  # Basically means two trees have the exact same values so the function thinks that you have a duplicate row
  res <- treespace(mp, nf=3, method = "wRF")
  # Use ggplot to plot the treespace using colours to sort loci as above
  # plotGroves(res$pco, lab.show=FALSE) # plotting treespace using plotGroves function = simple but hard to add colours for factors
  # Extract principal component vectors and add the sorted column so you know how to group the taxa
  pc_df <- res$pco$tab
  pc_df$sorted <- all_taxa_df$sorted
  pc_df$sorted_p_value <- all_taxa_df$sorted_p_value
  # sorted_p_value
  p <- ggplot(pc_df,aes(x = A1, y = A2, color = sorted_p_value)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace: A2 ~ A1", subtitle = "Grouped by 3seq and tree proportion p-values") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_pValue_PC_A1A2.png") , plot = p)
  p <- ggplot(pc_df,aes(x = A2, y = A3, color = sorted_p_value)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A2") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace: A3 ~ A2", subtitle = "Grouped by 3seq and tree proportion p-values") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_pValue_PC_A2A3.png") , plot = p)
  p <- ggplot(pc_df,aes(x = A1, y = A3, color = sorted_p_value)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace: A3 ~ A2", subtitle = "Grouped by 3seq and tree proportion p-values") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_pValue_PC_A1A3.png") , plot = p)
  # sorted
  p <- ggplot(pc_df,aes(x = A1, y = A2, color = sorted)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace: A2 ~ A1", subtitle = "Grouped by 3seq p-value and tree proportion") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_PC_A1A2.png") , plot = p)
  p <- ggplot(pc_df,aes(x = A2, y = A3, color = sorted)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A2") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace: A3 ~ A2", subtitle = "Grouped by 3seq p-value and tree proportion") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_PC_A2A3.png") , plot = p)
  p <- ggplot(pc_df,aes(x = A1, y = A3, color = sorted)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace: A3 ~ A2", subtitle = "Grouped by 3seq p-value and tree proportion") + 
    scale_color_brewer(palette = "Dark2")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_treelikeness_PC_A1A3.png") , plot = p)
  
  # Now find the groves
  groves <- findGroves(res, nclust=5)
  # Plot using plotGroves to see how the groves are organised
  plotGroves(groves)
  # Create a data frame to use for plotting
  groves_pc_df <- groves$treespace$pco$li
  groves_pc_df$groups <- groves$groups
  # Make and save the plots
  p <- ggplot(groves_pc_df,aes(x = A1, y = A2, color = groups)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Groves")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace with groves: A2 ~ A1") + 
    scale_color_brewer(palette = "Paired")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_groves_PC_A1A2.png") , plot = p)
  p <- ggplot(groves_pc_df,aes(x = A2, y = A3, color = groups)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Groves")) +
    scale_x_continuous(name = "A2") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace with groves: A3 ~ A2") + 
    scale_color_brewer(palette = "Paired")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_groves_PC_A2A3.png") , plot = p)
  p <- ggplot(groves_pc_df,aes(x = A1, y = A3, color = groups)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Groves")) +
    scale_x_continuous(name = "A1") + scale_y_continuous(name = "A3") + labs(title = "Vanderpool 2020 treespace with groves: A3 ~ A2") + 
    scale_color_brewer(palette = "Paired")
  ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treespace_groves_PC_A1A3.png") , plot = p)
  # Plot and save a 3D plot
  colours <- fac2col(groves$groups, col.pal=colorRampPalette(RColorBrewer::brewer.pal(5,"Paired")))
  plot3d(groves$treespace$pco$li[,1],
         groves$treespace$pco$li[,2],
         groves$treespace$pco$li[,3],
         col=colours, type="s", size=1.5,
         xlab="", ylab="", zlab="")
  rgl.snapshot(paste0(plot_dirs[[dataset]],dataset,"_treespace_groves_3dplot_snapshot.png"), fmt = "png")
  
  
}


