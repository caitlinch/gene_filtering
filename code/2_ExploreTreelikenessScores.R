### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to plot and explore results of the treelikeness simulations on empirical data
# Final result is a reformatted csv file and a number of graphs

# Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
# BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments



##### Step 1: Open packages #####
library(ggplot2)
library(reshape2)



##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
results_path <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/results/treelikeness_scores/Wu_2018_dnaLoci_Primates_completeResults.csv"
results_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/results/treelikeness_scores/"
plots_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/results/treelikeness_scores/plots/"
treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is



##### Step 3: Reformat dataframe #####
# Open data
ts_df <- read.csv(results_path, stringsAsFactors = FALSE)

# create a new column (% of parsimony informative sites)
ts_df$prop_parsimony_informative_sites = ts_df$num_parsimony_informative_sites/ts_df$n_sites

# Reformat data into long form
id_vars <- c("dataset","loci","alignment_file","n_taxa","n_sites","total_tree_depth","num_parsimony_informative_sites",
             "prop_parsimony_informative_sites","GC_content_mean","GC_content_variance", "GC_content_sd")
measure_vars <- c("X3SEQ_num_recombinant_triplets","X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences",
                  "X3SEQ_p_value","neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean","sCF_median",
                  "x3seq_numRecomSeq_sig","x3seq_propRecomSeq_sig","sCF_mean_sig","sCF_median_sig")
melt_df <- melt(ts_df, id = id_vars, measure.vars = measure_vars)
# Output melted data as a csv file
output_name <- gsub(".csv","_melted.csv",results_path)
write.csv(melt_df, file = output_name, row.names = FALSE)



##### Step 4: Plot histograms #####
# Create a facetted histogram to compare the distributions of different test statistics
# Plot histogram with all test statistics and p_values
e = melt_df[melt_df$variable %in% c("neighbour_net_trimmed","nn_trimmed_sig","neighbour_net_untrimmed","nn_untrimmed_sig",
                                    "X3SEQ_prop_recombinant_sequences","X3SEQ_p_value","sCF_mean","sCF_mean_sig","sCF_median","sCF_median_sig"),]
e$group = factor(e$variable,levels = c("neighbour_net_untrimmed","nn_untrimmed_sig","neighbour_net_trimmed","nn_trimmed_sig",
                                       "X3SEQ_prop_recombinant_sequences","X3SEQ_p_value","sCF_mean","sCF_mean_sig","sCF_median","sCF_median_sig"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("neighbour_net_trimmed" = "NeighborNet (trimmed)","nn_trimmed_sig" = "NeighborNet (trimmed) \n p value",
                    "neighbour_net_untrimmed" = "NeighborNet \n (untrimmed)","nn_untrimmed_sig" = "NeighborNet \n (untrimmed) p value",
                    "X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences","X3SEQ_p_value" = "3SEQ (inbuilt) p value",
                    "sCF_mean" = "Mean sCF","sCF_mean_sig" = "Mean sCF p value", "sCF_median" = "Median sCF","sCF_median_sig" = "Median sCF p value")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
  return(variable)
}

p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, scale = "free", ncol = 5, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(plots_dir,e$dataset[1],"_histogram_all_freey.png"), plot = p, units = "in")

p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, scale = "free_x", ncol = 5, nrow = 2, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(plots_dir,e$dataset[1],"_histogram_all_samey.png"), plot = p, units = "in")



##### Step 5: Plot potential predictors of treelikeness against test statistic values #####
# Plot the number of sites against the treelikeness test statistics
# Prepare the dataframe
e = melt_df[melt_df$variable %in% c("neighbour_net_trimmed","X3SEQ_prop_recombinant_sequences","sCF_mean","sCF_median"),]
e$group = factor(e$variable,levels = c("neighbour_net_trimmed","X3SEQ_prop_recombinant_sequences","sCF_mean","sCF_median"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("neighbour_net_trimmed" = "NeighborNet (trimmed)", "X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences",
                    "sCF_mean" = "Mean sCF", "sCF_median" = "Median sCF")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
  return(variable)
}
# Prepare the variables for each plot
x_axis <- c("n_sites","n_taxa","total_tree_depth","prop_parsimony_informative_sites","num_parsimony_informative_sites","GC_content_mean",
            "GC_content_variance","GC_content_sd")
x_axis_names <- c("Number of sites","Number of taxa","Total tree depth","Proportion of parsimony informative sites","Number of parsimony informative sites",
                  "Mean GC content","Variance in GC content","Standard deviation of GC content")
name_additions <- c("numSites","numTaxa","treeDepth","propParsimonyInformativeSites","numParsimonyInformativeSites","gcContentMean","gcContentVar",
                    "gcContentSD")
# Automate the plotting of all the predictors
for (i in 1:length(x_axis)){
  p <- ggplot(e, aes(x = e[x_axis[i]], y = e["value"])) + geom_point() +  
    facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
    scale_x_continuous(name = x_axis_names[i]) + 
    scale_y_continuous("\n Test statistic value \n")
  ggsave(filename = paste0(plots_dir,dataset,"_predictor_",name_additions[[i]],"testStatistic_points.png"), plot = p, units = "cm")
}

# Repeat this process for p values
# Prepare the dataframe
e = melt_df[melt_df$variable %in% c("X3SEQ_p_value","x3seq_propRecomSeq_sig","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean_sig","sCF_median_sig"),]
e$group = factor(e$variable,levels = c("X3SEQ_p_value","x3seq_propRecomSeq_sig","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean_sig","sCF_median_sig"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("X3SEQ_p_value" = "3SEQ (inbuilt) p value","x3seq_propRecomSeq_sig" = "Proportion of recombinant sequences \n p value (parametric bootstrap)",
                    "nn_untrimmed_sig" = "NeighborNet (untrimmed) p value","nn_trimmed_sig" = "NeighborNet (trimmed) p value",
                    "sCF_mean_sig" = "Mean sCF p value","sCF_median_sig" = "Median sCF p value")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
  return(variable)
}
# Prepare the variables for each plot
x_axis <- c("n_sites","n_taxa","total_tree_depth","prop_parsimony_informative_sites","num_parsimony_informative_sites","GC_content_mean",
            "GC_content_variance","GC_content_sd")
x_axis_names <- c("Number of sites","Number of taxa","Total tree depth","Proportion of parsimony informative sites","Number of parsimony informative sites",
                  "Mean GC content","Variance in GC content","Standard deviation of GC content")
name_additions <- c("numSites","numTaxa","treeDepth","propParsimonyInformativeSites","numParsimonyInformativeSites","gcContentMean","gcContentVar",
                    "gcContentSD")
# Automate the plotting of all the predictors
for (i in 1:length(x_axis)){
  p <- ggplot(e, aes(x = e[x_axis[i]], y = e["value"])) + geom_point(size=0.5) +  
    facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
    scale_x_continuous(name = x_axis_names[i]) + 
    scale_y_continuous("\n Test statistic value \n")
  ggsave(filename = paste0(plots_dir,dataset,"_predictor_",name_additions[[i]],"pValue_points.png"), plot = p, units = "cm")
}

##### Step 6: Plot test statistics against potential predictors #####
# Melt the original df using ts values and p values as the id_var (so you can use it on the x axis!)
# Therefore, make the potential predictors the measure_vars (so you can compare them all on the one figure)

# Create a new melted dataframe, where the id variables are the test statistics and p values rather than the predictors
seq_id_vars <- c("X3SEQ_num_recombinant_triplets","X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences",
                 "x3seq_numRecomSeq_sig","X3SEQ_p_value", "neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig",
                 "nn_trimmed_sig","sCF_mean","sCF_median","sCF_mean_sig","sCF_median_sig")
seq_measure_vars <- c("n_taxa","n_sites","total_tree_depth","num_parsimony_informative_sites","prop_parsimony_informative_sites",
                      "GC_content_mean","GC_content_variance", "GC_content_sd")
seq_df <- melt(ts_df, id = seq_id_vars, measure.vars = seq_measure_vars)

# Filter the new dataframe
e = seq_df[seq_df$variable %in% c("n_taxa","n_sites","total_tree_depth","num_parsimony_informative_sites","prop_parsimony_informative_sites",
                                  "GC_content_mean","GC_content_variance","GC_content_sd"),]
e$group = factor(e$variable,levels = c("n_taxa","n_sites","num_parsimony_informative_sites","prop_parsimony_informative_sites",
                                       "total_tree_depth","GC_content_mean","GC_content_variance","GC_content_sd"))

# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("n_taxa" = "Number of taxa","n_sites" = "Number of sites","num_parsimony_informative_sites" = "Num. parsimony \n informative sites",
                      "prop_parsimony_informative_sites" = "Prop. parsimony \n informative sites","total_tree_depth" = "Tree depth",
                      "GC_content_mean" = "GC content (mean)","GC_content_variance" = "GC content (variance)","GC_content_sd" = "GC content (st dev)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

# Set up the variables to plot and the labels for axes
plot_vars <- c("X3SEQ_p_value", "X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences","neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig",
           "nn_trimmed_sig","sCF_mean","sCF_median","sCF_mean_sig","sCF_median_sig")
label_vars <- list("X3SEQ_p_value" = "3SEQ (inbuilt) p value", "X3SEQ_num_distinct_recombinant_sequences" = "Number of  \n recombinant sequences",
                   "X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences", "neighbour_net_untrimmed" = "NeighborNet (untrimmed)",
                   "neighbour_net_trimmed" = "NeighborNet (trimmed)", "nn_untrimmed_sig" = "NeighborNet (untrimmed) p value", 
                   "nn_trimmed_sig" = "NeighborNet (trimmed) p value", "sCF_mean" = "Mean sCF","sCF_mean_sig" = "Mean sCF p value", "sCF_median" = "Median sCF",
                   "sCF_median_sig" = "Median sCF p value")
dataset = ts_df$dataset[1]

# Plot all of the test statistics/p values for a single dataset againt the predictors at once
for (var_to_plot in plot_vars){
  p <- ggplot(e, aes(x = e[[var_to_plot]], y = value)) +
    geom_point() +
    facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 4) +
    scale_x_continuous(name = paste0("\n ",label_vars[var_to_plot][[1]]," \n")) +
    ylab("\n Predictor value \n")
  plot_filename <- paste0(plots_dir,dataset,"_testStatistic_",var_to_plot,"_points.png")
  ggsave(filename = plot_filename, plot = p, units = "cm")
}



##### Step 7: Explore relationships between test statistics #####
# Plot a bunch of variables against each other to see if you encounter anything interesting
# This is automated so you need to collect the x and y axis, the names for these axes and a suitable file name for each plot
# filepath_addition is not the file name - it's a chunk of text added inside the file name
x_axis <- c("X3SEQ_p_value","nn_trimmed_sig","neighbour_net_trimmed","neighbour_net_trimmed","X3SEQ_prop_recombinant_sequences","X3SEQ_prop_recombinant_sequences")
y_axis <- c("x3seq_propRecomSeq_sig","X3SEQ_p_value","X3SEQ_p_value","X3SEQ_prop_recombinant_sequences","nn_trimmed_sig","X3SEQ_p_value")
x_axis_name <- c("3SEQ (inbuilt) p value","NeighborNet (trimmed) p value","NeighborNet (trimmed)","NeighborNet (trimmed)","Proportion of recombinant sequences",
                 "Proportion of recombinant sequences")
y_axis_name <- c("Proportion of recombinant sequences p value (parametric bootstrap)","3SEQ (inbuilt) p value","3SEQ (inbuilt) p value",
                 "Proportion of recombinant sequences","NeighborNet (trimmed) p value","3SEQ (inbuilt) p value")
filepath_addition <- c("3seq_pvalues","NeighborNetTrimmedSig_3seq","NeighborNetTrimmed_3seq","NeighborNetTrimmed_ProportionRecombinantSequences",
                       "NeighborNetTrimmedPValue_ProportionRecombinantSequences","3seq_pValue_ProportionRecombinantSequences")
# Automate the plotting
for (i in 1:length(x_axis)){
  p <- ggplot(ts_df, aes(x = ts_df[x_axis[i]], y = ts_df[y_axis[i]])) + geom_point() +  
    scale_x_continuous(minor_breaks = seq(0, 1, 0.05), name = x_axis_name[i]) + 
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05), name = y_axis_name[i])
  ggsave(filename = paste0(plots_dir,dataset,"_comparison_",filepath_addition[[i]],"_points.png"), plot = p, units = "cm")
}



