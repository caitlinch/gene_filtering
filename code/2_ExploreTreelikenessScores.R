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
treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is



##### Step 3: Reformat dataframe #####
# Open data
ts_df <- read.csv(results_path, stringsAsFactors = FALSE)

# create a new column (% of parsimony informative sites)
ts_df$prop_parsimony_informative_sites = ts_df$num_parsimony_informative_sites/ts_df$n_sites

# Reformat data into long form
id_vars <- c("dataset","loci","alignment_file","n_taxa","n_sites","total_tree_depth","num_parsimony_informative_sites",
             "prop_parsimony_informative_sites","GC_content_mean","GC_content_variance")
measure_vars <- c("X3SEQ_num_recombinant_triplets","X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences",
                  "X3SEQ_p_value","neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean","sCF_median",
                  "x3seq_numRecomSeq_sig","sCF_mean_sig","sCF_median_sig")
melt_df <- melt(ts_df, id = id_vars, measure.vars = measure_vars)
# Output melted data as a csv file
output_name <- gsub(".csv","_melted.csv",results_path)
write.csv(melt_df, file = output_name, row.names = FALSE)



##### Step 4: Plot histograms #####
# Create a facetted histogram to compare the distributions of different test statistics
# Plot the neighborNet test statistic and p-value distributions
# Extract the relevant columns of interest (for tree proportion) and structure properly
e = melt_df[melt_df$variable %in% c("neighbour_net_trimmed","nn_trimmed_sig","neighbour_net_untrimmed","nn_untrimmed_sig"),]
e$group = factor(e$variable,levels = c("neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("neighbour_net_trimmed" = "NeighborNet (trimmed)","nn_trimmed_sig" = "NeighborNet (trimmed) p value",
                    "neighbour_net_untrimmed" = "NeighborNet (untrimmed)","nn_untrimmed_sig" = "NeighborNet (untrimmed) p value")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
# Form the plots
p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, scale = "free_y", ncol = 2, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_treeProportion_freey_histograms.png"), plot = p, units = "in")

p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, ncol = 2, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_treeProportion_samey_histograms.png"), plot = p, units = "in")

# Repeat the plots using the existing test statistics (sCF and 3seq)
e = melt_df[melt_df$variable %in% c("X3SEQ_prop_recombinant_sequences","X3SEQ_p_value","sCF_mean","sCF_mean_sig","sCF_median","sCF_median_sig"),]
e$group = factor(e$variable,levels = c("X3SEQ_prop_recombinant_sequences","sCF_mean","sCF_median","X3SEQ_p_value","sCF_mean_sig","sCF_median_sig"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences","X3SEQ_p_value" = "3SEQ (inbuilt) p value",
                    "sCF_mean" = "Mean sCF","sCF_mean_sig" = "Mean sCF p value", "sCF_median" = "Median sCF","sCF_median_sig" = "Median sCF p value")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
# Form the plots
# Plot histograms of existing test statistics
p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, scale = "free", ncol = 3, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_existingStats_freey_histograms.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = value)) +
  geom_histogram(bins=20) +
  facet_wrap(~group, scale = "free_x", ncol = 3, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Value \n") +
  ylab("\n Count \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_existingStats_samey_histograms.png"), plot = p, units = "cm")



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
}

# Plot all the predictors
p <- ggplot(e, aes(x = n_sites, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Number of sites \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_numSites_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = n_taxa, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Number of taxa \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_numTaxa_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = total_tree_depth, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Total tree depth \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_treeDepth_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = prop_parsimony_informative_sites, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Proportion of parsimony informative sites \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_propParsimonyInformativeSites_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = num_parsimony_informative_sites, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Number of parsimony informative sites \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_numParsimonyInformativeSites_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = GC_content_mean, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Mean GC content \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_gcContentMean_points.png"), plot = p, units = "cm")

p <- ggplot(e, aes(x = GC_content_mean, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 2) +
  scale_x_continuous(name = "\n Mean GC content \n") +
  ylab("\n Test statistic value \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_gcContentVar_points.png"), plot = p, units = "cm")

##### Step 6: Plot 3SEQ inbuilt p value against potential predictors #####
# Plot the 3SEQ inbuilt test statistic p values against those calculated using a parametric bootstrap on the proportion of recombinant sequences
p <- ggplot(ts_df, aes(x = X3SEQ_p_value, y = x3seq_numRecomSeq_sig)) + geom_point() +  
  scale_x_continuous(minor_breaks = seq(0, 1, 0.05), name = "3SEQ (inbuilt) p value") + 
  scale_y_continuous(minor_breaks = seq(0, 1, 0.05), name = "\n Number of recombinant sequences p value (parametric bootstrap) \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_3seq_pvalues_examination.png"), plot = p, units = "cm")
