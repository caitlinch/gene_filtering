### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to examine and explore treelikeness scores 
## Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments

# Caitlin Cherryh 2019

# Script to extract useful information from alignments and explore results of treelikeness scores.

# Open packages
library(ggplot2)
library(reshape2)

# Specify file paths
results_path <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/results/treelikeness_scores/Wu_2018_dnaLoci_Primates_completeResults.csv"
results_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/results/treelikeness_scores/"
treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is

# Open data
ts_df <- read.csv(results_path, stringsAsFactors = FALSE)

# Reformat data into long form
id_vars <- c("dataset","loci","alignment_file","n_taxa","n_sites")
measure_vars <- c("X3SEQ_num_recombinant_triplets","X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences",
                  "X3SEQ_p_value","neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean","sCF_median",
                  "x3seq_numRecomSeq_sig","sCF_mean_sig","sCF_median_sig")
melt_df <- melt(ts_df, id = id_vars, measure.vars = measure_vars)
output_name <- gsub(".csv","_melted.csv",results_path)
write.csv(melt_df, file = output_name, row.names = FALSE)

##### Plot Histograms #####
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
# Form the plot
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

##### Plot potential predictors of treelikeness #####
# Plot the 3SEQ inbuilt test statistic p values against those calculated using a parametric bootstrap on the proportion of recombinant sequences
p <- ggplot(ts_df, aes(x = X3SEQ_p_value, y = x3seq_numRecomSeq_sig)) + geom_point() +  
  scale_x_continuous(minor_breaks = seq(0, 1, 0.05), name = "3SEQ (inbuilt) p value") + 
  scale_y_continuous(minor_breaks = seq(0, 1, 0.05), name = "\n Number of recombinant sequences p value (parametric bootstrap) \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_3seq_pvalues_examination.png"), plot = p, units = "cm")

# Plot the number of sites against the treelikeness test statistics
rep <- ggplot(ts_df, aes(x = nn_trimmed_sig, y = X3SEQ_p_value)) + geom_point()
p

p <- ggplot(ts_df, aes(x = nn_trimmed_sig, y = X3SEQ_p_value)) + geom_point()
p

p <- ggplot(ts_df, aes(x = neighbour_net_trimmed, y = X3SEQ_prop_recombinant_sequences)) + geom_point()
p                  

p <- ggplot(ts_df, aes(x = n_sites, y = neighbour_net_trimmed)) + geom_point()
p

p <- ggplot(ts_df, aes(x = n_sites, y = X3SEQ_p_value)) + geom_point()
p

p <- ggplot(ts_df, aes(x = n_taxa, y = neighbour_net_trimmed)) + geom_point()
p

p <- ggplot(ts_df, aes(x = n_taxa, y = X3SEQ_p_value)) + geom_point()
p


