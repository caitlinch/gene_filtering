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

# Create a facetted histogram to compare the distributions of different test statistics
# Plot the neighborNet test statistic and p-value distributions
# Extract the relevant columns of interest and structure properly
e = melt_df[melt_df$variable %in% c("neighbour_net_trimmed","nn_trimmed_sig","neighbour_net_untrimmed","nn_untrimmed_sig"),]
e$group = factor(e$variable,levels = c("neighbour_net_untrimmed","neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig"))
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("neighbour_net_trimmed" = "NeighborNet (trimmed)","nn_trimmed_sig" = "NeighborNet (trimmed)",
                    "neighbour_net_untrimmed" = "NeighborNet (untrimmed)","nn_untrimmed_sig" = "NeighborNet (trimmed) p value")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
# Form the plot
p <- ggplot(e, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~group, scales = "free_y", ncol = 2, labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Tree Proportion Values \n")
ggsave(filename = paste0(results_dir,e$dataset[1],"_treeProportion_histograms.png"), plot = p, units = "in")

#####
rep <- ggplot(ts_df, aes(x = neighbour_net_trimmed, y = X3SEQ_p_value)) + geom_point()
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

# histograms <- want to look at the distribution of scores and pick one that has a useable variation in treelikeness
p <- ggplot(ts_df, aes(x = X3SEQ_prop_recombinant_sequences)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = X3SEQ_p_value)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = neighbour_net_untrimmed)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = nn_untrimmed_sig)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = neighbour_net_trimmed)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = nn_trimmed_sig)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = sCF_mean)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = sCF_mean_sig)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = sCF_median)) + geom_histogram()
p
p <- ggplot(ts_df, aes(x = sCF_median_sig)) + geom_histogram()
p
