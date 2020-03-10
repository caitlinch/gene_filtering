### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to plot and explore results of the treelikeness simulations on empirical data
# Final result is a reformatted csv file and a number of graphs

# Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
# BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments



##### Step 1: Open packages #####
library(ggplot2) # data visualisation and better plotting
library(reshape2) # restructure data between long and wide formats
library(treespace) # phylogenetic tree exploration
library(adegenet) # toolkit for exploring genomic and genertic data
library(adegraphics) # improved graphical functionalities from ade4 (multivariate data analysis)
library(rgl) # for interactive 3D plots
library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
                    # these methods all assume an unrooted tree so trees can be used as is for this analysis




##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
results_path <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/06_results/treelikeness_scores/Wu_2018_dnaLoci_Primates_completeResults.csv"
results_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/06_results/treelikeness_scores/"
plots_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/06_results/treelikeness_scores/plots/"
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
e = melt_df[melt_df$variable %in% c("X3SEQ_prop_recombinant_sequences","neighbour_net_trimmed","sCF_mean","sCF_median"),]
e$group = factor(e$variable,levels = levels(e$variable)[c(3,6,9,10)])
# Set up the labeller so that the facet names will be formatted nicely (rather than using variable names in the plots)
facet_names <- list("X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences", "neighbour_net_trimmed" = "NeighborNet (trimmed)", 
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
  ggsave(filename = paste0(plots_dir,dataset,"_predictor_",name_additions[[i]],"_testStatistic_points.png"), plot = p, units = "cm")
}

# Repeat this process for p values
# Prepare the dataframe
e = melt_df[melt_df$variable %in% c("X3SEQ_p_value","x3seq_propRecomSeq_sig","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean_sig","sCF_median_sig"),]
e$group = factor(e$variable, levels(e$variable)[c(4,12,7,8,13,14)])
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
    scale_y_continuous("\n p value \n")
  ggsave(filename = paste0(plots_dir,dataset,"_predictor_",name_additions[[i]],"_pValue_points.png"), plot = p, units = "cm")
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
  var_name <- facet_names[variable]
  return(var_name)
}

# Set up the variables to plot and the labels for axes
plot_vars <- c("X3SEQ_p_value", "X3SEQ_num_distinct_recombinant_sequences","X3SEQ_prop_recombinant_sequences","neighbour_net_untrimmed",
               "neighbour_net_trimmed","nn_untrimmed_sig","nn_trimmed_sig","sCF_mean","sCF_median","sCF_mean_sig","sCF_median_sig")
label_vars <- list("X3SEQ_p_value" = "3SEQ (inbuilt) p value", "X3SEQ_num_distinct_recombinant_sequences" = "Number of  \n recombinant sequences",
                   "X3SEQ_prop_recombinant_sequences" = "Proportion of  \n recombinant sequences", "neighbour_net_untrimmed" = "NeighborNet (untrimmed)",
                   "neighbour_net_trimmed" = "NeighborNet (trimmed)", "nn_untrimmed_sig" = "NeighborNet (untrimmed) p value", 
                   "nn_trimmed_sig" = "NeighborNet (trimmed) p value", "sCF_mean" = "Mean sCF","sCF_mean_sig" = "Mean sCF p value", "sCF_median" = "Median sCF",
                   "sCF_median_sig" = "Median sCF p value")
dataset = ts_df$dataset[1]

# Plot all of the test statistics/p values for a single dataset againt the predictors at once
for (var_to_plot in plot_vars){
  p <- ggplot(e, aes(x = value, y = e[[var_to_plot]])) +
    geom_point() +
    facet_wrap(~group,scales = "free", labeller = labeller(group = facet_labeller), nrow = 2, ncol = 4) +
    scale_x_continuous("\n Predictor value \n") +
    scale_y_continuous(name = paste0("\n ",label_vars[var_to_plot][[1]]," \n"))
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



##### Step 8: Exploring treespace #####
# to open a tree: read.tree(text=t) if t is a newick tree or read.tree(file=f) if f is a filepath to a newick tree
t <- read.tree(text="((PONGO_ABE:0.00678774,((PAN_TRO:0.00097927,PAN_PAN:-0.00097827):0.00498971,HOMO_SAP:-0.00024584):0.00721851):0.00341861,
               ((SAIMI_BOL:0.00017936,CALLI_JAC:0.03474404):0.02185356,(CHLOR_SAB:0.05385828,(PAPIO_ANU:0.01147204,
               (MACAC_MUL:0.00000050,MACAC_FAS:0.00000050):0.00231226):0.00574556):0.01722714):0.01398251,((TARSI_SYR:0.12918976,
               (OTOLE_GAR:0.08554176,(MICRO_MUR:0.06490312,DAUBE_MAD:-0.03564042):0.02627855):0.00424960):0.04311964,
               (NOMAS_LEU:0.06105755,GORIL_GOR:0.05080475):0.00548585):0.04554407);")
# to make t a multiphylo (as joining a class phylo object and class multiphylo object doesn't work: you need this to join one tree to an existing multiphylo object):
t2<-list(t)
class(t2)<-"multiPhylo"
t2
# to read in all trees at once (from the newick tree column in the df):
m <- read.tree(text=ts_df$newick_tree)

# Try a treespace using just the trees with all 16 species
# get the trees with all 16 species of primate from the dataset
e = ts_df[(ts_df$n_taxa == 16),]
# Problem trees are 79, 120 
# These trees claim to have 16 taxa but only have 14 when you plot them out
e_trees <- e$newick_tree[c(1:78,80:120,122:128,130:157,159:194,196:305)] # pick 300 trees you know work 
e_tl <- e$neighbour_net_trimmed[c(1:78,80:120,122:128,130:157,159:194,196:305)] # get tree proportion values for those scores
m <- read.tree(text = e_trees)
m100 <- read.tree(text = e_trees)[1:100]
tl100 <- e$neighbour_net_trimmed[c(1:78,80:120,122:128,130:157,159:194,196:305)][1:100] # get tree proportion values for those scores
seqp100 <- e$X3SEQ_p_value[c(1:78,80:120,122:128,130:157,159:194,196:305)][1:100] # get 3seq p values for those scores
seqts100 <- e$X3SEQ_prop_recombinant_sequences[c(1:78,80:120,122:128,130:157,159:194,196:305)][1:100] # get 3seq test statistic values for those scores
scf100 <- e$sCF_mean[c(1:78,80:120,122:128,130:157,159:194,196:305)][1:100] # get sCF values for those scores

# Error in treespace(m, nf = 3) : Tree 79 has different tip labels from the first tree.
t1 <- read.tree(text = e$newick_tree[1])
t79 <- read.tree(text = e$newick_tree[79]) # lists 16 taxa but only has 14...
t121 <- read.tree(text = e$newick_tree[121]) # lists 16 taxa but only has 14...
t129 <- read.tree(text = e$newick_tree[129]) # lists 16 taxa but only has 15...
t158 <- read.tree(text = e$newick_tree[158]) # lists 16 taxa but only has 15...
t195 <- read.tree(text = e$newick_tree[195]) # lists 16 taxa but only has 15...

# Identifying why some trees list 16 species but have only 14 tips...
num_taxa <- ts_df$n_taxa


# use treespace <- have to pick a method that works on unrooted trees. Do this by picking one that assumes trees are unrooted (sneaky!)
# "Warning message: In is.euclid(distmat) : Zero distance(s)" may appear when doing "RF" method <- this means you have duplicate rows in your distance matrix
# Basically means two trees have the exact same values so the function thinks that you have a duplicate row

res <- treespace(m, nf=3, method = "patristic")
res100 <- treespace(m100, nf=3, method = "wRF")
names(res100)
str(res100)
table.image(res100$D, nclass=30) # make a table of the 100 trees. Black = more different. 
# table.value with some customization
table.value(res100$D, nclass=10, method="color", col=funky(10))
# plot trees using first 2 PCs
plotGroves(res100$pco, lab.show=TRUE, lab.cex=1.5)
plotGrovesD3(res100$pco, treeNames=1:100)
groves <- findGroves(res100, nclust=7)
# basic plot
plotGrovesD3(groves)
# alternative with improved legend and tooltip text, giving the tree numbers:
plotGrovesD3(groves, tooltip_text=paste0("Tree ",1:100), legend_width=50, col_lab="Cluster")
# plot axes 2 and 3. This helps to show why, for example, clusters 2 and 4 have been identified as separate, despite them appearing to overlap when viewing axes 1 and 2.
plotGrovesD3(groves, xax=2, yax=3, tooltip_text=paste0("Tree ",1:100), legend_width=50, col_lab="Cluster")
# we can also plot in 3D
# prepare a colour palette:
colours <- fac2col(groves$groups, col.pal=funky)
plot3d(groves$treespace$pco$li[,1],
       groves$treespace$pco$li[,2],
       groves$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")

# to incorporate treelikeness value into the 3D 
# Trying to attach the tree proportion scores to the spheres
groups100 <- tl100
names(groups100) <- 1:100
groves100 <- list("groups" = groups100, "treespace" = res100)
# results in colouring by the treelikeness values
plotGrovesD3(groves100, tooltip_text=paste0("Tree ",1:100), legend_width=50, col_lab="TL value") 


# Ways to get your numbers into colours for the plots
library(rgl)
# Define colours (there may be a better way to do this ...)
# Step 1: use colourRamp to make a spectrum from 2 colours using your number vector
# Step 2: divide by 255 so each colour is on scale from 0 to 1
# Step 3: convert matrix of rgb values to matrix of hex colours row by row
cvec <- apply((colorRamp(c("red","yellow","blue"))(tl100))/255,1, function(x) {rgb(x[1],x[2],x[3])})
# OR Step 1: use viridis to generate colours from treelikeness vector
library(viridis)
vir_col <- viridis(7)
vir_vec <- apply((colorRamp(c(vir_col))(tl100))/255,1, function(x) {rgb(x[1],x[2],x[3])})
# OR Step 1: or using any2col to convert numbers to colours
any2col_op <- any2col(tl100)
any2col_col <- any2col_op$col #extract the colours from the any2col output
# Plot in 3d using the colours that were just scaled and applied
plot3d(groves100$treespace$pco$li[,1],
       groves100$treespace$pco$li[,2],
       groves100$treespace$pco$li[,3],
       col=vir_vec, type="s", size=1.5,
       xlab="", ylab="", zlab="")

# 3 seq p value
vir_vec <- apply((colorRamp(c(vir_col))(seqp100))/255,1, function(x) {rgb(x[1],x[2],x[3])})
plot3d(groves100$treespace$pco$li[,1],
       groves100$treespace$pco$li[,2],
       groves100$treespace$pco$li[,3],
       col=vir_vec, type="s", size=1.5,
       xlab="", ylab="", zlab="")

# 3 seq prop recombinant seq
vir_vec <- apply((colorRamp(c(vir_col))(seqts100))/255,1, function(x) {rgb(x[1],x[2],x[3])})
plot3d(groves100$treespace$pco$li[,1],
       groves100$treespace$pco$li[,2],
       groves100$treespace$pco$li[,3],
       col=vir_vec, type="s", size=1.5,
       xlab="", ylab="", zlab="")

# scf values
vir_vec <- apply((colorRamp(c(vir_col))(scf100/100))/255,1, function(x) {rgb(x[1],x[2],x[3])})
plot3d(groves100$treespace$pco$li[,1],
       groves100$treespace$pco$li[,2],
       groves100$treespace$pco$li[,3],
       col=vir_vec, type="s", size=1.5,
       xlab="", ylab="", zlab="")


