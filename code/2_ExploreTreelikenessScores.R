### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs

# Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
# BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments



##### Step 1: Open packages #####
library(ggplot2) # data visualisation and better plotting
library(reshape2)
library(treespace) # phylogenetic tree exploration
library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
                  # these methods all assume an unrooted tree so trees can be used as is for this analysis

# library(adegenet) # toolkit for exploring genomic and genetic data
# library(adegraphics) # improved graphical functionalities from ade4 (multivariate data analysis)
# library(rgl) # for interactive 3D plots




##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
print("set filepaths")
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# output_dir        <- for collated output and results from treelikeness analysis. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# plots_dir         <- for saving plots and analyses. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# datasets          <- set name(s) for the dataset(s)

treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code isresults_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
plot_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_results/"
datasets <- c("Vanderpool2020")



#### Step 3: Open files
# Create a set of output folders
output_dirs <- paste0(output_dir,datasets,"/")
names(output_dirs) <- datasets
plot_dirs <- paste0(plot_dir,datasets,"/")
names(plot_dirs) <- datasets

##### Step 3: Data exploration ####
# Open data
for (dataset in datasets){
  # Make sure folder for plots and results exists
  if (!dir.exists(plot_dirs[[dataset]])){
    dir.create(plot_dirs[[dataset]])
  }
  
  # Open output from treelikeness analysis
  op_files <- list.files(output_dir)
  results_files <- grep("collated_results", op_files, value = TRUE)
  ds_result_file <- paste0(output_dir, grep(dataset, results_files, value = TRUE)[1])
  op_df <- read.csv(ds_result_file, stringsAsFactors = FALSE)
  
  # Calculate derived variables
  op_df$proportion_constant_sites <- round(op_df$num_constant_sites / op_df$n_sites, 3)
  op_df$proportion_invariant_sites <- round(op_df$num_invariant_sites / op_df$n_sites, 3)
  op_df$proportion_informative_sites <- round(op_df$num_parsimony_informative_sites / op_df$n_sites, 3)
  # Save op_df with derived variables
  write.csv(op_df, file = ds_result_file, row.names = FALSE)
  
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
  
  # Plot and save treelikeness histograms
   p <- ggplot(long_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free") + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treelikeness_histogram_freeScales.png") , plot = p)
   p <- ggplot(long_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free_x") + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treelikeness_histogram_freeX.png") , plot = p)
   ggplot(long_df, aes(x = n_sites)) + geom_histogram() + theme_bw()
   
   # Plot and save information about the alignments
   info_df <- melt(op_df, id = c("dataset","loci"), 
                   measure.vars = c("n_taxa", "n_sites", "proportion_constant_sites", "proportion_informative_sites", "num_site_patterns", "total_tree_length",
                                    "proportion_internal_branches", "GC_content_mean", "GC_content_sd")
                   )
   info_df[info_df$variable == "proportion_internal_branches",4] <- round(info_df[info_df$variable == "proportion_internal_branches",4]/100,3)
   p <- ggplot(info_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free") + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_alignmentInfo_histogram_freeScales.png") , plot = p)
   p <- ggplot(info_df, aes(x = value)) + geom_histogram() + facet_wrap(variable ~ ., scales = "free_x") + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_alignmentInfo_histogram_freeX.png") , plot = p)

  # plot some variables against each other to investigate relationships
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
   p <- ggplot(data = op_df, aes(x = tree_proportion, y = proportion_internal_branches, color = tree_proportion_p_value)) + geom_point() + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_proportionInternalBranches_comparison.png") , plot = p)
   # look for a correlation between the 3seq and the proportion of internal branches
   p <- ggplot(data = op_df, aes(x = X3SEQ_prop_recombinant_sequences, y = proportion_internal_branches, color = X3SEQ_p_value)) + geom_point() + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_proportionInternalBranches_comparison.png") , plot = p)
   # look for a correlation between the tree proportion test and the proportion of internal branches
   p <- ggplot(data = op_df, aes(x = tree_proportion, y = proportion_informative_sites, color = tree_proportion_p_value)) + geom_point() + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_treeProportion_proportionInformativeSites_comparison.png") , plot = p)
   p <- ggplot(data = op_df, aes(x = X3SEQ_prop_recombinant_sequences, y = proportion_informative_sites, color = X3SEQ_p_value)) + geom_point() + theme_bw()
   ggsave(filename = paste0(plot_dirs[[dataset]],dataset,"_3seq_proportionInformativeSites_comparison.png") , plot = p)
}




##### Step 4: Plot histograms #####
# Create a facetted histogram to compare the distributions of different test statistics
# Plot histogram with all test statistics and p_values




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
# Repeat this process for p values


##### Step 6: Plot test statistics against potential predictors #####


##### Step 7: Explore relationships between test statistics #####


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


