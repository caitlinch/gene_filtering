### empirical_treelikeness/2_ExploreTreelikenessScores.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
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


