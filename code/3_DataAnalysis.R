### empirical_treelikeness/3_DataAnalysis.R
## R program to plot and explore results of the treelikeness test statistics on empirical data
# Final result is a reformatted csv file and a number of graphs
# Caitlin Cherryh 2021



##### Step 1: Open packages #####
library(treespace)
library(phangorn) # using phangorn to get the KF.dist, RF.dist, wRF.dist, nNodes, and patristic methods for summarising trees as vectors 
# these methods all assume an unrooted tree so trees can be used as is for this analysis



##### Step 2: Set the file paths for input and output files, and necessary functions/directories #####
print("set filepaths")
# treedir           <- "treelikeness" repository location (github.com/caitlinch/treelikeness)
# maindir           <- "empirical_treelikeness" repository location (github.com/caitlinch/empirical_treelikeness)
# output_dir        <- for collated output and results from treelikeness analysis. This file should contain a folder for each input_name (where the folder name and corresponding input_name are identical)
# datasets          <- set name(s) for the dataset(s)
## Choose which sections of the script you want to run
# collect_warnings  <- Whether to collect warnings from the IQ-Tree log and iqtree files. Can be TRUE or FALSE.
# collate_results   <- Whether to collect results from the output dir. Can be TRUE of FALSE.
# run_analysis      <- Whether to perform data analysis. Can be TRUE or FALSE.

location = "local"
# location = "server"

if (location == "local"){
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the treelikeness code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/" # where the empirical treelikeness code is
  csv_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/03_output/"
  tree_data_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/04_trees/"
  output_dir <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/05_results/"
  
  datasets <- c("Vanderpool2020")
  dataset = "Vanderpool2020"
  
  collect_warnings  <- FALSE
  collate_results   <- FALSE
  run_analysis      <- TRUE
  
} else if (location == "server"){
  treedir <- "/data/caitlin/treelikeness/" # where the treelikeness repository/folder is
  maindir <- "/data/caitlin/empirical_treelikeness/" # where the empirical treelikeness repository/folder is 
  csv_data_dir <- "/data/caitlin/empirical_treelikeness/Output/"
  tree_data_dir <- "/data/caitlin/empirical_treelikeness/Output_treeEstimation/"
  output_dir <- "/data/caitlin/empirical_treelikeness/Output_DataAnalysis/"
  
  datasets <- c("Vanderpool2020")
  dataset = "Vanderpool2020"
  
  collect_warnings  <- FALSE
  collate_results   <- FALSE
  run_analysis      <- TRUE
}



#### Step 3: Open files
source(paste0(maindir,"code/func_analysis.R"))

# Construct the folder names for where the estimated species trees are
tree_data_dirs <- paste0(tree_data_dir,datasets,"/")
names(tree_data_dirs) <- datasets
for (f in tree_data_dirs){
  if (!file.exists(f)){
    dir.create(f)
  }
}

# Construct the folder names for where to save the data analysis results
output_dirs <- paste0(output_dir,datasets,"/")
names(output_dirs) <- datasets
for (f in output_dirs){
  if (!file.exists(f)){
    dir.create(f)
  }
}



##### Step 3: Collect any warnings from tree estimation ####
# Identify any warnings from the IQ-Tree loci tree estimation
if (collect_warnings == TRUE){
  for (dataset in datasets){
    # Open this dataset's raw output file from the treelikeness analysis 
    all_csv_files <- grep(".csv",list.files(csv_data_dir), value = TRUE)
    all_untrimmed_csv_files <- grep("trimmed",all_csv_files, value = TRUE, invert = TRUE)
    dataset_csv_file <- grep(dataset, all_untrimmed_csv_files, value = TRUE)
    dataset_df <- read.csv(paste0(csv_data_dir,dataset_csv_file), stringsAsFactors = FALSE)
    
    # Take list of alignments from the raw output file
    all_alignments <- dataset_df$alignment_file
    
    # Collect warnings and write out as a csv file
    warning_df <- as.data.frame(do.call(rbind, (lapply(all_alignments, check.for.IQTree.warnings))))
    warning_df_file <- paste0(csv_data_dir, "empiricalTreelikeness_", dataset, "_collated_IQ-Tree_warnings.csv")
    write.csv(warning_df, file = warning_df_file)
  }
  
  # Investigate warnings
  warnings_df <- read.csv(paste0(csv_data_dir, grep("warnings", all_csv_files, value = TRUE)))
  log_df <- warnings_df[warnings_df$file == ".log",]
  unique(log_df$loci)
  unique(log_df$warnings)
  
  iq_df <- warnings_df[warnings_df$file == ".iqtree",]
  length(unique(iq_df$loci))
}



##### Step 4: Collate results from tree estimation #####
if (collate_results == TRUE){
  for (dataset in datasets){
    # Get a list of all the files in the folder
    all_files <- list.files(tree_data_dirs[[dataset]])
    # Collect the names of the species trees from ASTRAL and IQ-Tree
    astral_tree_files <- grep("ASTRAL_species.tre", all_files, value = TRUE)
    partition_tree_files <- paste0(grep("IQ-Tree_partition", all_files, value = TRUE), "/partitions.nex.contree")
    all_trees <- paste0(tree_data_dirs[[dataset]],c(astral_tree_files, partition_tree_files))
    
    # Create a new folder to copy only the species trees into
    copy_folder <- paste0(tree_data_dirs[[dataset]],"all_species_trees/")
    copy_partition_tree_files <- gsub("/partitions","", partition_tree_files)
    copy_trees <- paste0(copy_folder,c(astral_tree_files, copy_partition_tree_files))
    # Create the copy folder if it doesn't already exist
    if (!dir.exists(copy_folder)){
      dir.create(copy_folder)
    }
    # Copy files into this new folder
    for (i in 1:length(all_trees)){
      tree_original_location <- all_trees[i]
      tree_copied_location <- copy_trees[i]
      file.copy(from = tree_original_location, to = tree_copied_location, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
    
    # Get a list of each loci involved in each analysis
    partition_folders <- paste0(tree_data_dirs[[dataset]], grep("IQ-Tree_partition", all_files, value = TRUE))
    lapply(partition_folders, get.loci.from.analysis, copy_folder)
  }
}



##### Step 5: Compare species trees #####
if (run_analysis == TRUE){
  for (dataset in datasets) {
    #### Comparing trees from p-value categories
    ## Want to see if loci trees within one category are more similar on average then loci trees from different categories.
    
  }
}


# 
# ### Classify trees into 4 groups: 3seq only, tp only, both, neither
# ## plot and compare treespace - do trees group together?
# ## plot and compare groves (change colour to be based on group!) - do trees group together?
# # Start by adding columns to the dataframe that classify each value as either TREELIKE or NON-TREELIKE
# t29_df <- classify.treelikeness.statistics(t29_df, 0.7)
# 
# # Now plot the treelikeness test statistics against each other and colour by group
# pretty_colours <- RColorBrewer::brewer.pal(5,"YlGnBu")[2:5]
# dark_colours <- RColorBrewer::brewer.pal(4,"Dark2")
# ggplot(data = t29_df, aes(x = tree_proportion_p_value, y = X3SEQ_p_value, color = sorted_p_value)) + geom_point() + theme_bw() + 
#   scale_color_manual(labels = c("Neither (n = 779)","3seq only  (n = 587)","Tree proportion only  (n = 121)","Both (n = 223)"), values = dark_colours) + 
#   guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test p-value results", subtitle = "Vanderpool 2020") +
#   scale_x_continuous(name = "Tree proportion p-value") + scale_y_continuous(name = "3seq p-value")
# # Now plot the tree proportion test statistics against the 3seq p-value and colour by group
# ggplot(data = t29_df, aes(x = tree_proportion, y = X3SEQ_p_value, color = sorted)) + geom_point() + theme_bw() + 
#   scale_color_manual(labels = c("Neither (n = 378)","3seq only  (n = 380)","Tree proportion only  (n = 542)","Both (n = 430)"), values = dark_colours) + 
#   guides(color = guide_legend(title = "Loci treelikeness")) + labs(title = "Grouping loci by treelikeness test statistic results", subtitle = "Vanderpool 2020") +
#   scale_x_continuous(name = "Tree proportion") + scale_y_continuous(name = "3seq p-value")
# 
# ### Examine treespace
# # Break down trees to those with all species
# t29_df <- t29_df[t29_df$n_taxa == 29,]
# #t29_df <- t29_df[1:20,] # small dataframe for testing - comment/uncomment as needed
# # Read in all trees at once
# mp <- read.tree(text = t29_df$tree)
# # use treespace <- have to pick a method that works on unrooted trees. Do this by picking one that assumes trees are unrooted
# # "Warning message: In is.euclid(distmat) : Zero distance(s)" may appear when doing "RF" method <- this means you have duplicate rows in your distance matrix
# # Basically means two trees have the exact same values so the function thinks that you have a duplicate row
# res <- treespace(mp, nf=3, method = "wRF")
# # Use ggplot to plot the treespace using colours to sort loci as above
# # plotGroves(res$pco, lab.show=FALSE) # plotting treespace using plotGroves function = simple but hard to add colours for factors
# # Extract principal component vectors and add the sorted column so you know how to group the taxa
# pc_df <- res$pco$tab
# pc_df$sorted <- t29_df$sorted
# pc_df$sorted_p_value <- t29_df$sorted_p_value
# # sorted_p_value plot
# ggplot(pc_df,aes(x = A1, y = A2, color = sorted_p_value)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Loci treelikeness")) +
#   scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace: A2 ~ A1", subtitle = "Grouped by 3seq and tree proportion p-values") + 
#   scale_color_brewer(palette = "Dark2")
# 
# 
# # Now find the groves
# groves <- findGroves(res, nclust=5)
# # Plot using plotGroves to see how the groves are organised
# plotGroves(groves)
# # Create a data frame to use for plotting
# groves_pc_df <- groves$treespace$pco$li
# groves_pc_df$groups <- groves$groups
# # plot
# ggplot(groves_pc_df,aes(x = A1, y = A2, color = groups)) + geom_point() + theme_bw() + guides(color = guide_legend(title = "Groves")) +
#   scale_x_continuous(name = "A1") + scale_y_continuous(name = "A2") + labs(title = "Vanderpool 2020 treespace with groves: A2 ~ A1") + 
#   scale_color_brewer(palette = "Paired")
# 
# # Plot and save a 3D plot
# colours <- fac2col(groves$groups, col.pal=colorRampPalette(RColorBrewer::brewer.pal(5,"Paired")))
# plot3d(groves$treespace$pco$li[,1],
#        groves$treespace$pco$li[,2],
#        groves$treespace$pco$li[,3],
#        col=colours, type="s", size=1.5,
#        xlab="", ylab="", zlab="")
# 



