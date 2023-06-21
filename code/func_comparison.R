### gene_filtering/code/func_comparison.R
## R functions for comparing two species trees
# Caitlin Cherryh 2023

library(phytools)
library(ape)

# Function to read in tree from ASTRAL and edit it to use in the QuartetNetworkGoodnessOfFit Julia package
# ASTRAL does not output terminal branch lengths, so these are arbitrarily given the length 1
# ASTRAL includes posterior probabilities: these are stripped (the Julia function cannot handle extra node values)
reformat.ASTRAL.tree.for.Julia <- function(tree_path, add.arbitrary.terminal.branches = FALSE, terminal.branch.length = 1){
  # Open tree
  t <- read.tree(tree_path)
  if (add.arbitrary.terminal.branches == TRUE){
    # Identify which edges correspond to terminal branches
    n = length(t$tip.label)
    missing_branch_indexes <- sapply(1:n,function(x,y)   which(y==x),y=t$edge[,2])
    # Set length of terminal branches to 1 (or to whatever the terminal.branch.length variable is)
    t$edge.length[missing_branch_indexes] <- terminal.branch.length
  }
  # Don't want node labels - for quartnetGoFtest, need just taxa and branch lengths
  t$node.label <- NULL
  # Write this tree back out
  write.tree(t, tree_path)
}



# Function to read in list of genetrees from ASTRAL and edit to use in the QuartetNetworkGoodnessOfFit Julia package
# ASTRAL does not output terminal branch lengths, so these are arbitrarily given the length 1
# ASTRAL includes posterior probabilities: these are stripped (the Julia function cannot handle extra node values)
reformat.gene.tree.list.for.Julia <- function(trees_path, add.arbitrary.terminal.branches = FALSE, terminal.branch.length = 1){
  # Open trees
  ts <- read.tree(trees_path)
  for (i in 1:length(ts)){
    temp_tree <- ts[[i]]
    temp_tree <- reformat.phylo.for.Julia(temp_tree, add.arbitrary.terminal.branches, terminal.branch.length)
    ts[[i]] <- temp_tree
  }
  # Write these trees back out
  write.tree(ts, trees_path)
}

# Function to take one tree from a multiphylo object, format it for QuartetNetworkGoodnessOfFit Julia package, and return it
reformat.phylo.for.Julia <- function(tree, add.arbitrary.terminal.branches = FALSE, terminal.branch.length = 1){
  # If gene trees from ASTRAL, the terminal branches are not output and need to be set to one
  if (add.arbitrary.terminal.branches == TRUE){
    # Identify which edges correspond to terminal branches
    n = length(tree$tip.label)
    missing_branch_indexes <- sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])
    # Set length of terminal branches to 1 (or to whatever the terminal.branch.length variable is)
    tree$edge.length[missing_branch_indexes] <- terminal.branch.length
  }
  # Gene trees from ASTRAL have posterior probabilities, and those from IQ-Tree have bootstraps.
  # Either way, don't want node labels for quarnetGoFtest - need just taxa and branch lengths
  tree$node.label <- NULL
  return(tree)
}




# Function to take one tree from a file path, read it in, and save it as un ultrametric tree
make.tree.ultrametric <- function(tree_path, root.tree = FALSE, outgroup = NA){
  # Read in the tree
  tree <- read.tree(file = tree_path)
  # If root.tree == FALSE, then reroot the tree at the provided outgroup
  if (root.tree == TRUE){
    # Select the outgroup species which are within the tree
    present_outgroup_species <- outgroup[which(outgroup %in% tree$tip.label)]
    # Root the tree at the desired outgroup
    tree <- root(tree, present_outgroup_species)
  }
  # Extend the terminal branches to make the tree ultrametric
  ultrametric_tree <- force.ultrametric(tree, method = "extend")
  # Write out the ultrametric tree to the same filepath
  write.tree(ultrametric_tree, file = tree_path, append = FALSE)
}




# This function takes in the locations of multiple files and writes a Julia script to apply the 
# quarnetGoFtest from the QuartetNetworkGoodnessOfFit Julia package
write.Julia.GoF.script <- function(test_name, dataset, directory, pass_tree, fail_tree, all_tree, gene_trees, root.species.trees = FALSE, tree_root = NA, 
                                   output_csv_file_path, number_of_simulated_replicates = 1000){
  # Make name of output files
  script_file <- paste0(directory, "apply_GoF_test.jl")
  gene_cf_file <- gsub(".txt", "_geneCF.txt", gene_trees)
  op_df_file <- output_csv_file_path
  # Select a random seed
  seed <- round(as.numeric(Sys.time()))
  # Add code lines 
  lines <- c("# Code to take in a list of species trees, estimate concordance factors, and compare them to three trees",
             "# One tree estimated with all the data, one tree estimated with the 'good' data and one with the 'bad' data",
             "# Open packages",
             "using PhyloNetworks",
             "using QuartetNetworkGoodnessFit",
             "using DataFrames, CSV",
             "",
             "# Set working directory")
  # Add working directory
  lines <- c(lines, paste0('cd("', directory, '")'), '')
  # Add code for converting trees to quartet concordance factors
  if (file.exists(gene_cf_file) == TRUE){
    # If the gene trees have already been converted into concordance factors, open the concordance factors csv file
    lines <- c(lines,
               '# Open quartet concordance factors file',
               paste0('genetrees_cf = readTableCF("', gene_cf_file, '")'), 
               '')
  } else if (file.exists(gene_cf_file) == FALSE){
    # If the quartet concordance factors have not been calculated, calculate them from the gene trees
    lines <- c(lines,
               '# Convert list of gene trees to concordance factors - using trees with bootstraps works fine',
               paste0('genetrees_cf = readTrees2CF("', gene_trees,
                      '", writeTab=true, CFfile="', gene_cf_file, '")'), 
               '')
  }
  # Open trees
  lines <- c(lines, 
             '# Read in the three trees',
             paste0('tree_test_pass = readTopology("', pass_tree, '");'),
             paste0('tree_test_fail = readTopology("', fail_tree, '");'),
             paste0('tree_all = readTopology("', all_tree, '");'),
             '')
  # Optional: root species trees by provided outgroup
  if (root.species.trees == TRUE){
    lines <- c(lines,
               '# Root the three trees at the same taxa',
               paste0('PhyloNetworks.rootatnode!(tree_test_pass, "', tree_root, '");'),
               paste0('PhyloNetworks.rootatnode!(tree_test_fail, "', tree_root, '");'),
               paste0('PhyloNetworks.rootatnode!(tree_all, "', tree_root, '");'),
               '') 
  }
  # Apply the QuartetNetwork Goodness of Fit test
  lines <- c(lines,
             "# Apply the QuartetNetworkGoodnessFit test",
             paste0("gof_test_pass = quarnetGoFtest!(tree_test_pass, genetrees_cf, false; quartetstat=:LRT, correction=:simulation, seed=", seed ,
                    ", nsim=", number_of_simulated_replicates , ", verbose=false, keepfiles=false)"),
             paste0("gof_test_fail = quarnetGoFtest!(tree_test_fail, genetrees_cf, false; quartetstat=:LRT, correction=:simulation, seed=", seed,
                    ", nsim=", number_of_simulated_replicates, ", verbose=false, keepfiles=false)"),
             paste0("gof_test_all = quarnetGoFtest!(tree_all, genetrees_cf, false; quartetstat=:LRT, correction=:simulation, seed=", seed ,
                    ", nsim=", number_of_simulated_replicates, ", verbose=false, keepfiles=false)"),
             "")
  # Create an output dataframe
  if (length(tree_root) > 1){
    tree_root_op <- paste(tree_root, collapse = ",")
  } else {
    tree_root_op <- tree_root
  }
  lines <- c(lines,
             '# Create a dataframe using the variables from the three gof tests',
             paste0('df = DataFrame(dataset = ["', dataset, '", "', dataset, '", "', dataset, '"],'),
             paste0('               concordance_factors = ["', basename(gene_cf_file), '", "', basename(gene_cf_file), '", "', basename(gene_cf_file), '"],'),
             paste0('               test = ["', test_name, '", "', test_name, '", "', test_name, '"],'),
             '               tree = ["test_pass", "test_fail", "no_test"],',
             paste0('               outgroup = ["', tree_root_op, '", "', tree_root_op, '", "', tree_root_op, '"],'),
             '               p_value_overall_GoF_test = [gof_test_pass[1], gof_test_fail[1], gof_test_all[1]],',
             '               uncorrected_z_value_test_statistic = [gof_test_pass[2], gof_test_fail[2], gof_test_all[2]],',
             '               estimated_sigma_for_test_statistic_correction = [gof_test_pass[3], gof_test_fail[3], gof_test_all[3]]',
             ')',
             '')
  # Save the output dataframe
  lines <- c(lines,
             '# Write output dataframe as .csv',
             paste0('CSV.write("', op_df_file, '",df)'))
  # Output script
  write(lines, file = script_file)
}





# This function takes in the locations of multiple files and writes a Julia script to apply the 
# quarnetGoFtest from the QuartetNetworkGoodnessOfFit Julia package
write.Julia.GoF.script.two.trees <- function(test_name, dataset, directory, pass_tree, all_tree, gene_trees, root.species.trees = FALSE, tree_root = NA, 
                                   output_csv_file_path, number_of_simulated_replicates = 1000){
  # Make name of output files
  script_file <- paste0(directory, "apply_GoF_test.jl")
  gene_cf_file <- gsub(".txt", "_geneCF.txt", gene_trees)
  op_df_file <- output_csv_file_path
  # Select a random seed
  seed <- round(as.numeric(Sys.time()))
  # Add code lines 
  lines <- c("# Code to take in a list of species trees, estimate concordance factors, and compare them to two trees",
             "# One tree estimated with all the data and one tree estimated with the 'good' data",
             "# Open packages",
             "using PhyloNetworks",
             "using PhyloPlots",
             "using QuartetNetworkGoodnessFit",
             "using DataFrames, CSV",
             "",
             "# Set working directory")
  # Add working directory
  lines <- c(lines, paste0('cd("', directory, '")'), '')
  # Add code for converting trees to quartet concordance factors
  if (file.exists(gene_cf_file) == TRUE){
    # If the gene trees have already been converted into concordance factors, open the concordance factors csv file
    lines <- c(lines,
               '# Open quartet concordance factors file',
               paste0('genetrees_cf = readTableCF("', gene_cf_file, '")'), 
               '')
  } else if (file.exists(gene_cf_file) == FALSE){
    # If the quartet concordance factors have not been calculated, calculate them from the gene trees
    lines <- c(lines,
               '# Convert list of gene trees to concordance factors - using trees with bootstraps works fine',
               paste0('genetrees_cf = readTrees2CF("', gene_trees,
                      '", writeTab=true, CFfile="', gene_cf_file, '")'), 
               '')
  }
  # Open trees
  lines <- c(lines, 
             '# Read in the two trees',
             paste0('tree_test_pass = readTopology("', pass_tree, '");'),
             paste0('tree_all = readTopology("', all_tree, '");'),
             '')
  # Optional: root species trees by provided outgroup
  if (root.species.trees == TRUE){
    lines <- c(lines,
               '# Root the two trees at the same taxa',
               paste0('PhyloNetworks.rootatnode!(tree_test_pass, "', tree_root, '");'),
               paste0('PhyloNetworks.rootatnode!(tree_all, "', tree_root, '");'),
               '') 
  }
  # Apply the QuartetNetwork Goodness of Fit test
  lines <- c(lines,
             "# Apply the QuartetNetworkGoodnessFit test",
             paste0("gof_test_pass = quarnetGoFtest!(tree_test_pass, genetrees_cf, false; quartetstat=:LRT, correction=:simulation, seed=", seed ,
                    ", nsim=", number_of_simulated_replicates , ", verbose=false, keepfiles=false)"),
             paste0("gof_test_all = quarnetGoFtest!(tree_all, genetrees_cf, false; quartetstat=:LRT, correction=:simulation, seed=", seed ,
                    ", nsim=", number_of_simulated_replicates, ", verbose=false, keepfiles=false)"),
             "")
  # Create an output dataframe
  if (length(tree_root) > 1){
    tree_root_op <- paste(tree_root, collapse = ",")
  } else {
    tree_root_op <- tree_root
  }
  lines <- c(lines,
             '# Create a dataframe using the variables from the three gof tests',
             paste0('df = DataFrame(dataset = ["', dataset, '", "', dataset, '"],'),
             paste0('               concordance_factors = ["', basename(gene_cf_file), '", "', basename(gene_cf_file), '"],'),
             paste0('               test = ["', test_name,  '", "', test_name, '"],'),
             '               tree = ["test_pass", "no_test"],',
             paste0('               outgroup = ["', tree_root_op, '", "', tree_root_op, '"],'),
             '               p_value_overall_GoF_test = [gof_test_pass[1], gof_test_all[1]],',
             '               uncorrected_z_value_test_statistic = [gof_test_pass[2], gof_test_all[2]],',
             '               estimated_sigma_for_test_statistic_correction = [gof_test_pass[3], gof_test_all[3]]',
             ')',
             '')
  # Save the output dataframe
  lines <- c(lines,
             '# Write output dataframe as .csv',
             paste0('CSV.write("', op_df_file, '",df)'))
  # Output script
  write(lines, file = script_file)
}


