### gene_filtering/code/func_comparison.R
## R functions for comparing two or more species trees
# Caitlin Cherryh 2024

library(phytools)
library(ape) # read.tree, Ntip, root
library(phangorn) # as.splits
library(dplyr) # bind_rows


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



#### Functions for comparing splits within trees ####
compare.splits.wrapper <- function(i, df){
  # Iterate through rows and return results of compare.splits.2.trees function
  
  temp_row  <- df[i, ]
  temp_results  <- compare.splits.2.trees(clean_tree_path = paste0(temp_row$tree_directory, temp_row$clean_tree),
                                          comparison_tree_path = paste0(temp_row$tree_directory, temp_row$comparison_tree))
  temp_results$dataset <- temp_row$dataset
  temp_results$clean_tree <- temp_row$clean_tree
  temp_results$comparison_tree <- temp_row$comparison_tree
  temp_results$clean_id <- temp_row$clean_id
  temp_results$comparison_id <- temp_row$comparison_id
  temp_results$recombination_test <- temp_row$recombination_test
  temp_results$comparison_gene_status <- temp_row$comparison_gene_status
  temp_results$tree <- factor(temp_results$tree,
                              levels = c("Clean", "Comp", "whole_clean_tree", "whole_comp_tree"),
                              labels = c(temp_row$clean_id, temp_row$comparison_id, "Clean_allTips", "Comp_allTips") )
  
  return(temp_results)
}



compare.splits.2.trees <- function(clean_tree_path, comparison_tree_path){
  ## Return qCF values for all splits within two trees, separated into conflicting and concordant splits
  
  ## Read in trees
  # Open both trees as phylo objects
  clean_tree  <- read.tree(clean_tree_path)
  comp_tree   <- read.tree(comparison_tree_path)
  
  ## Check whether the two trees have the same number of tips
  if (Ntip(clean_tree) == Ntip(comp_tree)){
    ## Trees have the same number of splits
    # Convert trees to splits
    clean_splits  <- as.splits(clean_tree)
    comp_splits   <- as.splits(comp_tree)
    # Remove trivial splits
    clean_splits  <- removeTrivialSplits(clean_splits)
    comp_splits   <- removeTrivialSplits(comp_splits)
    
    ## Create formatted dataframes of the splits for each tree
    # Convert to dataframes
    clean_df  <- as.data.frame(as.matrix(clean_splits))
    comp_df   <- as.data.frame(as.matrix(comp_splits))
    
    ## Set "whole tree" values - will not be needed as trees have identical tips
    whole_df <- NA
    whole_df_type <- "logical"
  } else {
    ## Trees have different numbers of splits
    ## We want to do 2 things: firstly extract all splits from the clean tree ("whole tree"), 
    #     and secondly keeping only taxa in the comp tree, compare splits in both trees ("clean tree" and "comp tree")
    
    ## Create formatted dataframes of the splits for each tree
    # Identify tips present in both trees
    both_tree_tips <- intersect(clean_tree$tip.label, comp_tree$tip.label)
    # Remove tips in the clean_tree that are not present in the comp_tree
    clean_tree    <- keep.tip(clean_tree, both_tree_tips)
    comp_tree     <- keep.tip(comp_tree, both_tree_tips)
    # Convert to splits
    clean_splits  <- as.splits(clean_tree)
    comp_splits   <- as.splits(comp_tree)
    # Remove trivial splits
    clean_splits  <- removeTrivialSplits(clean_splits)
    comp_splits   <- removeTrivialSplits(comp_splits)
    # Convert to dataframes
    clean_df      <- as.data.frame(as.matrix(clean_splits))
    comp_df       <- as.data.frame(as.matrix(comp_splits))
    
    ## Extract non-trivial splits for the whole trees
    whole_clean_tree <- read.tree(clean_tree_path)
    whole_comp_tree <- read.tree(comparison_tree_path)
    
    ## Extract splits present only in the whole clean tree
    if (Ntip(clean_tree) == Ntip(whole_clean_tree) & setequal(clean_tree$tip.label, whole_clean_tree$tip.label)){
      whole_clean_df      <- NA
      whole_clean_df_type <- "logical"
    } else {
      ## Nicely format the splits present only in the whole clean tree
      # Extract splits
      whole_clean_tree_splits   <- removeTrivialSplits(as.splits(whole_clean_tree))
      # Convert splits to dataframe
      whole_clean_df            <- as.data.frame(as.matrix(whole_clean_tree_splits))
      # Add details from the splits
      whole_clean_df$tree       <- "whole_clean_tree"
      whole_clean_df$confidence <- attr(whole_clean_tree_splits, "confidence")
      whole_clean_df$weights    <- attr(whole_clean_tree_splits, "weight")
      # Remove any rows with either "" or NA as the confidence
      remove_whole_rows         <- which( ("" == whole_clean_df$confidence) | is.na(whole_clean_df$confidence) )
      if (identical(remove_whole_rows, integer(0)) == FALSE){
        keep_whole_rows <- setdiff(1:nrow(whole_clean_df), remove_whole_rows)
        whole_clean_df  <- whole_clean_df[keep_whole_rows, ]
      } 
      # Extract qCF values
      whole_confidence_splits  <- gsub("\\[", "", gsub("\\]", "", gsub("'", "", whole_clean_df$confidence)))
      whole_confidence_vals    <- strsplit(whole_confidence_splits, ";")
      whole_confidence_vals    <- lapply(1:length(whole_confidence_vals), function(i){unlist(strsplit(whole_confidence_vals[[i]], "="))[c(FALSE,TRUE)]})
      whole_qcf_df             <- as.data.frame(do.call(rbind, whole_confidence_vals))
      names(whole_qcf_df)      <- c("q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")
      # Add extra columns
      whole_qcf_df$split_type <- "Clean_allTips"
      # Bind dataframes
      whole_clean_df      <- cbind(whole_clean_df, whole_qcf_df)
      # Rearrange whole_clean_df
      whole_clean_df      <- whole_clean_df[ , c("tree", "split_type", "weights", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")] 
      whole_clean_df_type <- "dataframe"
    }
    
    ## Extract splits present only in the whole comparison tree
    if (Ntip(comp_tree) == Ntip(whole_comp_tree) & setequal(comp_tree$tip.label, whole_comp_tree$tip.label)){
      whole_comp_df      <- NA
      whole_comp_df_type <- "logical"
    } else {
      ## Nicely format the splits present only in the whole comparison tree
      # Extract splits
      whole_comp_tree_splits  <- removeTrivialSplits(as.splits(whole_comp_tree))
      # Convert splits to dataframe
      whole_comp_df           <- as.data.frame(as.matrix(whole_comp_tree_splits))
      # Add details from the splits
      whole_comp_df$tree       <- "whole_comp_tree"
      whole_comp_df$confidence <- attr(whole_comp_tree_splits, "confidence")
      whole_comp_df$weights    <- attr(whole_comp_tree_splits, "weight")
      # Remove any rows with either "" or NA as the confidence
      remove_wcomp_rows <- which( ("" == whole_comp_df$confidence) | is.na(whole_comp_df$confidence) )
      if (identical(remove_wcomp_rows, integer(0)) == FALSE){
        keep_wcomp_rows <- setdiff(1:nrow(whole_comp_df), remove_wcomp_rows)
        whole_comp_df   <- whole_comp_df[keep_wcomp_rows, ]
      } 
      # Extract qCF values
      wcomp_confidence_splits  <- gsub("\\[", "", gsub("\\]", "", gsub("'", "", whole_comp_df$confidence)))
      wcomp_confidence_vals    <- strsplit(wcomp_confidence_splits, ";")
      wcomp_confidence_vals    <- lapply(1:length(wcomp_confidence_vals), function(i){unlist(strsplit(wcomp_confidence_vals[[i]], "="))[c(FALSE,TRUE)]})
      wcomp_qcf_df             <- as.data.frame(do.call(rbind, wcomp_confidence_vals))
      names(wcomp_qcf_df)      <- c("q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")
      # Add extra columns
      wcomp_qcf_df$split_type  <- "Comp_allTips"
      # Bind dataframes
      whole_comp_df            <- cbind(whole_comp_df, wcomp_qcf_df)
      # Rearrange whole_comp_df
      whole_comp_df            <- whole_comp_df[ , c("tree", "split_type", "weights", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")]
      whole_comp_df_type       <- "dataframe"
    }
    
    ## Collate the two datasets from the whole trees
    if (whole_clean_df_type == "dataframe" & whole_comp_df_type == "dataframe"){
      # Both trees had tips removed
      whole_df <- rbind(whole_clean_df, whole_comp_df)
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "dataframe" & whole_comp_df_type == "logical"){
      # Only the "clean" tree had tips removed
      whole_df <- whole_clean_df
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "logical" & whole_comp_df_type == "dataframe"){
      # Only the "comparison" tree had tips removed
      whole_df <- whole_comp_df
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "logical" & whole_comp_df_type == "logical"){
      # Neither tree had tips removed
      whole_df <- NA
      whole_df_type <- "logical"
    } 
    
  }
  
  ## Format dataframes
  # Reorder columns
  col_order <- colnames(clean_df)
  comp_df <- comp_df[, col_order]
  # Add details to each dataframe
  clean_df$tree       <- "Clean"
  clean_df$confidence <- attr(clean_splits, "confidence")
  clean_df$weights    <- attr(clean_splits, "weight")
  comp_df$tree        <- "Comp"
  comp_df$confidence  <- attr(comp_splits, "confidence")
  comp_df$weights     <- attr(comp_splits, "weight")
  
  ## Identify the two types of splits: splits that are present in both trees, and splits that are present in one tree only
  all_df <- rbind(clean_df, comp_df)
  all_tree_splits <- all_df |>
    group_by_at(vars(-tree,-confidence, -weights)) %>% 
    filter(n() == 2) |>
    ungroup()
  unique_splits <- all_df |>
    group_by_at(vars(-tree,-confidence, -weights)) %>% 
    filter(n() < 2) |>
    ungroup()
  
  ## Extract qCF values for the splits that occur in both trees
  if (nrow(all_tree_splits) > 0){
    # Remove rows without qCF values
    shared_qcf_inds <- which(is.na(all_tree_splits$confidence) == FALSE & identical("", all_tree_splits$confidence) == FALSE)
    all_tree_splits <- all_tree_splits[shared_qcf_inds, ]
    # Create small dataframe for confidence intervals
    shared_confidence_splits  <- gsub("\\[", "", gsub("\\]", "", gsub("'", "", all_tree_splits$confidence)))
    shared_confidence_vals    <- strsplit(shared_confidence_splits, ";")
    shared_confidence_vals    <- lapply(1:length(shared_confidence_vals), function(i){unlist(strsplit(shared_confidence_vals[[i]], "="))[c(FALSE,TRUE)]})
    shared_df                 <- as.data.frame(do.call(rbind, shared_confidence_vals))
    names(shared_df)          <- c("q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")
    # Bind dataframes together
    all_tree_splits <- cbind(all_tree_splits, shared_df)
    all_tree_splits_type <- "dataframe"
    # Add new column for the type of split
    all_tree_splits$split_type      <- "Concordant"
  } else {
    all_tree_splits <- NA
    all_tree_splits_type <- "logical"
  }
  
  ## Extract qCF values for the splits that occur in only one tree
  if (nrow(unique_splits) > 0){
    # Remove rows without qCF values
    unique_qcf_inds <- which(is.na(unique_splits$confidence) == FALSE & identical("", unique_splits$confidence) == FALSE)
    unique_splits   <- unique_splits[unique_qcf_inds, ]
    # Create small dataframe for confidence intervals
    unique_confidence_splits  <- gsub("\\[", "", gsub("\\]", "", gsub("'", "", unique_splits$confidence)))
    unique_confidence_vals    <- strsplit(unique_confidence_splits, ";")
    unique_confidence_vals    <- lapply(1:length(unique_confidence_vals), function(i){unlist(strsplit(unique_confidence_vals[[i]], "="))[c(FALSE,TRUE)]})
    unique_df                 <- as.data.frame(do.call(rbind, unique_confidence_vals))
    names(unique_df)          <- c("q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")
    # Bind dataframes together
    unique_splits   <- cbind(unique_splits, unique_df)
    unique_splits_type <- "dataframe"
    # Add new column for the type of split
    unique_splits$split_type      <- "Conflicting"
  } else {
    unique_splits <- NA
    unique_splits_type <- "logical"
  }
  
  ## Bind dataframes of split types together, for the dataframes of split types that exist
  if (all_tree_splits_type == "dataframe" & unique_splits_type == "dataframe"){
    results_df   <- rbind(all_tree_splits, unique_splits)
  } else if (all_tree_splits_type == "dataframe" & unique_splits_type == "logical"){
    results_df  <- all_tree_splits
  } else if (all_tree_splits_type == "logical" & unique_splits_type == "dataframe"){
    results_df  <- unique_splits
  } else if (all_tree_splits_type == "logical" & unique_splits_type == "logical"){
    results_df <- as.data.frame(matrix(ncol = 14, nrow = 0))
    names(results_df) <- c("tree", "split_type", "weights", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")
  }
  # Trim columns
  results_df <- results_df[ , c("tree", "split_type", "weights", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN")]
  
  ## Add any splits from the whole tree, if the two trees have different numbers of tips
  if (whole_df_type == "dataframe"){
    results_df <- rbind(results_df, whole_df)
  }
  
  # Return the number of different splits
  return(results_df)
}


compare.branch.wrapper <- function(i, df){
  # Iterate through rows and return results of compare.branches.2.trees function
  
  temp_row  <- df[i, ]
  temp_results  <- compare.branches.2.trees(clean_tree_path = paste0(temp_row$tree_directory, temp_row$clean_tree),
                                            comparison_tree_path = paste0(temp_row$tree_directory, temp_row$comparison_tree),
                                            dataset = temp_row$dataset,
                                            analysis_method = temp_row$analysis_method)
  temp_results$dataset <- temp_row$dataset
  temp_results$clean_tree <- temp_row$clean_tree
  temp_results$comparison_tree <- temp_row$comparison_tree
  temp_results$clean_id <- temp_row$clean_id
  temp_results$comparison_id <- temp_row$comparison_id
  temp_results$recombination_test <- temp_row$recombination_test
  temp_results$comparison_gene_status <- temp_row$comparison_gene_status
  temp_results$tree <- factor(temp_results$tree,
                              levels = c("Clean", "Comp", "whole_clean_tree", "whole_comp_tree"),
                              labels = c(temp_row$clean_id, temp_row$comparison_id, "Clean_allTips", "Comp_allTips") )
  temp_results$analysis_method <- temp_row$analysis_method
  
  return(temp_results)
}


compare.branches.2.trees <- function(clean_tree_path, comparison_tree_path, dataset, analysis_method){
  ## Return qCF values for all splits within two trees, separated into conflicting and concordant splits
  
  ## Read in trees
  # Open both trees as phylo objects
  clean_tree  <- read.tree(clean_tree_path)
  comp_tree   <- read.tree(comparison_tree_path)
  
  ## Check whether the two trees have the same number of tips
  if (Ntip(clean_tree) == Ntip(comp_tree)){
    ## Trees have the same number of splits
    # Convert trees to splits
    clean_splits  <- as.splits(clean_tree)
    comp_splits   <- as.splits(comp_tree)
    # Remove trivial splits
    clean_splits  <- removeTrivialSplits(clean_splits)
    comp_splits   <- removeTrivialSplits(comp_splits)
    
    ## Create formatted dataframes of the splits for each tree
    # Convert to dataframes
    clean_df  <- as.data.frame(as.matrix(clean_splits))
    comp_df   <- as.data.frame(as.matrix(comp_splits))
    
    ## Set "whole tree" values - will not be needed as trees have identical tips
    whole_df <- NA
    whole_df_type <- "logical"
  } else {
    ## Trees have different numbers of splits
    ## We want to do 2 things: firstly extract all splits from the clean tree ("whole tree"), 
    #     and secondly keeping only taxa in the comp tree, compare splits in both trees ("clean tree" and "comp tree")
    
    ## Create formatted dataframes of the splits for each tree
    # Identify tips present in both trees
    both_tree_tips <- intersect(clean_tree$tip.label, comp_tree$tip.label)
    # Remove tips in the clean_tree that are not present in the comp_tree
    clean_tree    <- keep.tip(clean_tree, both_tree_tips)
    comp_tree     <- keep.tip(comp_tree, both_tree_tips)
    # Convert to splits
    clean_splits  <- as.splits(clean_tree)
    comp_splits   <- as.splits(comp_tree)
    # Remove trivial splits
    clean_splits  <- removeTrivialSplits(clean_splits)
    comp_splits   <- removeTrivialSplits(comp_splits)
    # Convert to dataframes
    clean_df      <- as.data.frame(as.matrix(clean_splits))
    comp_df       <- as.data.frame(as.matrix(comp_splits))
    
    ## Extract non-trivial splits for the whole trees
    whole_clean_tree <- read.tree(clean_tree_path)
    whole_comp_tree <- read.tree(comparison_tree_path)
    
    ## Extract splits present only in the whole clean tree
    if (Ntip(clean_tree) == Ntip(whole_clean_tree) & setequal(clean_tree$tip.label, whole_clean_tree$tip.label)){
      whole_clean_df      <- NA
      whole_clean_df_type <- "logical"
    } else {
      ## Nicely format the splits present only in the whole clean tree
      # Extract splits
      whole_clean_tree_splits   <- removeTrivialSplits(as.splits(whole_clean_tree))
      # Convert splits to dataframe
      whole_clean_df            <- as.data.frame(as.matrix(whole_clean_tree_splits))
      # Add details from the splits
      whole_clean_df$tree       <- "whole_clean_tree"
      if (dataset == "Plants" & analysis_method == "IQTREE"){
        whole_clean_df$confidence <- NA
      } else {
        whole_clean_df$confidence <- attr(whole_clean_tree_splits, "confidence")
      }
      whole_clean_df$weights    <- attr(whole_clean_tree_splits, "weight")
      # Add extra columns
      whole_clean_df$split_type  <- "Clean_allTips"
      # Rearrange whole_clean_df
      whole_clean_df      <- whole_clean_df[ , c("tree", "split_type", "confidence", "weights")] 
      whole_clean_df_type <- "dataframe"
    }
    
    ## Extract splits present only in the whole comparison tree
    if (Ntip(comp_tree) == Ntip(whole_comp_tree) & setequal(comp_tree$tip.label, whole_comp_tree$tip.label)){
      whole_comp_df      <- NA
      whole_comp_df_type <- "logical"
    } else {
      ## Nicely format the splits present only in the whole comparison tree
      # Extract splits
      whole_comp_tree_splits  <- removeTrivialSplits(as.splits(whole_comp_tree))
      # Convert splits to dataframe
      whole_comp_df           <- as.data.frame(as.matrix(whole_comp_tree_splits))
      # Add details from the splits
      whole_comp_df$tree       <- "whole_comp_tree"
      if (dataset == "Plants" & analysis_method == "IQTREE"){
        whole_comp_df$confidence <- NA
      } else {
        whole_comp_df$confidence <- attr(whole_comp_tree_splits, "confidence")
      }
      whole_comp_df$weights    <- attr(whole_comp_tree_splits, "weight")
      # Add extra columns
      whole_comp_df$split_type  <- "Comp_allTips"
      # Rearrange whole_comp_df
      whole_comp_df            <- whole_comp_df[ , c("tree", "split_type", "confidence", "weights")]
      whole_comp_df_type       <- "dataframe"
    }
    
    ## Collate the two datasets from the whole trees
    if (whole_clean_df_type == "dataframe" & whole_comp_df_type == "dataframe"){
      # Both trees had tips removed
      whole_df <- rbind(whole_clean_df, whole_comp_df)
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "dataframe" & whole_comp_df_type == "logical"){
      # Only the "clean" tree had tips removed
      whole_df <- whole_clean_df
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "logical" & whole_comp_df_type == "dataframe"){
      # Only the "comparison" tree had tips removed
      whole_df <- whole_comp_df
      whole_df_type <- "dataframe"
    } else if (whole_clean_df_type == "logical" & whole_comp_df_type == "logical"){
      # Neither tree had tips removed
      whole_df <- NA
      whole_df_type <- "logical"
    } 
    
  }
  
  ## Format dataframes
  # Reorder columns
  col_order <- colnames(clean_df)
  comp_df <- comp_df[, col_order]
  # Add details to each dataframe
  clean_df$tree       <- "Clean"
  if (dataset == "Plants" & analysis_method == "IQTREE"){
    clean_df$confidence <- NA
  } else {
    clean_df$confidence <- attr(clean_splits, "confidence")
  }
  clean_df$weights    <- attr(clean_splits, "weight")
  comp_df$tree        <- "Comp"
  if (dataset == "Plants" & analysis_method == "IQTREE"){
    comp_df$confidence <- NA
  } else {
    comp_df$confidence <- attr(comp_splits, "confidence")
  }
  comp_df$weights     <- attr(comp_splits, "weight")
  
  ## Identify the two types of splits: splits that are present in both trees, and splits that are present in one tree only
  all_df <- rbind(clean_df, comp_df)
  all_tree_splits <- all_df |>
    group_by_at(vars(-tree,-confidence, -weights)) %>% 
    filter(n() == 2) |>
    ungroup()
  unique_splits <- all_df |>
    group_by_at(vars(-tree,-confidence, -weights)) %>% 
    filter(n() < 2) |>
    ungroup()
  
  ## Extract qCF values for the splits that occur in both trees
  if (nrow(all_tree_splits) > 0){
    # Add new column for the type of split
    all_tree_splits$split_type <- "Concordant"
    all_tree_splits_type = "dataframe"
  } else {
    all_tree_splits <- NA
    all_tree_splits_type <- "logical"
  }
  
  ## Extract qCF values for the splits that occur in only one tree
  if (nrow(unique_splits) > 0){
    # Add new column for the type of split
    unique_splits$split_type      <- "Conflicting"
    unique_splits_type = "dataframe"
  } else {
    unique_splits <- NA
    unique_splits_type <- "logical"
  }
  
  ## Bind dataframes of split types together, for the dataframes of split types that exist
  if (all_tree_splits_type == "dataframe" & unique_splits_type == "dataframe"){
    results_df   <- rbind(all_tree_splits, unique_splits)
  } else if (all_tree_splits_type == "dataframe" & unique_splits_type == "logical"){
    results_df  <- all_tree_splits
  } else if (all_tree_splits_type == "logical" & unique_splits_type == "dataframe"){
    results_df  <- unique_splits
  } else if (all_tree_splits_type == "logical" & unique_splits_type == "logical"){
    results_df <- as.data.frame(matrix(ncol = 14, nrow = 0))
    names(results_df) <- c("tree", "split_type", "confidence", "weights")
  }
  # Trim columns
  results_df <- results_df[ , c("tree", "split_type", "confidence", "weights")]
  
  ## Add any splits from the whole tree, if the two trees have different numbers of tips
  if (whole_df_type == "dataframe"){
    results_df <- rbind(results_df, whole_df)
  }
  
  # Return the number of different splits
  return(results_df)
}


