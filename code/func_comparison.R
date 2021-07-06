### empirical_treelikeness/code/func_comparison.R
## R functions for comparing two species trees
# Caitlin Cherryh 2021

# Function to read in tree from ASTRAL and edit it to use in the QuartetNetworkGoodnessOfFit Julia package
# ASTRAL does not output terminal branch lengths, so these are arbitrarily given the length 1
# ASTRAL includes posterior probabilities: these are stripped (the Julia function cannot handle extra node values)
reformat.ASTRAL.tree.for.Julia <- function(tree_path){
  # Open tree
  t <- read.tree(tree_path)
  # Identify which edges correspond to terminal branches
  n = length(t$tip.label)
  missing_branch_indexes <- sapply(1:n,function(x,y)   which(y==x),y=t$edge[,2])
  # Set length of terminal branches to 1
  t$edge.length[missing_branch_indexes] <- 1
  # Don't want node labels - for qaurnetGoFtest, need just taxa and branch lengths
  t$node.label <- NULL
  # Write this tree back out
  write.tree(t, tree_path)
}



# Function to read in list of genetrees from ASTRAL and edit to use in the QuartetNetworkGoodnessOfFit Julia package
# ASTRAL does not output terminal branch lengths, so these are arbitrarily given the length 1
# ASTRAL includes posterior probabilities: these are stripped (the Julia function cannot handle extra node values)
reformat.gene.tree.list.for.Julia <- function(trees_path, gene.tree.source = "IQ-TREE"){
  # Open trees
  ts <- read.tree(trees_file)
  for (i in 1:length(ts)){
    temp_tree <- ts[[i]]
    temp_tree <- reformat.phylo.for.Julia(temp_tree, gene.tree.source)
    ts[[i]] <- temp_tree
  }
  # Write these trees back out
  write.tree(ts, gsub(".txt", "_test.txt", trees_path))
}

# Function to take one tree from a multiphylo object, format it for QuartetNetworkGoodnessOfFit Julia package, and return it
reformat.phylo.for.Julia <- function(tree, gene.tree.source = "IQ-TREE"){
  # If gene trees from ASTRAL, the terminal branches are not output and need to be set to one
  if (gene.tree.source == "ASTRAL"){
    # Identify which edges correspond to terminal branches
    n = length(tree$tip.label)
    missing_branch_indexes <- sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])
    # Set length of terminal branches to 1
    tree$edge.length[missing_branch_indexes] <- 1
  }
  # Gene trees from ASTRAL have posterior probabilities, and those from IQ-Tree have bootstraps.
  # Either way, don't want node labels for quarnetGoFtest - need just taxa and branch lengths
  tree$node.label <- NULL
  return(tree)
}
