### empirical_treelikeness/code/func_comparison.R
## R functions for comparing two species trees
# Caitlin Cherryh 2021

# Function to read in tree from ASTRAL and edit it to use in the QuartetNetworkGoodnessOfFit julia package
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
