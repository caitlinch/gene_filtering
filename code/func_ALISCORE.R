library(ape)

# trialling aliscore function
library(ips)
alipath <- "/Users/caitlincherryh/Documents/Honours/Executables/Aliscore_v.2.0/Aliscore.02.2.pl"
alipath <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/Aliscore.02.2.pl"
p <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/16S.nex"
n <- read.nexus.data(p)
d <- as.DNAbin(n)


# ALISCORE is only able to work in the same folder as itself
aliscore(alignment_path = p, w = 6, r = 100, aliscore_path = alipath)

# This function runs ALISCORE on a given alignment
aliscore <- function(alignment_path, output_path, gaps = "5char", w, r, tree_path, l, s, o, aliscore_path){
  # convert alignment to a fasta if it isn't already
  fasta_path <- check.fasta.format(alignment_path)
  output_path <- ifelse(missing(output_path),paste0(dirname(alignment_path),"/"),output_path) # if output_path is missing, set it to the same directory as the alignment
  #setwd(output_path)
  N <- ifelse(gaps == "5char","","-N") # gaps treated as 5th character if the -N option is not invoked - -N means gaps are ambiguous character
  w <- ifelse(missing(w), "-w 6", paste("-w", w)) # window size
  r <- ifelse(missing(r), "", paste("-r", r)) # number of random pairs to use
  t <- ifelse(missing(tree_path), "", paste("-t",tree_path)) # the path to the tree if using a tree
  l <- ifelse(missing(l), "", paste("-l", l)) # the node level to restrict iterating through the tree to
  s <- ifelse(missing(s), "", "-s") # generate strict profile - 1 no = no (not median)
  o <- ifelse(missing(o), "", paste("-o", paste(o, collapse = ","))) # outgroups
  command <- paste("perl",aliscore_path,"-i",fasta_path,N,w,r,t,l,s,o)
  print(command)
  system(command)
}

# This function makes sure an alignment is in fasta format
check.fasta.format <- function(alignment_path){
  # This alignment returns a fasta format of a given alignment
  # If an alignment is in fasta format, the alignment name is returned
  # If an alignment is nexus format, it is converted to fasta and the fasta name returned
  filetype = tail(strsplit(alignment_path,"\\.")[[1]],n=1) # extract file format
  if (filetype == "fasta"){
    # if the alignment is already in fasta format, return the file
    return(alignment_path)
  } else if (filetype == "nexus" || filetype == "nex"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = as.DNAbin(read.nexus.data(alignment_path)) # read in nexus format alignment
    fasta.name <- paste0(alignment_path,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.FASTA(x = data, file = fasta.name) # output alignment as a fasta format
    return(fasta.name)
  }
}
