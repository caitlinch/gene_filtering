# R functions to check the quality of a given alignment

library(ape)

# This function evaluates whether a given locus has a high or low quality alignment
# This function runs ALISCORE on a given alignment and outputs whether the alignment is high or low quality based on the selected quality threshold (proportion of randomly similar sites from 0 to 1)
aliscore <- function(alignment_path, output_path, gaps = "5char", w, r, tree_path, l, s, o, aliscore_paths, quality_threshold = 0.5, redo = FALSE){
  print(paste0("checking ",alignment_path))
  # Extract the names of the dataset and the loci name so you can output them into a nice dataframe at the end
  loci_name <- gsub(".nex","",basename(alignment_path))
  dataset <- basename(dirname(alignment_path))
  output_path <- ifelse(missing(output_path),paste0(dirname(alignment_path),"/"),output_path) # if output_path is missing, set it to the same directory as the alignment
  # Make a file name to save the output csv with 
  output_csv_name <- paste0(output_path,loci_name,"_ALISCORE_locusAlignmentQuality.csv")
  if ((file.exists(output_csv_name) == TRUE) && (redo == FALSE)){
    # If the csv file already exists and you don't want to redo the analysis, this just skips that particular alignment
    print(paste0(dataset, " - ", loci_name," - alignment skipped"))
  } else {
    # If you do want to redo the analysis OR the analysis hasn't been done yet, this will run the alignment and create the quality check output csv for the alignment
    print(paste0(dataset, " - ", loci_name," - alignment processing"))
    # Copy the ALISCORE perl file and the module file into this folder if they are not already here.
    aliscore_file <- paste0(dirname(alignment_path),"/Aliscore.02.2.pl")
    module_file <- paste0(dirname(alignment_path),"/Aliscore_module.pm")
    if (file.exists(aliscore_file) == FALSE){
      file.copy(aliscore_paths[1],aliscore_file)
    }
    if (file.exists(module_file) == FALSE){
      file.copy(aliscore_paths[2],module_file)
    }
    # Check if this alignment has already been processed: if it has already been run, and redo = FALSE, skip the alignment
    # Convert alignment to a fasta (with no "?" - replace them with "N" so ALISCORE will run) if it isn't already
    fasta_path <- check.fasta.format(alignment_path)
    # Change the wd to where you want the files created by ALISCORE to go
    setwd(output_path)
    # Form the ALISCORE command, collect the information to feed into the ALISCORE program
    N <- ifelse(gaps == "5char","","-N") # gaps treated as 5th character if the -N option is not invoked (here if value for gaps is set a "5char") -> -N means gaps are ambiguous character
    w <- ifelse(missing(w), "-w 6", paste("-w", w)) # window size
    r <- ifelse(missing(r), "", paste("-r", r)) # number of random pairs to use - if -r used with no argument or -r not used, the program selectes 4*N pairs where N = number of taxa
    t <- ifelse(missing(tree_path), "", paste("-t",tree_path)) # the path to the tree if using a tree
    l <- ifelse(missing(l), "", paste("-l", l)) # the node level to restrict iterating through the tree to
    s <- ifelse(missing(s), "", "-s") # generate strict profile - 1 no = no (not median)
    o <- ifelse(missing(o), "", paste("-o", paste(o, collapse = ","))) # outgroups
    # Run the ALISCORE program
    command <- paste("perl -I ",module_file,aliscore_file,"-i",fasta_path,N,w,r,t,l,s,o)
    print(command)
    system(command)
    # Open the file that contains the list of sites that are randomly similar
    files <- list.files(dirname(alignment_path))
    files <- files[grep(basename(fasta_path),files)]
    op_file <- files[grep("List",files)]
    print(paste0("opening output: ",op_file))
    # Extract the number of sites and taxa from the nexus file
    params <- get.nexus.parameters(alignment_path)
    nchar <- params$n_chars
    # If the file exists, create a nice row of results
    if (length(op_file) >= 1){
      # If multiple hits, select the first one
      if (length(op_file)>1){
        op_file <- op_file[1]
      }
      # Open the list of randomly similar sites
      rss_sites <- as.numeric(strsplit(readLines(op_file),split=" ")[[1]])
      # Calculate the proportion of sites that are randomly similar
      proportion_rss <- length(rss_sites)/nchar
      # If the proportion of sites that are randomly similar is higher than the quality threshold, the alignment passes the quality check
      if (proportion_rss > quality_threshold){
        quality_check <- "FAIL"
        print(paste0(dataset, " - ", loci_name," - loci FAILED quality test : ",proportion_rss))
      } else if (proportion_rss <= quality_threshold){
        quality_check <- "PASS"
        print(paste0(dataset, " - ", loci_name," - loci PASSED quality test : ",proportion_rss))
      }
    } else {
      # If no output file exists, record that into the csv
      print(paste0(dataset, " - ", loci_name," - loci run error"))
      # If the file doesn't exist, create a nice row of missing values
      rss_sites <- "RUN_ERROR"
      proportion_rss <- "RUN_ERROR"
      quality_check <- "NA"
    }
    # Output the information about the alignment and the results of the quality check using ALISCORE
    df <- data.frame(dataset,loci_name,alignment_path,params$n_taxa,nchar,length(rss_sites),proportion_rss,quality_check)
    names(df) <- c("dataset","loci","alignment_file","n_taxa","n_sites","n_rss_sites","proportion_rss_sites","loci_quality_assessment")
    # Save the little dataframe in the right place
    write.csv(df, file = output_csv_name, row.names = FALSE)
  }
}



# This function makes sure an alignment is in fasta format
check.fasta.format <- function(alignment_path){
  # This alignment returns a fasta format of a given alignment
  # If an alignment is in fasta format, the alignment name is returned
  # If an alignment is nexus format, it is converted to fasta and the fasta name returned
  filetype = tail(strsplit(alignment_path,"\\.")[[1]],n=1) # extract file format
  if (filetype == "fasta"|| filetype == "fas"){
    # if the alignment is already in fasta format, return the file
    fasta.name <- paste0(gsub(paste0(".",filetype),"",alignment_path),"_noq.fasta")
    fasta.noq(fasta.name)
    return(alignment_path)
  } else if (filetype == "nexus" || filetype == "nex"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = as.DNAbin(read.nexus.data(alignment_path)) # read in nexus format alignment
    fasta.name <- paste0(alignment_path,"_noq.fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.FASTA(x = data, file = fasta.name) # output alignment as a fasta format
    fasta.noq(fasta.name)
    return(fasta.name)
  }
}



# This function takes the "?" in a file and changes them to "N" so that ALISCORE will run
fasta.noq <- function(alignment_path){
  fasta_file <- readLines(alignment_path)
  fasta_file <- gsub("\\?","N",fasta_file)
  writeLines(fasta_file,alignment_path)
}



# Given a nexus file, this function extracts the taxa, number of characters and number of taxa
get.nexus.parameters <- function(nexus_file){
  # Open the nexus file and get the information you want from it
  n <- read.nexus.data(nexus_file)
  taxa <- names(n)
  n_taxa <- length(taxa)
  n_chars <- length(n[[1]])
  # stick the params together and return them
  params <- list(taxa,n_taxa,n_chars)
  names(params) <- c("taxa","n_taxa","n_chars")
  return(params)
}



# Function to check whether an alignment has been processed by ALISCORE, and return the alignment path if it hasn't
ALISCORE.output.exists <- function(al){
  # construct the file name
  loci_name <- gsub(".nex","",basename(al))
  dataset <- basename(dirname(al))
  output_path <- paste0(dirname(al),"/")
  output_csv_name <- paste0(output_path,loci_name,"_ALISCORE_locusAlignmentQuality.csv")
  # Check whether that file exists
  if (file.exists(output_csv_name) == TRUE){
    # return nothing because the file has already been run
    file_result <- character(0)
  } else if (file.exists(output_csv_name) == FALSE){
    # return the alignment because it hasn't been run
    file_result <- al
  }
  return(file_result)
}
  