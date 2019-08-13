### empirical_treelikeness/code/1_TestStatistics_EmpiricalData.R
## R program to apply treelikeness statistics to transcriptomes from empirical data
## Specifically designed to work with Rob Lanfear's "BenchmarkAlignments"
## BenchmarkAlignments and metadata have CC0 or CCBY licenses and are available here: https://github.com/roblanf/BenchmarkAlignments
## A number of additional software packages are required, specifically:
##     - IQTREE (Nguyen et al 2015) (http://www.iqtree.org/)
##     - 3SEQ (Lam et al 2018) (http://mol.ax/software/3seq/)
##     - PHI test (Bruen et al 2006) (https://www.maths.otago.ac.nz/~dbryant/software.html)
##     - Splitstree (Huson and Bryant 2006) (http://www.splitstree.org/)
# Caitlin Cherryh 2019

# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

print("opening packages")
library(ape)
library(parallel)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)
library(TreeSim)

print("initialising namespace")

run_location = "mac"
# run_location = "soma"


if (run_location == "mac"){
  input_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/empiricalData/Wu_2018_dna_loci/"
  output_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/output/output_20190806/"
  treedir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the code is
  maindir <- "/Users/caitlincherryh/Documents/Repositories/empirical_treelikeness/"
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree-1.7-beta13_sCF/bin/iqtree","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
} else if (run_location=="soma"){
  input_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet_Results/"
  treedir <- "/data/caitlin/treelikeness/" # where the code is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","SplitsTree")
}

# Source files for functions
source(paste0(treedir,"code/func_test_statistic.R"))
source(paste0(treedir,"code/func_process_data.R"))
source(paste0(treedir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_empirical.R"))

# Initialise list of species of interest
mammals <- c("CALLI_JAC", "MACAC_FAS", "MACAC_MUL", "PAPIO_ANU", "CHLOR_SAB", "DAUBE_MAD", "GORIL_GOR", "HOMO_SAP", "PAN_PAN", "PAN_TRO", "PONGO_ABE",
             "MICRO_MUR", "NOMAS_LEU", "OTOLE_GAR", "SAIMI_BOL", "TARSI_SYR")
mammals_commonNames <- c("Marmoset", "Crab-eating macaque", "Rhesus macaque", "Olive Baboon", "Green monkey", "Aye-aye", "Gorilla"," Human", "Bonobo",
                         "Chimpanzee", "Orangutan", "Mouse Lemur", "Gibbon", "Galago (Bushbaby)"," Squirrel monkey", "Tarsier")
mammals_speciesNames <- c("Callithrix jacchus", "Macaca fascicularis", "Macaca mulatta", "Papio anubis", "Chlorocebus sabaeus", "Daubentonia madagascariensis",
                          "Gorilla gorilla", "Homo sapiens", "Pan paniscus", "Pan troglodytes", "Pongo abelii", "Microcebus murinus", "Nomascus leucogenys",
                          "Otolemur garnettiiz", "Saimiri boliviensis", "Tarsius syrichta"  )
# Initialise list of outgroups in case rooting is needed later on
outgroups <- c("GALLU_GAL", "MELEA_GAL", "ANOLI_CAR", "PELOD_SIN", "XENOP_TRO", "LATIM_CHA", "DANIO_RER", "GASTE_ACU")
outgroups_commonName <- c("Chicken", "Turkey", "Anole lizard", "Chinese softshell turtle", "Frog", "Coelacanth", "Zebrafish", "Stickleback fish")
outgroups_speciesName <- c("Gallus gallus", "Meleagris gallopavo", "Anolis carolinensis", "Pelodiscus sinensis", "Xenopus tropicalis", "Latimeria chalumnae",
                           "Danio rerio", "Gasterosteus aculeatus")

# Extract the loci from the folder
loci <- list.files(input_dir)
loci <- paste0(input_dir, loci)
# Remove the three alignment files from the folder
loci <- loci[grep("alignment",loci,invert=TRUE)]
loci <- loci[grep("README",loci,invert=TRUE)]
# Create path for whole alignment file
whole_alignment <- paste0(input_dir, "alignment.nex")

# Trim the unecessary species from the alignments
print("trimming alignments")
# Call the trim function to remove all the species you want to remove from the alignment (unless you want to keep ALL the species in the alignment)
# cutSpecies(alignment_path, keep, output_path_provided = "FALSE")
#lapply(loci, cutSpecies, keep = mammals, output_path_provided = FALSE) # trim all the unwanted species from the collected loci
#cutSpecies(alignment_path = whole_alignment, keep = mammals, output_path_provided = FALSE) # trim all the unwanted species from the whole alignment file


# Calculate the test statistics and run the bootstraps
#print("starting analysis")
#print("apply treelikeness test statistics")
# To run for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths,
#                                                        number_of_replicates = 9, iqtree.num_threads = AUTO, iqtree.num_quartets = 100)
# Parameter choice:
#       ~ iqtree.num_thread = 1: allows parallelisation higher up in the workflow
#       ~ iqtree.num_quartets = 1000: greater than 100 quartets necessary for stable sCF values

if (run_location=="soma"){
  #lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 99, iqtree.num_threads = 1, iqtree.num_quartets = 1000) 
} else if (run_location=="mac"){
  #lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 5, iqtree.num_threads = 1, iqtree.num_quartets = 1000) 
}

#print("collate resulst")
# Collate all the dataframes together
#results_file <- paste0(output_dir,basename(BA_dir),"_completeResults.csv")
#df <- collate.bootstraps(directory = BA_dir, file.name = "pValues", id = "", output.file.name = results_file)

# test alignment for sampling
test_al <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/tests/test_03_sCFAndMore/ENSG00000000419dna.nex"
test_als <- paste0("/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/tests/test_03_sCFAndMore/",list.files("/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/tests/test_03_sCFAndMore/"))
test_mldist <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/tests/test_03_sCFAndMore/ENSG00000000419dna.nex" # don't include mldist

# test mst
# pdm <- mldist.pdm(test_mldist) 

# for testing non-parallel
#sCF(iqtree_path = exec_paths["IQTree"], alignment_path = test_al, num_threads = "AUTO", num_quartets = 1000)

# for testing whole implementation
#lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 5, iqtree.num_threads = 1, iqtree.num_quartets = 1000) 

scf_file <- '/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/tests/test_03_sCFAndMore/ENSG00000000419dna.nex.treefile.cf.stat'

# testing bootstrap etc
#lapply(test_als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 5, iqtree.num_threads = 1, iqtree.num_quartets = 1000) 
#empirical.runTS(alignment_path = test_al, program_paths = exec_paths, bootstrap_id = "alignment", iqtree.num_threads = 1, iqtree.num_quartets = 1000)
#empirical.bootstraps.wrapper(empirical_alignment_path = test_al, program_paths = exec_paths, number_of_replicates = 9, iqtree.num_threads = 1, iqtree.num_quartets = 1000)

