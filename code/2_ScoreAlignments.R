# trialling aliscore function
library(ips)
alipath <- "/Users/caitlincherryh/Documents/Honours/Executables/Aliscore_v.2.0/Aliscore.02.2.pl"
#alipath <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/Aliscore.02.2.pl"
p <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/COI.nex"
op <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/16S.nex.fasta_List_random.txt"
n <- read.nexus.data(p)
d <- as.DNAbin(n)


aliscore(alignment_path = p, w = 6, r = 10, aliscore_path = alipath, quality_threshold = 0.5, redo = TRUE)