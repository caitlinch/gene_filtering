# trialling aliscore function
library(ips)
alipath <- "/Users/caitlincherryh/Documents/Honours/Executables/Aliscore_v.2.0"
p <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/test_01_aliscoretest/Anderson_2013/16S.nex"
n <- read.nexus.data(p)
d <- as.DNAbin(n)
aliscore(x = d, gaps = "ambiguous", w = 10, o = c("8CH"), path = alipath)

# ALISCORE is only able to work in the same folder as itself