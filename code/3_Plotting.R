# Open Wu_2018 melted results
df_f <- "/Users/caitlincherryh/Documents/C1_EmpiricalTreelikeness/TestStatistics_BenchmarkAlignments/06_results/treelikeness_scores/Wu_2018_dnaLoci_Primates_completeResults_melted.csv"
df <- read.csv(df_f, stringsAsFactors = FALSE)

library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(patchwork)

# change sCF from 0-100 to 0-1.0 by dividing by 100(divide by 100 so can be easily plotted with treelikeness, 3seq)
inds <- which(df$variable == "sCF_mean")
df$value[inds] <- df$value[inds]/100
inds <- which(df$variable == "sCF_median")
df$value[inds] <- df$value[inds]/100

p_df <- df
p_df <- p_df[p_df$variable %in% c("neighbour_net_trimmed","X3SEQ_prop_recombinant_sequences","sCF_mean"),]
p_df$group <- factor(p_df$variable, levels = c("X3SEQ_prop_recombinant_sequences","neighbour_net_trimmed","sCF_mean"), ordered = TRUE, 
                     labels = c("Proportion of \n recombinant sequences","Tree proportion", "Site concordance \n factor"))
p <- ggplot(p_df, aes(x = group, y = value)) + geom_violin(scale = "area", fill = "grey70") + theme_minimal() + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12)) + 
  labs(x = "\n Test statistic", y = "Value")
p

s_df <- df
s_df <- s_df[s_df$variable %in% c("nn_trimmed_sig","X3SEQ_p_value","sCF_median_sig"),]
s_df$group <- factor(s_df$variable, levels = c("X3SEQ_p_value","nn_trimmed_sig","sCF_median_sig"), ordered = TRUE, 
                     labels = c("3SEQ \n (native p-value)","Tree proportion \n (parametric bootstrap)", "Site concordance factor \n (parametric bootstrap)"))
s <- ggplot(s_df, aes(x = group, y = value)) + geom_violin(scale = "area", fill = "grey70") + theme_minimal() + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12)) + 
  labs(x = "\n p-value", y = "Value")
s

# title = "p-values for treelikeness test statistics in an empirical dataset", subtitle = "16 species of primate from Wu et al (2018)"
p + s + plot_annotation(title = "p-values for treelikeness test statistics in an empirical dataset", subtitle = "16 species of primate from Wu et al (2018)",
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 16)),
                      tag_levels = "a")

cairo_pdf(filename = "/Users/caitlincherryh/Desktop/M1_Wu2018_sample_plot.pdf", width = 13, height = 7, fallback_resolution = 300)
p + s + plot_annotation(title = "p-values for treelikeness test statistics in an empirical dataset", subtitle = "16 species of primate from Wu et al (2018)",
                          theme = theme(plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 16)),
                          tag_levels = "a")
dev.off()


tl_df <- df[df$variable %in% c("neighbour_net_trimmed"),]
seq_df <- df[df$variable %in% c("X3SEQ_prop_recombinant_sequences"),]
scf_df <- df[df$variable %in% c("sCF_mean"),]
tl_sum <- summary(tl_df$value)
seq_sum <- summary(seq_df$value)
scf_sum <- summary(scf_df$value)



tl_p_df <- df[df$variable %in% c("nn_trimmed_sig"),]
seq_p_df <- df[df$variable %in% c("X3SEQ_p_value"),]
scf_p_df <- df[df$variable %in% c("sCF_mean_sig"),]
tl_p_sum <- summary(tl_p_df$value)
seq_p_sum <- summary(seq_p_df$value)
scf_p_sum <- summary(scf_p_df$value)

length(which(tl_p_df$value < 0.05))/length(tl_p_df$value)
# 25% of loci have a statistically significant tree-proportion value

length(which(seq_p_df$value < 0.05))/length(seq_p_df$value)
# 63% of loci have a statistically significant p-value for 3SEQ

length(which(scf_p_df$value < 0.05))/length(scf_p_df$value)
# 5% of scf values have a significant p-value
