library(stringr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(ggplotify)
library(RColorBrewer)

# Compare two analyses

# data neded per project:
# 1) gene name / ID
# 2) log fold change
# 3) pvalue for cutoff (adjusted or not)

results_comparison1 <- read.csv("DESeq_results_comparison1.csv", row.names = 1, check.names = F)
results_comparison1_filt <- results_comparison1[complete.cases(results_comparison1$padj) & 
                                      results_comparison1$pvalue < 0.05,]

results_comparison2 <- read.csv("DESeq_results_comparison2.csv", row.names = 1, check.names = F)
results_comparison2_filt <- results_comparison2[complete.cases(results_comparison2$padj) & 
                                          results_comparison2$pvalue < 0.05,]

# step 1) get union of gene names between both comparisons
genes_union <- union(rownames(results_comparison1_filt),
                     rownames(results_comparison2_filt))
subset_comparison1 <- results_comparison1[rownames(results_comparison1) %in% genes_union,]
subset_comparison2 <- results_comparison2[rownames(results_comparison2) %in% genes_union,]

genes_intersect <- intersect(rownames(subset_comparison1),rownames(subset_comparison2))

subset_comparison1 <- subset_comparison1[rownames(results_comparison1) %in% genes_intersect,]
subset_comparison2 <- subset_comparison2[rownames(results_comparison2) %in% genes_intersect,]
nrow(subset_comparison1)
nrow(subset_comparison2)

# step 2) creating new table with all data
data <- data.frame(genes = subset_comparison1$gene_name,
                   comparison1_logFC = subset_comparison1$log2FoldChange,
                   comparison2_logFC = subset_comparison2$log2FoldChange,
                   comparison1_pval = subset_comparison1$pvalue,
                   comparison2_pval = subset_comparison2$pvalue)

# step 3) assign zone to each gene (for colors)
data$coodinate_group = ifelse(data$comparison1_logFC > 0 & data$comparison2_logFC > 0,
                        yes = "up_both",
                        no = ifelse(data$comparison1_logFC < 0 & data$comparison2_logFC < 0,
                                    yes = "down_both",
                                    no = ifelse(data$comparison1_logFC < 0 & data$comparison2_logFC > 0,
                                                yes = "down_comparison1_up_comparison2",
                                                no = "up_comparison1_down_comparison2")))

# step 4) assign significances 
data$pval_color = ifelse(data$comparison1_pval < 0.1 & data$comparison2_pval < 0.1,
                        yes = "sig_both",
                        no = ifelse(data$comparison1_pval < 0.1 & data$comparison2_pval > 0.1,
                                    yes = "sig_comparison1",
                                    no = ifelse(data$comparison1_pval > 0.1 & data$comparison2_pval < 0.1,
                                                yes = "sig_comparison2",
                                                no = "n.s.")))

head(data)
table(data$pval_color)

write.csv(data, file="comparison1_vs_comparison2_table.csv")

p1 <- ggplot(data, aes(y = comparison1_logFC, x = comparison2_logFC, label = genes)) +
  geom_point(aes(color = factor(pval_color)), size = 1.25, alpha = 0.5, na.rm = TRUE) + 
  theme_bw(base_size = 16) + 
  #theme(legend.position = "none") + # remove legend
  ggtitle(label = "Comparison Plot", 
          subtitle= "logFC of sig genes in comparison1 vs. comparison2") +  
  ylab("log2 (Fold Change) in comparison1") + 
  xlab("log2 (Fold Change) in comparison2") + 
  geom_vline(xintercept = 0, colour = "black") + 
  geom_hline(yintercept = 0, colour = "black") + 
  scale_color_manual("Sections",values = c("sig_comparison1"="#a6cee3",
                                           "sig_comparison2"="#fb9a99",
                                           "sig_both"="#33a02c")) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_text(data=data.frame(x=1.3,y=1,label="up_both"),  # edit x and y positions manually depending on the data
            aes(x,y,label=as.character(label)),size=5) +
  geom_text(data=data.frame(x=-3,y=-1,label="down_both"), # edit x and y positions manually depending on the data
             aes(x,y,label=as.character(label)),size=5) +
  geom_text(data=data.frame(x=-10,y=9,label="up_comparison1 \n down_comparison2"),  # edit x and y positions manually depending on the data
            aes(x,y,label=as.character(label)),size=4) +
  geom_text(data=data.frame(x=7.5,y=-7,label="down_comparison1 \n up_comparison2"),  # edit x and y positions manually depending on the data
            aes(x,y,label=as.character(label)),size=4)
p1
ggsave("comparison_volcano_plot_comparison1_vs_comparison2.pdf",plot=p1, width=9, height=7)

# save as interactive file
htmlwidgets::saveWidget(as_widget(ggplotly(p1)), 
                        "comparison_volcano_plot_comparison1_vs_comparison2_interactive.html")


# ***
# Florencia Schlamp, PhD  
# Florencia.Schlamp@nyulangone.org  
# Assistant Director of Bioinformatics
# NYU Cardiovascular Research Center
# New York University Langone Health
# 435 East 30th Street | Science Bldg 604 | New York, NY | 10016
