# Sample script using 'airway' test data
# set working directory to source file location

library(ggplot2)
library(iCellR)
library(ggplotify)
library(ggrepel)

load("DESeq_results.Rdata")

#### === make volcano plots === #####

DEresults <- as.data.frame(res)
## merge DE results and counts table ##
# this only works if both tables have the same rows (same number, same name, same order)
all_data <- cbind(DEresults,counts_data_norm_log2)

# data to plot:
head(all_data)

# get rid of genes with NAs for pvalue 
# (or padj, depending on what you want to use)
filtered_data <- all_data[complete.cases(all_data$pvalue),]


## OPTIONAL: filtering steps ##
# 1) filter by basemeans:
filtered_data2 <- filtered_data[filtered_data$baseMean > 30,]

# 2) filter by norm counts:
colnames(filtered_data) # the column numbers that correspond to norm counts (in this case) are 7 to 14
keep <- rowSums((filtered_data[,c(7:14)])>6) >= 2
filtered_data3 <- filtered_data[keep,]

# compare distributions:
hist(as.matrix(filtered_data[,c(7:14)])) # all data
hist(as.matrix(filtered_data2[,c(7:14)])) # filtered by basemeans
hist(as.matrix(filtered_data3[,c(7:14)])) # filtered by norm counts


## volcano plot data
Pval = -log10(0.05)
dot.col = c("#E64B35","#3182bd","#636363")

data <- data.frame(gene = rownames(filtered_data), # choose here if you want geneID or gene_name
                   pvalue = -log10(filtered_data$pvalue), # choose here if you want pvalue or padj
                   lfc = filtered_data$log2FoldChange)
head(data) # <- this is what will be plotted

# assign up and down-regulted colors
data <- data %>% mutate(color = ifelse(data$lfc > 0 & data$pvalue > Pval,
                                       yes = "upregulated",
                                       no = ifelse(data$lfc < 0 & data$pvalue > Pval,
                                                   yes = "downregulated",
                                                   no = "not-significant")))

## Volcano plot with genes of interest labeled:

genes_of_interest <- c("ENSG00000167641","ENSG00000131242","ENSG00000001084") 

subset <- data[data$gene %in% genes_of_interest,]
subset$color <- "highlight"
head(subset)

ggplot(data, aes(x = lfc, y = pvalue, text = gene)) +
  geom_point(aes(color = factor(color)), size = 1.25, alpha = 0.5, na.rm = TRUE) + 
  theme_bw(base_size = 16) + 
  theme(legend.position = "none") +
  ggtitle(label = "Volcano Plot", subtitle= "treatment vs control") +  # add title
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (p-value)") + # 
  geom_vline(xintercept = 0, colour = "black") + 
  geom_hline(yintercept = Pval, colour = "black") + 
  scale_color_manual("Color Key", values = c("upregulated" = dot.col[1],
                                             "downregulated" = dot.col[2],
                                             "not-significant" = dot.col[3],
                                             "highlight"="black")) +
  scale_y_continuous(trans = "log1p") + 
  geom_text_repel(data = subset, mapping = aes(label = gene, color=factor(color)),
                  size = 5,
                  fontface = 'bold',
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines")) + # add labels to chosen genes
  ggsave(paste0("volcano_plot_with_labels.pdf"), width=6, height=6)




# ***
# Florencia Schlamp, PhD  
# Florencia.Schlamp@nyulangone.org  
# Assistant Director of Bioinformatics
# NYU Cardiovascular Research Center
# New York University Langone Health
# 435 East 30th Street | Science Bldg 604 | New York, NY | 10016

