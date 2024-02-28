# Sample script using 'airway' test data
# set working directory to source file location

library(ggplot2)
library(stringr)
library(edgeR)
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(iCellR)
library(ggplotify)
library(ggrepel)

##### 'manual' exploration of data #####
#1) read in data and explore
raw_counts <- read.csv("airway_scaledcounts.csv",row.names = 1)
head(raw_counts)
dim(raw_counts)

sample_table <- read.csv("airway_metadata.csv")
sample_table
dim(sample_table)

# lets use better sample names
better_sample_names <- c("C1","T1","C2","T2","C3","T3","C4","T4")
sample_table$sample_name <- better_sample_names
sample_table

head(raw_counts)
colnames(raw_counts)
colnames(raw_counts) <- better_sample_names

# explore library size
libsizes=colSums(raw_counts)
summary(libsizes/1000000)

data_to_plot <- data.frame(sample_name = sample_table$sample_name,
                           libsize = libsizes)

ggplot(data_to_plot, aes(x=sample_name, y=libsizes/1000000)) +
  geom_bar(stat="identity") + 
  labs(y="Library size (in million of reads)", title="Library size per sample (in million of reads)")

ggplot(data_to_plot, aes(x = reorder(sample_name, -libsizes/1000000), y=libsizes/1000000)) +
  geom_bar(stat="identity") + coord_flip() +
  labs(y="Library size (in million of reads)", title="Library size per sample (in million of reads)") +
  geom_hline(yintercept=summary(libsizes/1000000)[[4]]) #as.numeric(summary(libsizes/1000000)["Mean"])

#2) remove genes that have 0 counts across all samples
dim(raw_counts)
table(rowSums(raw_counts) != 0)

sel.rn=rowSums(raw_counts) != 0
raw_counts_filt=raw_counts[sel.rn,]
dim(raw_counts_filt)

# save data for later:
save(raw_counts_filt,sample_table,file="StartData.Rdata")

#3) normalize by library size
head(raw_counts_filt)

libsizes=colSums(raw_counts_filt)
size.factor=libsizes/exp(mean(log(libsizes)))
norm_counts=t(t(raw_counts_filt)/size.factor)

head(norm_counts)

#4) log scaled
range(raw_counts_filt)
range(norm_counts)

norm_log_counts = log2(norm_counts+1)

head(norm_log_counts)
range(norm_log_counts)
hist(norm_log_counts)

plotMDS(norm_log_counts)

hist(norm_log_counts,breaks=30, main="histogram of normalized counts")

keep=rowSums((norm_log_counts)>6) >= 2
norm_log_counts_filt <- norm_log_counts[keep,]

dim(norm_log_counts_filt)

hist(norm_log_counts_filt,breaks=30, main="histogram of normalized counts")

plotMDS(norm_log_counts_filt)

# distances
sampleDists <- dist(t(norm_log_counts_filt))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###### DESeq2 #######
#1) exploration of data

ddsMat <- DESeqDataSetFromMatrix(countData = raw_counts_filt,
                                 colData = sample_table,
                                 design = ~ 1)

rld <- rlog(ddsMat) # takes some time to run

# basic PCA
plotPCA(rld, intgroup=c("dex","celltype"))

# better PCA plot:
data <- plotPCA(rld, intgroup = c("dex","celltype"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=dex)) + 
  geom_point(aes(shape=celltype), size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        axis.title.x=element_text(size=18, margin=margin(20,0,0,0)),
        axis.title.y=element_text(size=18, margin=margin(0,20,0,0)))

# 2) differential expression
load("StartData.Rdata")

ddsMat <- DESeqDataSetFromMatrix(countData = raw_counts_filt,
                                 colData = sample_table,
                                 design = ~ dex) # choose comparison formula here

ddsMat = DESeq(ddsMat) # normalization and differential expression steps are done here
resultsNames(ddsMat) # print result terms

res = results(ddsMat, name = "dex_treated_vs_control") # extract results for desired term
res # <--- save this table for DE results

table(res$padj < 0.05)
# 2175 genes have a significance of padj < 0.05

head(res[order(res$padj),], 10) # there are the top 10 most significant genes

## top 10 genes
list_of_interest <- rownames(head(res[order(res$padj),],10))

# extract norm counts
counts_data_norm <- counts(ddsMat, normalized=TRUE)
counts_data_norm_log2 <- log2(1 + counts_data_norm)
dim(counts_data_norm_log2)
head(counts_data_norm_log2)

# save data for later
save(res,counts_data_norm_log2,file="DESeq_results.Rdata")

# plot top 10 genes:
start_table <- sample_table
datalist = list()

for(gene_name in list_of_interest){
  if(gene_name %in% rownames(counts_data_norm_log2)){
    gene_counts <- melt(counts_data_norm_log2[rownames(counts_data_norm_log2) == gene_name,])$"value"
    start_table["counts"] <- gene_counts
    start_table["Gene"] <- c(rep(gene_name,length(gene_counts)))
    start_table["gene_ID"] <- c(rep(gene_name,length(gene_counts)))
    start_table["type"] <- c(rep("norm_counts",length(gene_counts)))
    datalist[[gene_name]] <- start_table
  }
}

big_data = do.call(rbind,datalist)

ggplot(big_data, aes(x=dex, y=counts, color=Gene, group=Gene)) +
  geom_point(size=2) + 
  stat_summary(fun.y = mean, geom="line", mapping = aes(group=Gene), size=1) +
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title.x=element_text(size=12, margin=margin(5,0,0,0)),
        axis.title.y=element_text(size=12, margin=margin(0,5,0,0))) +
  ylab("Normalized Counts") +
  xlab("treatment") +
  labs(title="Normalized counts of top 10 genes")


## filter NAs
table(complete.cases(res))
# 9917 genes have NAs, should be removed

res_filtered <- res[complete.cases(res),]
dim(res)
dim(res_filtered)

## volcano plot
Pval = -log10(0.05)
dot.col = c("#E64B35","#3182bd","#636363")

# remove genes with NA pvalues
filtered_data <- res_filtered

data <- data.frame(gene = rownames(filtered_data),
                   pvalue = -log10(filtered_data$padj),
                   lfc = filtered_data$log2FoldChange)

data <- data %>% mutate(color = ifelse(data$lfc > 0 & data$pvalue > Pval,
                                       yes = "upregulated",
                                       no = ifelse(data$lfc < 0 & data$pvalue > Pval,
                                                   yes = "downregulated",
                                                   no = "not-significant")))

# Color corresponds to fold change directionality
ggplot(data, aes(x = lfc, y = pvalue, text = gene)) +
  geom_point(aes(color = factor(color)), size = 1.25, alpha = 0.5, na.rm = TRUE) + # add gene points
  theme_bw(base_size = 16) + # clean up theme
  #theme(legend.position = "none") + # remove legend
  ggtitle(label = "Volcano Plot", subtitle= "treatment vs control") +  # add title
  xlab("log2 (Fold Change)") + # x-axis label
  ylab("-log10 (adjusted p-value)") + # y-axis label
  geom_vline(xintercept = 0, colour = "black") + # add line at 0
  geom_hline(yintercept = Pval, colour = "black") + # p(0.05) = 1.3
  scale_color_manual("Color Key", values = c("upregulated" = dot.col[1],
                                             "downregulated" = dot.col[2],
                                             "not-significant" = dot.col[3])) +
  scale_y_continuous(trans = "log1p") + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggsave(paste0("volcano_plot.pdf"), width=8, height=6) # add this line to save plot as PDF

# run plot again but save as variable instead of plotting
p1 <- ggplot(data, aes(x = lfc, y = pvalue, text = gene)) +
  geom_point(aes(color = factor(color)), size = 1.25, alpha = 0.5, na.rm = TRUE) + 
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle= "treatment vs control") + 
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") + 
  geom_vline(xintercept = 0, colour = "black") + 
  geom_hline(yintercept = Pval, colour = "black") + 
  scale_color_manual("Color Key", values = c("upregulated" = dot.col[1],
                                             "downregulated" = dot.col[2],
                                             "not-significant" = dot.col[3])) +
  scale_y_continuous(trans = "log1p") + 
  guides(colour = guide_legend(override.aes = list(size=5)))

# save plot as interactive HTML
htmlwidgets::saveWidget(as.widget(ggplotly(p1)), "volcano_plot_interactive.html")


## Filter genes of interest
# lets start with significant genes (padj < 0.05), with a logFC cutoff of < -2 or > 2
res_of_interest <- as.data.frame(res_filtered[res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 2,])
dim(res_of_interest)
# we end up with 169 genes
genes_of_interest <- rownames(res_of_interest)

# subset norm counts for genes of interest:
norm_counts_subset <- counts_data_norm_log2[genes_of_interest,]

## Row Z-score scaled normalized counts heatmap
# TEST: simple heatmap:
pheatmap(norm_counts_subset)

# better heatmap:
pheatmap(norm_counts_subset,
         show_rownames=FALSE,  # too many genes to have the row names
         scale="row",    # row Z-score transformation to show max/min expression
         angle_col=0) # change angle of sample labels

# add column membership
sample_table
column_annotation <- data.frame(row.names=sample_table$sample_name,
                                treatment=sample_table$dex)

pheatmap(norm_counts_subset,
         show_rownames=FALSE,  
         scale="row",    
         angle_col=0,
         annotation_col=column_annotation)

# TEST: could add cell type annotation as well
column_annotation2 <- data.frame(row.names=sample_table$sample_name,
                                treatment=sample_table$dex,
                                cell_type=sample_table$celltype)
pheatmap(norm_counts_subset,
         show_rownames=FALSE,  
         scale="row",    
         angle_col=0,
         annotation_col=column_annotation2)

# TEST: how does it look like without column clustering?
pheatmap(norm_counts_subset,
         show_rownames=FALSE,  
         scale="row",    
         angle_col=0,
         annotation_col=column_annotation2,
         cluster_cols=FALSE)

## final heatmap for now:
pheatmap(norm_counts_subset,
         show_rownames=FALSE,  
         scale="row",    
         angle_col=0,
         annotation_col=column_annotation)

# another way of saving plots as PDFs:
pdf("heatmap_plot.pdf",width=5,height=6,paper="special")
pheatmap(norm_counts_subset,
         show_rownames=FALSE,  
         scale="row",    
         angle_col=0,
         annotation_col=column_annotation)
dev.off() # some times you have to run this line multiple times, until it says 'null device'



# get subsets of genes for functional enrichment:
subset_of_genes <- rownames(res_of_interest)

write.csv(subset_of_genes,file="top_expressed_genes.csv")

# 1) copy paste IDs to website here: http://www.pantherdb.org/
# 2) choose 'homo sapiens'
# 3) choose 'Statistical annotation set' with the annotation set 'GO biological process complete'

# in next page, it asks for select reference list. you can upload the whole list of genes before the subset,
# or you can just choose on the left "Use a Reference List that includes all genes from a whole genome", and pick human genome



### translate IDs to gene names ###
input_IDs <- as.character(rownames(raw_counts_filt))

head(input_IDs)

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") <-- if you have mouse

## alternatives specifiying host:
#marts <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org") # host="useast" or "uswest"
#human <- useDataset("hsapiens_gene_ensembl", mart=marts) <-- if you have human
#mouse <- useDataset("mmusculus_gene_ensembl", mart=marts) <-- if you have mouse

# if you want to see all possible attributes and filters you can use:
mouseatts <- listAttributes(mouse)
mousefilters <- listFilters(mouse)

gene_convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)
# gene_convert is now your dictionary

# to create new column with gene names, you can map values using this code:
res$gene_name <- as.character(plyr::mapvalues(x=res$geneID,
                                              from=gene_convert$ensembl_gene_id,
                                              to=gene_convert$external_gene_name))


# you can get other information (check listAttributes/listFilters options above)
# for example getting chromosome numbers as well:
gene_convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)

### get human homolog genes for mouse genes (or viceversa) ###
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_convert = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                      values = input_IDs, mart = mouse, 
                      attributesL = c("hgnc_symbol"), martL = human)



##### plot correlation between genes #####
library(corrplot)

load("DESeq_results.Rdata")

## pick genes of interest
list_of_interest <- rownames(head(res,25))

# get norm counts for chosen genes
data <- counts_data_norm_log2[rownames(counts_data_norm_log2) %in% list_of_interest,]

# transpose data
data <- t(data)

# function to get p-values
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# get matrix of correlation values
R <- cor(data)

# get matrix of the p-value of the correlation
p.mat <- cor.mtest(data)

# simple correlation plot
corrplot(R)

# ordered correlation plot
corrplot(R, order="hclust") # much easier to interpret!

# more options:
corrplot(R, order="hclust", tl.col="black", type="upper", diag=FALSE)

# if you want to know if the correlation is significant:
corrplot(R, order="hclust", tl.col="black", type="upper", diag=FALSE,
         p.mat=p.mat, sig.level=0.05, insig="blank")
# try changing the sig.level (p-value cutoff)

# NOTE: I recommend using gene names for this instead



# ***
# Florencia Schlamp, PhD  
# Florencia.Schlamp@nyulangone.org  
# Sr. Bioinformatics Programmer  
# Division of Cardiology  
# New York University Langone Health  
# 435 East 30th Street | Science Bldg 604 | New York, NY | 10016

