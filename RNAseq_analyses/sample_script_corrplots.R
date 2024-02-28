# Sample script using 'airway' test data
# set working directory to source file location

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

# save as PDF:
pdf("corrplot.pdf",width=10,height=10)
corrplot(R, order="hclust", tl.col="black")
dev.off()


# NOTE: I recommend using gene names for this instead




# ***
# Florencia Schlamp, PhD  
# Florencia.Schlamp@nyulangone.org  
# Assistant Director of Bioinformatics
# NYU Cardiovascular Research Center
# New York University Langone Health
# 435 East 30th Street | Science Bldg 604 | New York, NY | 10016

