# Sample script using 'airway' test data
# set working directory to source file location

# load data from previous analyses
load("DESeq_results.Rdata")

### (1) Translate IDs to gene names ###
input_IDs <- rownames(res)

head(input_IDs)

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") <-- if you have mouse

## alternatives specifiying host:
#marts <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org") # host="useast" or "uswest"
#human <- useDataset("hsapiens_gene_ensembl", mart=marts) <-- if you have human
#mouse <- useDataset("mmusculus_gene_ensembl", mart=marts) <-- if you have mouse

# if you want to see all possible attributes and filters you can use:
humanatts <- listAttributes(human)
humanfilters <- listFilters(human)

gene_convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)
# gene_convert is now your dictionary

head(res)

# to create new column with gene names, you can map values using this code:
res$gene_name <- as.character(plyr::mapvalues(x=rownames(res),
                                              from=gene_convert$ensembl_gene_id,
                                              to=gene_convert$external_gene_name))

head(res)
# NOTE: if there is no translation available, this 'plyr::mapvalues' function
# will keep the original name. example:
tail(res) # look at gene "ENSG00000283111"


# you can get other information (check listAttributes/listFilters options above)
# (2) For example getting chromosome numbers as well:
gene_convert <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name"),
                      filters="ensembl_gene_id",
                      values=input_IDs, mart=human)

head(gene_convert)

res$chr <- as.character(plyr::mapvalues(x=rownames(res),
                                        from=gene_convert$ensembl_gene_id,
                                        to=gene_convert$chromosome_name))
head(res)
tail(res) # WARNING: same as before, if it can't find a gene,
# it will keep the original ID (look at "ENSG00000283111"). 
# In this case it doesn't make sense to keep the gene ID in the chr column
# better alternative in next example


### (3) Get human homolog genes for mouse genes (or viceversa) ###
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mouseatts <- listAttributes(mouse)

gene_convert = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                      values = input_IDs, mart = human, 
                      attributesL = c("external_gene_name"), martL = mouse)

head(gene_convert)

# this is a different way of creating a new column
# with this method, if the value is not there, it will add 'NA' instead
results_table <- as.data.frame(res)
results_table$geneID <- rownames(results_table)
colnames(gene_convert)
colnames(gene_convert) <- c("geneID","mouse_gene_name")
results_table <- merge(results_table,gene_convert,all=T)
head(results_table)
tail(results_table) # lots of NAs
# NOTE: this method is better for adding chromosome information in previous example

# IMPORTANT NOTE: the reason we don't want NAs when translating from geneID to gene_name
# is because this translation if very useful to create new data frames 
# (for volcano plots for example!) with gene_name as rownames. And rownames cannot be NAs.

# FINAL NOTE: this translation is not final/static. As the biomart database is updated,
# the translations can change too. So always save the 'gene_convert' dictionary as a file
# ideally with the date when you queried biomart on the name.


# ***
# Florencia Schlamp, PhD  
# Florencia.Schlamp@nyulangone.org  
# Sr. Bioinformatics Programmer  
# Division of Cardiology  
# New York University Langone Health  
# 435 East 30th Street | Science Bldg 604 | New York, NY | 10016

