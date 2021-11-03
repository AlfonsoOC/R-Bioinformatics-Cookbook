
count_dataframe <- readr::read_tsv(file.path(getwd(),"/..","datasets", "ch1", "modencodefly_count_table.txt" ))

genes <- count_dataframe[['gene']]
count_dataframe[['gene']] <- NULL

count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes


pheno_data <- readr::read_table2(file.path(getwd(),"/..", "datasets", "ch1", "modencodefly_phenodata.txt")) #here they use readtable2 from readr function, nut sure why

experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which( pheno_data[['stage']] %in% experiments_of_interest ) #gives the numeric indices of the above experiments in the matrix

library("dplyr")
library("magrittr")
grouping <- pheno_data %>% 
  dplyr::filter(stage %in% experiments_of_interest )

counts_of_interest <-  count_matrix[,columns_of_interest] #using the indexing to extract the columns we want to analize


library("DESeq2")
# convert the matrix in DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_of_interest,
                              colData = grouping,
                              design = ~ stage)

#Performing the Data Analysis
dds <- DESeq(dds)

#For Results: This shows the mean counts, the log2 fold change between samples for a gene, the standard error of the log2 fold change, the Wald statistic, and the raw and adjusted P value.
#The padjcolumn for adjusted P values is the one most commonly used for concluding about significance.
res <- results(dds, contrast=c("stage","L2Larvae","L1Larvae"))



library(SummarizedExperiment)
load(file.path(getwd(), "/../datasets/ch1/modencodefly_eset.RData"))
summ_exp <- makeSummarizedExperimentFromExpressionSet(modencodefly.eset)
ddsSE <- DESeqDataSet(summ_exp, design= ~ stage)

ddsSE <- DESeq(ddsSE)
resSE <- results(ddsSE, contrast=c("stage","L2Larvae","L1Larvae"))
