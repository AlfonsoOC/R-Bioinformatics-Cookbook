## Recipe 1 - Estimating differential expression with edgeR

## Count table

### Load data
count_dataframe <- readr::read_tsv(file.path(getwd(), "R-Bioinformatics-Cookbook", "datasets", "ch1", "modencodefly_count_table.txt" ))
genes <- count_dataframe[['gene']]
count_dataframe[['gene']] <- NULL
count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes
pheno_data <- readr::read_table(file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch1", "modencodefly_phenodata.txt"))


### Specify experiments of interest 

experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which( pheno_data[['stage']] %in% experiments_of_interest ) #gives the numeric indices of the above experiments in the matrix

### Form grouping factor
library(magrittr)
grouping <- pheno_data[['stage']][columns_of_interest] %>% 
  forcats::as_factor()

### Form subset count data
counts_of_interest <-  count_matrix[,columns_of_interest]

### Form DGE
library(edgeR)
count_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)


### Perform differential expression 
design <- model.matrix(~ grouping) #we create the experimental design descriptor design object with the base model.matrix() function.
eset_dge <- edgeR::estimateDisp(count_dge, design) #We estimate the dispersions of each gene

#Finally, a generalized linear model is fit and the quasi-likelihood F-test is applied with the two uses of glmQLFTest(), first with the dispersal estimates, eset_dge, then with the resulting fit object.
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)

topTags(result) #The columns show the gene name, the logFC value of the gene, the F value, the P value and the False Detection Rate (FDR). Usually, the column we want to make statistical conclusions from is FDR.


##eset
# Here we are using an "ExpressionSet" object 
library(Biobase)
load(file.path(getwd(), "R-Bioinformatics-Cookbook/datasets/ch1/modencodefly_eset.RData")) #loading the eset file

experiments_of_interest <- c("L1Larvae", "L2Larvae") #check getting our interest experiments
columns_of_interest <- which(phenoData(modencodefly.eset)[['stage']] %in% experiments_of_interest ) #this give us the pointer of interest experiments 

grouping <- droplevels(phenoData(modencodefly.eset)[['stage']][columns_of_interest] ) #making the factors, but in here we use "droplevels" - phenoData(modencodefly.eset)[['stage']] -> give usall the levels in the ExpresionSet

counts_of_interest <- exprs(modencodefly.eset)[, columns_of_interest] # Extract the columns of interest of our data 

#the next steps are has the first example
eset_dge <- edgeR::DGEList(
  counts = counts_of_interest,
  group = grouping 
  )

design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(eset_dge, design)
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)
topTags(result)

