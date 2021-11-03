##Estimating batch effects using SVA
#High throughput data such as RNAseq is often complicated by technical errors that are not explicitly modeled in the experimental design and can confound the detection of differential expression.

##Load the libraries and data:
library(sva)
arab <- readRDS(file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "arabidopsis.RDS")) #this is a  experiment with three replicates of mock and three of hrcc treatment.

##Filter out rows with too few counts in some experiments:
#we create a vector of row indices that we wish to retain, which we do by testing whether the row has at least two columns with a count of over 3—this is done by using apply() 
#and an anonymous function over the rows of the count matrix.
keep <- apply(arab, 1, function(x) { length(x[x>3])>=2 } )
arab_filtered <- arab[keep,]

#Create the initial design:
groups <- as.factor(rep(c("mock", "hrcc"), each=3)) #make groups as factors

#Set up the test and null models and run SVA:
test_model <- model.matrix(~groups) #we use the groups factor in model.matrix() to create the model design we want to use.
null_model <- test_model[,1]#We also need a null model, which, in this experimental design, is equivalent to the first column, so we extract that from the test_model design object.

svar <- svaseq(arab_filtered, test_model, null_model, n.sv=1) #We can then use the key svaseq() function to estimate the surrogate variable to add to our design. We add in test_model and null_model and select the number of surrogate variables to use with n.sv, which should be one for a simple design like this.


#Extract the surrogate variables to a new design for downstream use:
design <- cbind(test_model, svar$sv)#add the surrogate variable to the design model, which we do by binding test_model and the sv column of svar (svsar$sv) together. The final design object can then be used in packages such as edgeR and DESeq2 as any other and those methods will use the surrogate variable to deal with batch effects.


