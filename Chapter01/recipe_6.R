#Differential peak analysis

#When you've discovered unannotated transcripts you may want to see whether they are differentially expressed between experiments. We've already looked at how we might do that with 
#edgeR and DESeq, but one problem is going from an object such as a RangedSummarizedExperiment, comprised of the data and a GRanges object that describes the peak regions, 
#to the internal DESeq object. In this recipe, we'll look at how we can summarise the data in those objects and get them into the correct format.


library(SummarizedExperiment)

#Load data and set up a function that creates region tags:
arab_rse <- readRDS(file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "arabidopsis_rse.RDS") ) #RDS and RDATA R Data File Format, in this case is a GRanges Object

#this is a fuction that make the rownames as follows: cromosoma:start-end
make_tag <- function(grange_obj){ 
    paste0(
        grange_obj@seqnames, 
        ":", 
        grange_obj@ranges@start, 
        "-", 
        (grange_obj@ranges@start + grange_obj@ranges@width)
    ) 
}

#Extract data and annotate rows:
#the code pulls out the actual data from RangedSummarizedExperiment using the assay() function.
counts <- assay(arab_rse)

#then we use our function to give rownames to our data 
if ( ! is.null(names(rowRanges(arab_rse))) ){
  rownames(counts) <- names(rowRanges(arab_rse))
} else {
  rownames(counts) <- make_tag(rowRanges(arab_rse))
}

head(counts) #This will give the following output—a count matrix that has the location tag as the row name that can be used in DESeq and edgeR as described in Recipes 1 and 2 in this chapter:


