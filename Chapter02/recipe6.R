##Extracting information in genomic regions of interest

#Very often, you'll want to look in more detail at data that falls in a particular genomic region of interest, whether that be the SNPs and variants in a gene or the genes in a particular locus. 
#This extremely common task is handled very well by the extremely powerful GRanges and SummarizedExperiment objects, which are a little fiddly to set up but have very flexible subsetting 
#operations that make the effort well worth it. We'll look at a few ways to set up these objects and a few ways we can manipulate them to get interesting information.

#1 Load in packages and define some functions that create GRanges from common files:
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)

#The three functions we create all load in information from different file types, namely, .gff, .bed, and a tab-delimited .txt file, and return the necessary GRanges object.
get_granges_from_gff <- function(file_name) {
  gff <- rtracklayer::import.gff(file_name)
  as(gff, "GRanges")
}

get_granges_from_bed <- function(file_name){
  bed <- rtracklayer::import.bed(file_name)
  as(bed, "GRanges")
}

get_granges_from_text <- function(file_name){
  df <- readr::read_tsv(file_name, col_names = TRUE )
  GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
}

get_annotated_regions_from_gff <- function(file_name) {
  gff <- rtracklayer::import.gff(file_name)
  as(gff, "GRanges")
}

#2 Actually create some GRanges objects using those functions:
#we make use of the GFF and text functions to create two GRanges objects: gr_from_gff and gr_from_txt. These are then used in subsetting.
gr_from_gff <- get_annotated_regions_from_gff(file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch2", "arabidopsis_chr4.gff"))
gr_from_txt <- get_granges_from_text(file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch2", "arabidopsis_chr4.txt"))


##3 Extract a region by filtering on attributes; in this case—the seqnames and metadata columns:
## Extract by seqname or metadata
#we subset on feature attributes. The code finds features of type gene on chromosome 4. Note the difference in syntax between finding genes and features in Chr4. The base columns in the GRanges object—namely, seqnames, width, and start—all have accessor functions that return vectors.
genes_on_chr4 <- gr_from_gff[ gr_from_gff$type == "gene" & seqnames(gr_from_gff) %in% c("Chr4") ]


##4 Manually create a region of interest:
## By range
#we create a specific region in a custom minimal GRanges object. This contains only one region but more could be added just by putting more seqnames, start, and width in the manually specified vectors.

region_of_interest_gr <- GRanges(
  seqnames = c("Chr4"), 
  IRanges(c(10000,80000), width= c(1000,600))
)


##Use the region of interest to subset the larger object:
#we use the findOverlaps() function to get the indices of features in the gr_from_gff object that overlap the manually created region_of_interest and use those indices to subset the larger gr_from_gff object.
overlap_hits <- findOverlaps(region_of_interest_gr, genes_on_chr4)
features_in_region <- genes_on_chr4[subjectHits(overlap_hits) ]
features_in_region






#It's also possible to extract subsets of dataframes or matrices in the same way by taking advantage of GRanges that are part of other objects. In the following example, we create a matrix of 
#random data and use that to build a SummarizedExperiment object that uses a GRanges object to describe its rows:
set.seed(4321)
experiment_counts <- matrix( runif(4308 * 6, 1, 100), 4308)
sample_names <- c(rep("ctrl",3), rep("test",3) )
se <- SummarizedExperiment::SummarizedExperiment(rowRanges = gr_from_txt, assays = list(experiment_counts), colData = sample_names)
#Then, we can subset in the same way as before and get back a subset of the data as well as a subset of the ranges. The assay() function returns the actual data matrix:
overlap_hits <- findOverlaps(region_of_interest_gr, se)
data_in_region <- se[subjectHits(overlap_hits) ]
assay(data_in_region)
