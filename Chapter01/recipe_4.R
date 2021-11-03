
#Set up a loading function:
#We create a function that load gene information in a GFF file and convert it into a Bioconductor GRanges object using the rtracklayer package.
get_annotated_regions_from_gff <- function(file_name) {
  gff <- rtracklayer::import.gff(file_name)
  as(gff, "GRanges")
}

get_annotated_regions_from_bed <- function(file_name){ # This is a function to work with BedFiles
  bed <- rtracklayer::import.bed(file_name)
  as(bed, "GRanges")
}

#get counts in windows across whole genome
whole_genome <- csaw::windowCounts( 
    file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch1", "windows.bam"),
    bin = TRUE,
    filter = 0,
    width = 500, #windows size of 500
    param = csaw::readParam(
        minq = 20, #minimum phred quality is 20 (phredScore)
        dedup = TRUE, #remove PCR duplicates
        pe = "both" #require that both of the pairs of a read are aligned
        )
)
colnames(whole_genome) <- c("small_data") # We set the name of the single data column in whole_genome to small_data.

#make a GRanges object of the annotated features
annotated_regions <- get_annotated_regions_from_gff(file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "genes.gff"))


library(IRanges)
library(SummarizedExperiment)

#find the overlaps between the windows we made 
#we use the IRanges overlapsAny() function to check whether the window locations overlap at all with the gene regions. This give us a binary result (False, True) that we use to do the subseting
windows_in_genes <-IRanges::overlapsAny(
  SummarizedExperiment::rowRanges(whole_genome),
  annotated_regions
)

#Subset the windows into those in annotated and non-annotated regions:
annotated_window_counts <- whole_genome[windows_in_genes,]
non_annotated_window_counts <- whole_genome[ ! windows_in_genes,]

#Get the data out to a count matrix:
assay(annotated_window_counts)

#we can also use "assay(non_annotated_window_counts)" to get the non anotated regions
