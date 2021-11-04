##Finding SNPs and indels from sequence data using VariantTools

#A key bioinformatics task is to take an alignment of high-throughput sequence reads, typically stored in a BAM file, and compute a list of variant positions.
#we'll use a set of synthetic reads on the first 83 KB or so of the human genome chromosome 17.

##Importing libraries
library(GenomicRanges)
library(gmapR)
library(rtracklayer)
library(VariantAnnotation)
library(VariantTools)



##load the datasets:

bam_folder <- file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch2")
bam_file <- file.path( bam_folder, "hg17_snps.bam") #this is the aligment with the reference sequence or genome

fasta_file <- file.path(bam_folder,"chr17.83k.fa") #this is the reference sequence or genome

##Set up the genome object and the parameter objects:
fa <- rtracklayer::FastaFile(fasta_file)

genome <- gmapR::GmapGenome(fa, create=TRUE) #we create a GmapGenome object using the gmapR::GmapGenome() function with the fasta object—this describes the genome to the later variant-calling function.

#The next two functions we use, TallyVariantParams() and VariantCallingFilters(), are vital for the correct calling and filtering of candidate SNPs.
#These are the functions in which you can set the parameters that define an SNP or indel. The options here are deliberately very poor. As you can see from the output, there are 6 SNPs called, when we created 64.
qual_params <- TallyVariantsParam(
  genome = genome,
  minimum_mapq = 20)


var_params <- VariantCallingFilters(read.count = 19,
                                    p.lower = 0.01
)

##Call the variants:
#Once the parameters are defined, we use the callVariants() function with all of the information we set up to get a vranges object of variants.
called_variants <- callVariants(bam_file, 
                                qual_params, 
                                calling.filters = var_params
)

head(called_variants)

##Now, we move on to annotation and load in the feature position information from a .gff or .bed file:
#We can then set up the GRanges object of the GFF file of annotations (I also provided a function for getting annotations from BED files).
VariantAnnotation::sampleNames(called_variants) <- "sample_name"
vcf <- VariantAnnotation::asVCF(called_variants)
VariantAnnotation::writeVcf(vcf, "hg17.vcf")

get_annotated_regions_from_gff <- function(file_name) {
  gff <- rtracklayer::import.gff(file_name) 
  as(gff, "GRanges")
}

get_annotated_regions_from_bed <- function(file_name){
  bed <- rtracklayer::import.bed(file_name)
  as(bed, "GRanges")
}

genes <- get_annotated_regions_from_gff(file.path( bam_folder, "chr17.83k.gff3"))


##Now we calculate which variants overlap which genes:
#The final step is to use the powerful overlapping and subsetting capability of the XRanges objects. We use GenomicRanges::findOverlaps() to find the actual overlap—the returned overlaps object actually contains the indices in each input object of the overlapped object.
overlaps <- GenomicRanges::findOverlaps(called_variants, genes)
overlaps

##Finally, we subset the genes with the list of overlaps.
#Hence, we can use subjectHits(overlaps) to directly subset the genes with SNPs inside and get a very non-redundant list.
genes[subjectHits(overlaps)]



#When we're happy with the filters and the set of variants we called, we can save a VCF file of the variants using the following code:
VariantAnnotation::sampleNames(called_variants) <- "sample_name"
vcf <- VariantAnnotation::asVCF(called_variants)
VariantAnnotation::writeVcf(vcf, "hg17.vcf")

