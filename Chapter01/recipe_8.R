##Finding allele-specific expressions with AllelicImbalance
#An allele-specific expression is a situation that occurs when there is a differential abundance of different allelic variants of a transcript.

##Load libraries and set up an import folder:

library(AllelicImbalance)
library(VariantAnnotation)

region_of_interest <- GRanges(seqnames = c("17"), ranges = IRanges(79478301, 79478361)) #create a GRanges Objectof our region of interest  .

bam_folder <- file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "allele_expression") #the folder holding the .bam files of reads

##Load reads and variants in regions of interest:
reads <- impBamGAL(bam_folder, region_of_interest, verbose = FALSE) #2, the impBamGAL() function loads in reads in the region of interest.

vcf_file <-file.path( getwd(), "R-Bioinformatics-Cookbook","datasets", "ch1", "allele_expression","ERP000101.vcf" )#The variant information is loaded into variant_positions—
variant_positions <- granges(VariantAnnotation::readVcf(vcf_file), "hg19" )

allele_counts <- getAlleleCounts(reads, variant_positions, verbose=FALSE) #another GRanges object and the reads and variants are used to make allele counts with getAlleleCounts().

##Build the ASEset object:
ase.vcf <- ASEsetFromCountList(rowRanges = variant_positions, allele_counts) #we can build the ASESet object, ase.vcf (a class that inherits from RangedSummarizedExperiment), using the constructor function, ASEsetFromCountList();

reference_sequence <- file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "allele_expression", "hg19.chr17.subset.fa") #we then use the setter functions, ref() and alt(), to apply the reference and alternative base identities.

ref(ase.vcf) <- refAllele(ase.vcf,fasta=reference_sequence)
alt(ase.vcf) <- inferAltAllele(ase.vcf)


##Run the test on all variants:
#apply tests. binom.test() carries out binomial per position per sample (.bam file) tests for deviations from equality in counts in reference and alternative alleles.
binom.test(ase.vcf, n="*") #The parameter n tells the test which strand to consider—in this example, we haven't set up per-strand information, so we use "*" to ignore strandedness.

#The preceding analysis can be extended to carry out per strand and per phenotype tests if required. The script would need amending to introduce strand information in the ASESet object 
#construction step. Doing so usually requires that the RNAseq experiment and alignment steps were performed with strandedness in mind and the bioinformatics pipeline up to here configured 
#accordingly. Phenotype information can be added in the construction step using the colData parameter and a vector of phenotype or sample types for columns in the ASESet object.

