##Finding phenotype and genotype associations with GWAS

#A powerful application of being able to find many thousands of genetic variants in many samples using high-throughput sequencing is genome-wide association studies (GWAS) of genotype and phenotypes.
#GWAS is a genomic analysis set of genetic variants in different individuals or genetic lines to see whether any particular variant is associated with a trait. There are numerous techniques for doing this, but all rely on gathering data on variants in particular samples and working out each sample's genotype before cross-referencing with the phenotype in some way or other.


## 1 Load in the libraries and get the VCF file:
library(VariantAnnotation)
library(rrBLUP)
set.seed(1234)

vcf_file <- file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch2", "small_sample.vcf")
vcf <- readVcf(vcf_file, "hg19")

## 2 Extract the genotype, sample, and marker position information:
gts <- geno(vcf)$GT


samples <- samples(header(vcf))
markers <- rownames(gts)
chrom <- as.character(seqnames(rowRanges(vcf)))
pos <- as.numeric(start(rowRanges(vcf)))

## 3 Create a custom function to convert VCF genotypes into the convention used by the GWAS function:
#Briefly, in VCF, "0/0" means AA (homozygous), which is encoded as 1 in GWAS(), "0/1" and "1/0" is heterozygous Aa or 0 in GWAS(), and "1/1" is homozygous aa or -1 in GWAS().
convert <- function(v){
  v <- gsub("0/0", 1, v)
  v <- gsub("0/1", 0, v)
  v <- gsub("1/0", 0, v)
  v <- gsub("1/1",-1, v)
  return(v)
}

## 4 Call the function and convert the result into a numeric matrix:
gt_char<- apply(gts, convert, MARGIN = 2)

genotype_matrix <- matrix(as.numeric(gt_char), nrow(gt_char) )
colnames(genotype_matrix)<- samples

## 5 Build a dataframe describing the variant:
variant_info <- data.frame(marker = markers,
                           chrom = chrom,
                           pos = pos)

## 6 Build a combined variant/genotype dataframe:
genotypes <-  cbind(variant_info, as.data.frame(genotype_matrix))
genotypes

## 7 Build a phenotype dataframe:
phenotypes <- data.frame(
  line = samples,
  score = rnorm(length(samples))
)

phenotypes

## 8 Run GWAS:
GWAS(phenotypes, genotypes,plot=FALSE)

#GWAS(phenotypes, genotypes,plot=T)

