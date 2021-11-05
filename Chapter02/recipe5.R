##Selecting and classifying variants with VariantAnnotation

#In pipelines where we've called variants, we'll often want to do subsequent analysis steps that 
#need further filtering or classification based on features of the individual variants, such as 
#the depth of coverage in the alternative allele. This is best done from a VCF file, and a common 
#protocol is to save a VCF of all variants from the actual calling step and then experiment with 
#filtering that. In this section, we'll look at taking an input VCF and filtering it to retain 
#variants in which the alternative allele is the major allele in the sample.

#Selecting and classifying variants with VariantAnnotation can be done using the following steps:

##1 Load lybrary
library(VariantAnnotation)

#The general outline is that we need to define two sets of filtering rules—prefilter and filter.


#C#2 reate a prefilter function:
#Our first line of code defines a is_not_microsat() function that, when passed a character string,
#uses the grepl() function to work out whether the line contains the word microsat and returns 
#TRUE if it doesn't.
is_not_microsat <- function(x){ !grepl("microsat", x, fixed = TRUE)}

##3 Load up the prefilter function into a FilterRules object:
prefilters <- FilterRules(list(microsat = is_not_microsat) )

##4 Create a filter function to keep variants where the reference allele is in less than half the reads:
#Our first line of code defines a is_not_microsat() function that, when passed a character string, uses the grepl() function to work out whether the line contains the word microsat and returns TRUE if it doesn't.
#To iterate over those elements, we use the lapply() function to apply an anonymous function that returns TRUE if the reference allele has a proportion lower than 0.5
#To iterate over those elements, we use the lapply() function to apply an anonymous function that returns TRUE if the reference allele has a proportion lower than 0.5
major_alt <- function(x){
  af <- info(x)$AF ## also geno() fixed() 
  result <- unlist(lapply(af, function(x){x[1] < 0.5}))
  return(result)
}

##5 Load the filter function into a FilterRules object:
filters <- FilterRules(list(alt_is_major = major_alt))

##6 Load the input VCF file and apply filters:
vcf_file <- file.path(getwd(), "/..","datasets", "ch2", "sample.vcf.gz")
filterVcf(vcf_file, "hg17", "filtered.vcf", prefilters = prefilters, filters = filters)


#TIPS
#In filter functions, we can take advantage of other accessor functions to get at different parts of the VCF record. There are the geno() and fixed() functions, which will return data structures describing these parts of the VCF record. You can use these to create filters in the same way we used info().

