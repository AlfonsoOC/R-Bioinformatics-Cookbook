##Aligning genomic length sequences with DECIPHER

## 1 Load in the libraries and genome sequences:
library(DECIPHER)
long_seqs <- readDNAStringSet(file.path(getwd(),"GitHub/R-Bioinformatics-Cookbook/", "datasets", "ch3", "plastid_genomes.fa"))
long_seqs

## 2 Prepare the sequences in a local database:
##make a sequence db on disk
Seqs2DB(long_seqs, "XStringSet", "long_db", names(long_seqs))

## 3 Find the blocks of synteny:
## Find Syntenic bloocks
synteny <- FindSynteny("long_db")

##view syntenic blocks
pairs(synteny)
plot(synteny) ## 4 Plot the syntetic blocks


## 5 make an actual alignment
alignment <- AlignSynteny(synteny, "long_db")

## 6 And save the pairwise alignments one-by-one:
##Is a structure holding all aligned syntenic blocks for each pair of sequences here genomes 1 and 2
blocks <- unlist(alignment[[1]])

##therefore write one alignment at a time
writeXStringSet(blocks, "genome_blocks_out.fa")
```