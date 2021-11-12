#Finding DNA motifs with universalmotif

##1 First, load the libraries and a motif of interest:

library(universalmotif)
library(Biostrings)


#this is the description of the motif 
motif <- read_matrix(file.path(getwd(), "GitHub","R-Bioinformatics-Cookbook","datasets", "ch3","simple_motif.txt"))
view_motifs(motif)

##2 Then, load in sequences to scan with the motif:

sequences <- readDNAStringSet(file.path(getwd(), "GitHub","R-Bioinformatics-Cookbook", "datasets", "ch3", "promoters.fa"))

## 3 Perform a scan of the sequences:
#The scan_sequences() function searches each of the sequences and reports where it finds motifs—check out the motif_hits object to see where they are.

motif_hits <- scan_sequences(motif, sequences = sequences)
head(motif_hits)

## 4 Calculate whether the motif is enriched in the sequences:
#It searches the sequences to find likely instances of motifs and counts them, performing Fisher's exact test to compare the frequencies of motifs in our set of sequences with their frequencies in an automatically generated background set.
motif_info <- enrich_motifs(motif, sequences, shuffle.k = 3, verbose = 0, RC = TRUE) #progress = FALSE
motif_info

## 5 Run MEME to find novel motifs:
#To find novel motifs, we run the external software MEME
meme_path = "~/Documents/GitHub/R-Bioinformatics-Cookbook/SET/TO/YOUR/MEME/PATH"
meme_run <- run_meme(sequences, bin = meme_path, output = "meme_out", overwrite.dir = TRUE)

#the meme function return a result out of R environment, then we have to import that result to R
motifs <- read_meme("meme_out/meme.txt")
view_motifs(motifs)

#Loading in motifs from pre-existing databases such as JASPAR and TRANSFAC is very easy with universalmotif as there are straightforward replacements for the read_matrix() function. Look at the following functions to load in motifs from various formats: read_cisbp(), read_homer(), read_jaspar(), read_matrix(), read_meme(), read_motifs(), and read_uniprobe().

