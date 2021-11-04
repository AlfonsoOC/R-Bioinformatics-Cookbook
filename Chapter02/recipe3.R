##Predicting open reading frames in long reference sequences

#A draft genome assembly of a previously unsequenced genome can be a rich source of biological knowledge, but when genomics resources such as gene annotations aren't available, it can be tricky 
#to proceed. Here, we'll look at a first stage pipeline for finding potential genes and genomic loci of interest absolutely de novo and without information beyond the sequence. sequence. 
#We'll use a very simple set of rules to find open reading frames—sequences that begin with a start codon and end with a stop codon. The tools for doing this are encapsulated within a single 
#function in the Bioconductor package, systemPipeR. We'll end up with yet another GRanges object that we can integrate into processes downstream that allow us to cross-reference other data, 
#such as RNAseq, as we saw in the Finding unannotated transcribed regions recipe of Chapter 1, Performing Quantitative RNAseq. As a final step, we'll look at how we can use a genome simulation 
#to assess which of the open reading frames are actually likely to be real and not just occurring by chance.

#1.-Load the libraries and input genome:
library(Biostrings)
library(systemPipeR)

dna_object <- readDNAStringSet(file.path(getwd(), "R-Bioinformatics-Cookbook","datasets","ch2", "arabidopsis_chloroplast.fa")) #we load in the DNA sequence as a DNAStringSet object using readDNAStringSet() from Biostrings.

##2.-Predict the ORFs (open reading frames):
#The predORF() function from systemPipeR uses this object as input and actually predicts open reading frames according to the options set. Here, we're returning all ORFs on both strands.
predicted_orfs <- predORF(dna_object, n = 'all', type = 'gr', mode='ORF', strand = 'both', longest_disjoint = TRUE)
length(predicted_orfs)
predicted_orfs
#We receive a GRanges object in return, with 2,501 open reading frames described in "predicted_orfs". This is far too many, so we need to filter out those; in particular, we can work out which are ORFs that occurred by chance from the sequence. To do this, we need to do a little simulation and that's what happens in the next section of code.


##3.-Calculate the properties of the reference genome:
#To estimate the length that random ORFs can reach, we're going to create a series of random genomes of a length equal to our input sequence and with the same base proportion and see what the longest ORF that can be predicted is. We do a few iterations of this and we get an idea of what the longest ORF occurring by chance could be. This length serves as a cut-off we can use to reject the predicted ORFs in the real sequence.
bases <- c("A", "C", "T", "G") #First, we define the bases we will use as a simple character vector.

raw_seq_string <- strsplit(as.character(dna_object), "") #Then, we get a character vector of the original DNA sequence by splitting the as.character version of dna_object.

#We use this information to work out the proportions of each base in the input sequence by first counting the number of each base (resulting in counts ),then dividing it by the sequence length, resulting in probs.
#In both these steps, we use lapply() to loop over the vector bases and the list counts and apply an anonymous function that uses these two variables to give lists of results. unlist() is used on our final list to reduce it to a simple vector.
seq_length <- width(dna_object[1])
counts <- lapply(bases, function(x) {sum(grepl(x, raw_seq_string))}  )
probs <- unlist(lapply(counts, function(base_count){signif(base_count / seq_length, 2) }))

##4.-Create a function that finds the longest ORF in a simulated genome:
#Once we have the setup done, we can build our get_longest_orf_in_random_genome() simulation function. This generates a random genome by sampling length characters from the selection in bases 
#with probabilities given in probs. The vector is paste0() into a single string and then converted into a DNAStringSet object for the predORF() function. This time, we ask for only the longest 
#ORF using n = 1 and return the length of that.
get_longest_orf_in_random_genome <- function(x,
                                             length = 1000, 
                                             probs = c(0.25, 0.25, 0.25, 0.25), 
                                             bases = c("A","C","T","G")){
  
  random_genome <- paste0(sample(bases, size = length, replace = TRUE, prob = probs), collapse = "")
  random_dna_object <- DNAStringSet(random_genome)
  names(random_dna_object) <- c("random_dna_string")
  orfs <- predORF(random_dna_object, n = 1, type = 'gr', mode='ORF', strand = 'both', longest_disjoint = TRUE)
  return(max(width(orfs)))
}

##5.-Run the function on 10 simulated genomes:
#Now, we can run the function, which we do 10 times using lapply() and the length, probs, and bases information we calculated before. unlist() turns the result into a simple vector and we extract the longest of the 10 runs with max().
random_lengths <- unlist(lapply(1:10, get_longest_orf_in_random_genome, length = seq_length, probs = probs, bases = bases))

##6.-Get the length of the longest random ORF:
#Here we take the longest random ORF generated in by the past function
longest_random_orf <- max(random_lengths)

##7.-Keep only predicted ORFs longer than the longest random ORF:
# We can use subsetting on our original predicted_orfs GRanges object to keep the ORFs longer than the ones generated by chance.
keep <- width(predicted_orfs) > longest_random_orf
orfs_to_keep <- predicted_orfs[keep]
orfs_to_keep

##8.-writing to file
extracted_orfs <- BSgenome::getSeq(dna_object, orfs_to_keep) 
names(extracted_orfs) <- paste0("orf_", 1:length(orfs_to_keep))
writeXStringSet(extracted_orfs, "saved_orfs.fa")
