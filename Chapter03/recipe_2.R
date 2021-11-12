##Finding protein domains with PFAM and bio3d


print("Arabidopsis")
library(org.At.tair.db)
columns(org.At.tair.db)

print("E.coli")
library(org.EcK12.eg.db)
columns(org.EcK12.eg.db)

print("Human")
library(org.Hs.eg.db)


##inspect the keytype available - different from columns in that these are the things we can actually use to access
#check if "PFAM" is available in here
keytypes(org.Hs.eg.db)

##make a vector of keys (here ensembl IDs) to work with, could come from other source, just pick 3 genes
#we can use other keys has: ALIAS, EVIDENCE, GENENAME... insetad of ENSEMBL 
k <- head(keys(org.Hs.eg.db, keytype = "ENSEMBL"), n = 3 ) # here we are selecting the first 3 genes
#We also can do this wit k:
#to check the names of the genes keys(org.Hs.eg.db, keytype = "ENSEMBL")
#then select 3 or more:
#k <- c("ENSG00000138435", "ENSG00000226651", "ENSG00000135902")

#select the data and get the annotation
result <- select(org.Hs.eg.db, keys = k, keytype="ENSEMBL", columns = c("PFAM")) #if we change "ENSEMBL in top we have to changed in here too 
result # this give us the PFAM IDs

## now get the descriptions

library(PFAM.db)

## Get the description object - other mappings are available
descriptions <- PFAMDE
##get all keys IE all PFAM IDs that exist in the database
all_ids <- mappedkeys(descriptions)

#get all the descriptions for the mappings as a data frame

id_description_mapping <- as.data.frame(descriptions[all_ids])

dplyr::left_join(result, id_description_mapping, by = c("PFAM" = "ac") ) #we join the descriptions with the PFAM IDs


##There's more...
##procedure. If it doesn't show up, then as an alternative procedure, it is possible to run PFAM and make the annotations yourself.
library(bio3d)
sequence <- read.fasta(file.path(getwd(),"GitHub/R-Bioinformatics-Cookbook/", "datasets", "ch3", "ecoli_hsp.fa") )
# run pfamseq on protein
result <- hmmer(sequence, type="hmmscan", db="pfam")
result$hit.tbl

