##Finding regions showing high expression ab initio with bumphunter

#The aim here is to group read alignments together in such a way that we will be able to mark regions that have significant coverage and then go on to compare samples for differences in expression levels.


library(Rsamtools)
library(bumphunter)

#Load data and get per-position coverage:
#Rsamtools::pileup function hel us to get a per-base coverage dataframe.
pileup_df <- Rsamtools::pileup(file.path(getwd(),"R-Bioinformatics-Cookbook", "datasets", "ch1", "windows.bam"))

#Find preliminary clusters:
clusters <- bumphunter::clusterMaker(pileup_df$seqnames, pileup_df$pos, maxGap = 100) # We give it the sequence names, positions, and a maximum distance parameter (maxGap). The function returns a vector of cluster numbers of equal length to the dataframe, indicating the cluster membership of each row in the dataframe.
table(clusters)

#Find the bumps with a minimum cutoff:
bumphunter::regionFinder(pileup_df$count, pileup_df$seqnames, pileup_df$pos, clusters, cutoff=1)




#If you have multiple experiments to analyze, try the bumphunter() function. This will operate over multiple data columns in a matrix and perform linear modeling to assess uncertainty 
#about the position and existence from the replicates; it is very similar to regionFinder() in operation.

