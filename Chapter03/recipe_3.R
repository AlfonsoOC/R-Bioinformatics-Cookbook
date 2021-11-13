##Finding InterPro domains

#InterPro is a database of predictive models, or signatures, provided by multiple protein databases. InterPro aggregates information from multiple sources to reduce redundancy in annotations and aid interpretability.
#In this recipe, we'll extend the approach we used for just PFAM domains and look at getting annotations of InterPro domains on sequences of interest. We'll start with Ensembl core databases.

# #1 Load the libraries and double-check whether our database package carries the protein data we need:

library(ensembldb)
library(EnsDb.Rnorvegicus.v79)
## check has proper data
# * We use the EnsemblDB package-specific hasProteinData() function to check whether the database has the information we need. If the output is TRUE, we're good:
hasProteinData(EnsDb.Rnorvegicus.v79)

## 2 Build a list of genes to query with—note the keytype I need here is GENEID:
listTables(EnsDb.Rnorvegicus.v79) #lloking for column with ipro accession
e <- EnsDb.Rnorvegicus.v79

# get a list of geneids 
k <- head(keys(e, keytype = "GENEID"), n = 3 )

## 3 Use the select() function to pull the relevant data:
select(e, keys = GeneIdFilter(k),
       columns = c("TXBIOTYPE", "UNIPROTID", "PROTEINID","INTERPROACCESSION"))


##There's more...
#The approach we used in this recipe works well for Ensembl core databases, but there are other non-Ensembl core databases that we might want to search; for that, there is biomaRt. biomaRt allows us to define connections to other databases we may know of.


##with biomart for non ensembl organisms
library(biomaRt)

biomart_athal <- useMart(biomart = "plants_mart", host = "plants.ensembl.org", dataset = "athaliana_eg_gene")

##use below to view which marts are available, see written notes

#listMarts(host="plants.ensembl.org")
#listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))
#get data from the biomart, first arg is what to retrieve, filters is what to filter on, values is values to filter on mart is connection objet

getBM( c("tair_locus", "interpro"), filters=c("tair_locus"), values = c("AT5G40950", "AT2G40510"), mart = biomart_athal)

##See also...

#If you're having trouble working out the names of marts and columns, try the listMarts() and listDatasets() functions from bioMart, which will provide lists of currently available marts and the data they contain.
listMarts()
listDatasets(mart = ,verbose = T)
