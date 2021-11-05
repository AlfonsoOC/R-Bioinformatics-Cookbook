##Estimating the copy number at a locus of interest

#It is often of interest to know how often a sequence occurs in a sample of interest—that is, to estimate whether, in your particular sample, a locus has been duplicated or its copy number has increased. The locus could be anything from a gene at Kbp scale or a large section of DNA at Mbp scale.
#scale. Our approach in this recipe will be to use HTS read coverage after alignment to estimate a background level of coverage and then inspect the coverage of our region of interest. The ratio of the coverage in our region of interest to the background level will give us an estimate of the copy number in the region.

## 1Load the library and get counts in windows across the genome:
library(csaw)

whole_genome <- csaw::windowCounts( 
  file.path(getwd(), "R-Bioinformatics-Cookbook","datasets", "ch2", "hg17_snps.bam"),
  bin = TRUE,
  filter = 0,
  width = 100, #100bp windows over our Chromosome 17
  param = csaw::readParam(
    minq = 20,
    dedup = TRUE,
    pe = "both"
  )
)
colnames(whole_genome) <- c("h17")

## 2Extract the data from SummarizedExperiment:
counts <- assay(whole_genome)[,1] # "counts" is a vector

## 3 Work out a low count threshold and set windows with lower counts to NA:
#we use the quantile() function to get the min_count value in the lower 10th percentile of the counts vector. The min_count value will act as a cut-off. All values in the counts vector lower than this are set to NA to remove them from the analysis—
min_count <- quantile(counts, 0.1)[[1]]
counts[counts < min_count] <- NA

## 4 Double the counts of a set of windows in the middle—these will act as our high copy number region:
#we add some regions with doubled coverage—so that we can detect them. We select a number of windows to double the counts in and then create a multiplier vector of equal length to counts that contains 1 where we don't wish to change counts and 2 where we wish to double them. We then apply the multiplication. Step 4 will likely be left out in your own analysis as it is a synthetic data-generation step.
n <- length(counts)

doubled_windows <- 10

left_pad <- floor( (n/2) - doubled_windows )
right_pad <- n - left_pad -doubled_windows
multiplier <- c(rep(1, left_pad ), rep(2,doubled_windows), rep(1, right_pad) )
counts <- counts * multiplier

#  5 Calculate the mean coverage and the ratio in each window to that mean coverage, and inspect the ratio vector with a plot:
mean_cov <- mean(counts, na.rm=TRUE) 

ratio <- matrix(log2(counts / mean_cov), ncol = 1)
#We quickly use plot() to inspect ratio and can clearly see the count doubled windows in the middle of the data.
plot(ratio)

## 6 Build SummarizedExperiment with the new data and the row data of the old one:
#we build a new SummarizedExperiment object, se, to hold the window ranges and the new ratio data. We take the GRanges and colData objects from window_counts and add our new ratio matrix. We can now start to subset this and see what coverage is in our regions of interest.
se <- SummarizedExperiment(assays=list(ratio), rowRanges= rowRanges(whole_genome), colData = c("CoverageRatio"))

## 7 Create a region of interest and extract coverage data from it:
#we construct a manual GRanges object for an arbitrary region we're interested in, helpfully called region_of_interest, and use that to find the overlapping windows in our se object using findOverlaps().

region_of_interest <- GRanges(
  seqnames = c("NC_000017.10"),
  IRanges(c(40700), width = c(1500) )
)

overlap_hits <- findOverlaps(region_of_interest, se)
data_in_region <- se[subjectHits(overlap_hits)]
#In the output, we can see the region has roughly a log2 ratio of 1 (twofold) coverage relative to the background, which we can interpret as a copy number of 2.
assay(data_in_region)




##TIPS

#For bether Backgroung analisys
#The calculation for the background level in this recipe is really simple—which is great for learning the recipe, but might be quickly underpowered in your own real data. There are numerous options you could take to modify the way you calculate the background level for your own data. Check out the rollmeans() and rollmedians() functions in the zoo package—these give the mean and median in rolling windows of arbitrary step length and can give you a moving window background average that may be more appropriate.

#Other way to check SNP allele frequency
#A related analysis to copy number is the estimation of ploidy from SNP allele frequencies. You can check out the vcfR package's freq_peaks() function as a starting place to estimate ploidy from variant information in BAM files.

