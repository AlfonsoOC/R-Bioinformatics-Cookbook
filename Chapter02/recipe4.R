##Plotting features on genetic maps with karyoploteR
#The package takes as input the familiar GRanges objects and creates detailed plots from configuration.

##1 First, we load the libraries:
library(karyoploteR)
library(GenomicRanges)

##2 Then, set up the genome object that will be the base for our karyotype:
genome_df <- data.frame(
  chr = paste0("chr", 1:5),
  start = rep(1, 5),
  end = c(34964571, 22037565, 25499034, 20862711, 31270811)
)
genome_gr <- makeGRangesFromDataFrame(genome_df)

##3 Set up the SNP positions we will draw on as markers:
snp_pos <- sample(1:1e7, 25)
snps <- data.frame(
  chr = paste0("chr", sample(1:5,25, replace=TRUE)),
  start = snp_pos,
  end = snp_pos
)
snps_gr <- makeGRangesFromDataFrame(snps)

##4 Create some labels for the markers:
snp_labels <- paste0("snp_", 1:25)

##5 Set the plot margins:
#Now we can set up our plot. First, we get the default plot parameter object from inside the package using getDefaultPlotParams(). We can modify this object to make any changes to the default settings in our plot.
#Note we have selected plot.type = 1—this is a simple plot with one data track directly above each chromosome region. We'll need to change the margin height of the data track to stop our marker labels pouring out over the top—this is done with plot.params$data1outmargin <- 600. Finally, we can draw our plot; we create the base plot object, kp, by calling plotKaryotype() and passing in the genome_gr object, plot.type, and the parameters in the modified plot.params object.
plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$data1outmargin <- 600

##6 Create the base plot and add tracks:
kp <- plotKaryotype(genome=genome_gr, plot.type = 1, plot.params = plot.params)
kpPlotMarkers(kp, snps_gr, labels = snp_labels)


### extrapanels
#The following example shows how to draw some numeric data onto a plot as a simple line.
#Here, we create a data.frame that has 100 random numbers that map into 100 windows of chromosome 4 and, as before, we create a GRanges object.
numeric_data <- data.frame(
  y = rnorm(100,mean = 1,sd = 0.5  ),
  chr = rep("chr4", 100),
  start = seq(1,20862711, 20862711/100),
  end = seq(1,20862711, 20862711/100)
)


numeric_data_gr <- makeGRangesFromDataFrame(numeric_data)

#This time, we'll have a data track above and below our chromosome—one for SNP markers and the other for the new data (note that this is plot.type = 2).
plot.params <- getDefaultPlotParams(plot.type=2)

##Set up plot margins:
#We then need to set the parameters for the plo—in particular, the margins, to stop labels and data overlapping; but after that, it's the same plot calls, this time adding a kpLines() call.
plot.params$data1outmargin <- 800
plot.params$data2outmargin <- 800
plot.params$topmargin <- 800

##Create a plot and add tracks:
kp <- plotKaryotype(genome=genome_gr, plot.type = 2, plot.params = plot.params)
kpPlotMarkers(kp, snps_gr, labels = snp_labels)
kpLines(kp, numeric_data_gr, y = numeric_data$y, data.panel=2)

#for more information use:
vignette("karyoploteR")

#A quirk of karyoploteR means that it only draws chromosomes horizontally. For vertical maps, there is also the chromPlot package in Bioconductor.

