#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.
# Input text file of chromosome size
# Chromosome size, start and end based on Oncoscan coverage


# Packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges",  lib="/mnt/c/Users/prac/Documents/R_lib")
#biocLite("IRanges")
#BiocManager::install("Repitools", version = "3.8")


# Libraries
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(ggplot2)


# Functions
source("functions_arm_level_alt.R")

###########################################################
#Reads in the Oncoscan.na33.r1.chromStats.tsv file
###########################################################

chromstats <- read.table("OncoScan.na33.r1.chromStats.tsv", header = TRUE, na.strings = 'None', stringsAsFactors = FALSE)
chr_table <- data.frame(Chromosome = c(chromstats$Chrom, chromstats$Chrom), 
                        Arm = c(rep('p', dim(chromstats)[1]), rep('q', dim(chromstats)[1])),
                        Length = c((chromstats$P_End-chromstats$P_Start+1), (chromstats$Q_End-chromstats$Q_Start+1)),
                        Arm_str = c(chromstats$P_Start, chromstats$Q_Start),
                        Arm_end = c(chromstats$P_End, chromstats$Q_End))
rownames(chr_table) <- paste0(chr_table$Chromosome, chr_table$Arm) #rownames are the full arm name (e.g. '12p')

#Remove arms not covered in Oncoscan
chr_table <- chr_table[!is.na(chr_table$Length),]

##############################################################################################################
# Read User input file
###########################################################
args <- commandArgs(trailingOnly = TRUE)

# Input file name.txt
in_file<-args[1]

sample_id<-gsub(".txt","",basename(in_file))

oncoscan <- read.table(in_file, 
               sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"), stringsAsFactors = FALSE)

# Filter out not broad enough segments defined by threshold segment_thr, mega base pair range
segment_thr<-args[2]

# minimum gap allowed between adjcent segments
min_gap<-args[3]

##############################################################################################################
# Temporary hard coded values to be replaced by user input 
oncoscan <- read.table("test.txt", 
                       sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"), stringsAsFactors = FALSE)

sample_id<-"test"
segment_thr<-1000000
min_gap<-2000

##############################################################################################################
# Data extraction and transformation
#############################################

seg_gain<-getSegments(oncoscan, "Gain")
seg_loss<-getSegments(oncoscan, "Loss")
seg_loh<-getSegments(oncoscan, "LOH")

tup_seg_gain<-tuple(seg_gain)
tup_seg_loss<-tuple(seg_loss)
tup_seg_loh<-tuple(seg_loh)

# make GR ranges from data frame

gain_GR<-makeGRangesFromDataFrame(tup_seg_gain, keep.extra.columns = TRUE)
loss_GR<-makeGRangesFromDataFrame(tup_seg_loss, keep.extra.columns = TRUE)
loh_GR<-makeGRangesFromDataFrame(tup_seg_loh, keep.extra.columns = TRUE)

###############
# hasOverlaps (segs)
# Reduce
###############
gain_GR_reduced<-reduce(gain_GR, ignore.strand=T)
loss_GR_reduced<-reduce(loss_GR, ignore.strand=T)
loh_GR_reduced<-reduce(loh_GR, ignore.strand=T)

#################
# Get width of the merge fragments
################
gain_GR_reduced$dist<-width(gain_GR_reduced)
loss_GR_reduced$dist<-width(loss_GR_reduced)
loh_GR_reduced$dist<-width(loh_GR_reduced)

#################
#Filter segments smaller than x kb
#################

gain_GR_reduced_filt<-gain_GR_reduced[gain_GR_reduced$dist >=segment_thr]
loss_GR_reduced_filt<-loss_GR_reduced[loss_GR_reduced$dist >=segment_thr]
loh_GR_reduced_filt<-loh_GR_reduced[loh_GR_reduced$dist >=segment_thr]

#################
#Perform smoothing
#################
gain_GR_reduced_filt_smooth<-smoothing(gain_GR_reduced_filt, min_gap)
loss_GR_reduced_filt_smooth<-smoothing(loss_GR_reduced_filt, min_gap)
loh_GR_reduced_filt_smooth<-smoothing(loh_GR_reduced_filt, min_gap)

# Grange to dataframe conversion
gain_df<-annoGR2DF(gain_GR_reduced_filt_smooth)
loss_df<-annoGR2DF(loss_GR_reduced_filt_smooth)
loh_df<-annoGR2DF(loh_GR_reduced_filt_smooth)

#################
# % Gain, Loss and LOH calculations
#################

gain_per<-percent_alt(gain_df,'GAIN')
loss_per<-percent_alt(loss_df,'LOSS')
loh_per<-percent_alt(loh_df,'LOH')

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, loss_per, loh_per)

print (oncoscan_summary)

#################
# segment visualization
#################

gir2 = loss_GR[seqnames(loss_GR) == '13']
plotRanges(gir2,xlim=c(19084823,115103150))

##########################

