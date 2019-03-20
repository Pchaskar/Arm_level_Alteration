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
segment_thr<-1000
min_gap<-2000

# Data extraction and transformation

###############
# getSegments
# list of GenomicRanges with one extra column for the copy number and another for the LOH.
###############

segments<-getSegments(oncoscan)

###############
# Segment list based on Alteration
###############

tup_seg_gain<-segments_alt(segments, "Gain")
tup_seg_ampli<-segments_alt(segments, "Ampli")
tup_seg_loss<-segments_alt(segments, "Loss")
tup_seg_loh<-segments_alt(segments, "LOH")

# make GR ranges from data frame

gain_GR<-makeGRangesFromDataFrame(tup_seg_gain, keep.extra.columns = TRUE)
ampli_GR<-makeGRangesFromDataFrame(tup_seg_ampli, keep.extra.columns = TRUE)
loss_GR<-makeGRangesFromDataFrame(tup_seg_loss, keep.extra.columns = TRUE)
loh_GR<-makeGRangesFromDataFrame(tup_seg_loh, keep.extra.columns = TRUE)

###############
# hasOverlaps (segs)
# takes as input a list of segments and returns true if there exists two 
# overlapping segments within the list. Returns false otherwise.
# TRUE: No overlap
# False: Has overlap
###############
olaps_gain<-hasOverlaps(gain_GR)
olaps_ampli<-hasOverlaps(ampli_GR)
olaps_loss<-hasOverlaps(loss_GR)
olaps_loh<-hasOverlaps(loh_GR)

###############
# sum(segs)
# takes as input a list of segments and returns the sum of the length of all segments. 
# Segments are expected to be non-overlapping
###############
gain_GR_sum<-sum_seg(gain_GR, olaps_gain)
ampli_GR_sum<-sum_seg(ampli_GR, olaps_ampli)
loss_GR_sum<-sum_seg(loss_GR, olaps_loss)
loh_GR_sum<-sum_seg(loh_GR, olaps_loh)

#################
# trim(segs, x)
# takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.
#################

gain_GR_filt<-trim(gain_GR_sum,segment_thr)
ampli_GR_filt<-trim(ampli_GR_sum,segment_thr)
loss_GR_filt<-trim(loss_GR_sum,segment_thr)
loh_GR_filt<-trim(loh_GR_sum,segment_thr)

#################
#smooth(segs, x)
#################
gain_GR_smooth<-smoothing(gain_GR_filt, min_gap)
ampli_GR_smooth<-smoothing(ampli_GR_filt, min_gap)
loss_GR_smooth<-smoothing(loss_GR_filt, min_gap)
loh_GR_smooth<-smoothing(loh_GR_filt, min_gap)

##################
#longest(segs)
##################

gain_GR_smooth_dist<-longest(gain_GR_smooth)
ampli_GR_smooth_dist<-longest(ampli_GR_smooth)
loss_GR_smooth_dist<-longest(loss_GR_smooth)
loh_GR_smooth_dist<-longest(loh_GR_smooth) 

# Grange to dataframe conversion
gain_df<-annoGR2DF(gain_GR_smooth_dist)
ampli_df<-annoGR2DF(ampli_GR_smooth_dist)
loss_df<-annoGR2DF(loss_GR_smooth_dist)
loh_df<-annoGR2DF(loh_GR_smooth_dist)

#################
# % Gain, Loss and LOH calculations
#################

gain_per<-percent_alt(gain_df,'GAIN')
ampli_per<-percent_alt(ampli_df,'Ampli')
loss_per<-percent_alt(loss_df,'LOSS')
loh_per<-percent_alt(loh_df,'LOH')

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

print (oncoscan_summary)

#################
# segment visualization
#################

gir2 = ampli_GR_smooth_dist[seqnames(ampli_GR_smooth_dist) == '14']
plotRanges(gir2,xlim=c(19265147,107282024))

##########################
