#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.


# Packages
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("GenomicRanges", version = "3.8")
#BiocManager::install("IRanges", version = "3.8")
#BiocManager::install("Repitools", version = "3.8")
#BiocManager::install("Repitools", version = "3.8")


# Libraries
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(ggplot2)

# Functions
source("functions_arm_level_alt.R")

###########################################################
#getChrTable() 
#Reads in the Oncoscan.na33.r1.chromStats.tsv file
###########################################################

chr_table <- getChrTable()

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

##############
#Temporary inputs
sample_id<-"test"
in_file<-"test.txt"
segment_thr<-1000
min_gap<-2000

###############
# getSegments
# list of GenomicRanges with one extra column for the copy number and another for the LOH.
###############

segments<-getSegments(in_file, chr_table)

###############
# segments_alt()
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

gain_GR_smooth_longest<-longest(gain_GR_smooth)
ampli_GR_smooth_longest<-longest(ampli_GR_smooth)
loss_GR_smooth_longest<-longest(loss_GR_smooth)
loh_GR_smooth_longest<-longest(loh_GR_smooth) 

#################
# % Gain, Ampli, Loss and LOH calculations based on sum of segments altered
#################

gain_per<-percent_alt(gain_GR_smooth,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_smooth,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_smooth,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_smooth,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

print (oncoscan_summary)

#################
# % Gain, Ampli, Loss and LOH calculations based on longest segment
#################

gain_per<-percent_alt(gain_GR_smooth_longest,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_smooth_longest,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_smooth_longest,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_smooth_longest,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

print (oncoscan_summary)

#################
# segment visualization
#################

#gir2 = gain_GR_smooth[seqnames(gain_GR_smooth) == '3q']
#plotRanges(gir2,xlim=c(93517443,197852564))

##########################