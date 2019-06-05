<<<<<<< HEAD
# Package
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

#################
#getChrTable() 
#Reads in the Oncoscan.na33.r1.chromStats.tsv file
#################

chr_table <- getChrTable()


#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
chr_alt_summary <- function(in_file,segment_thr,min_gap){
  

#####################
# Read User input file
#####################
#args <- commandArgs(trailingOnly = TRUE)

# Input file name.txt
#in_file<-args[1]

sample_id<-gsub("_gene_list_full_location.txt","",basename(in_file))

# Filter out not broad enough segments defined by threshold segment_thr, mega base pair range
#segment_thr<-args[2]

# minimum gap allowed between adjcent segments
#min_gap<-args[3]

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
# output: GR range object
###############

gain_GR<-c(segments_alt(segments, "Gain"), segments_alt(segments, "Ampli"))
ampli_GR<-segments_alt(segments, "Ampli")
loss_GR<-segments_alt(segments, "Loss")
loh_GR<-segments_alt(segments, "LOH")

###############
# sum_seg_all
# takes as input a list of all the segments and returns the sum of their lengths. 
###############
gain_GR_sum<-sum_seg_all(gain_GR)
ampli_GR_sum<-sum_seg_all(ampli_GR)
loss_GR_sum<-sum_seg_all(loss_GR)
loh_GR_sum<-sum_seg_all(loh_GR)

###############
# sum_seg_longest
# takes as input a list of all the segments and returns the length of longest segment. 
###############
gain_GR_longest<-sum_seg_longest(gain_GR, segment_thr, min_gap)
ampli_GR_longest<-sum_seg_longest(ampli_GR, segment_thr, min_gap)
loss_GR_longest<-sum_seg_longest(loss_GR, segment_thr, min_gap)
loh_GR_longest<-sum_seg_longest(loh_GR, segment_thr, min_gap)

#################
# % Gain, Ampli, Loss and LOH calculations based on sum of segments altered
#################

gain_per<-percent_alt(gain_GR_sum,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_sum,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_sum,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_sum,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

#print (oncoscan_summary)

# Write output file
csvFileName <- file.path("./chr_alt_out/",paste0(sample_id,"_","sum_seg",".csv",sep=""))
write.csv(oncoscan_summary, csvFileName)


#################
# % Gain, Ampli, Loss and LOH calculations based on longest segment
#################

gain_per<-percent_alt(gain_GR_longest,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_longest,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_longest,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_longest,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

#print (oncoscan_summary)
# Write output file
csvFileName <- file.path("./chr_alt_out/",paste0(sample_id,"_","long_seg",".csv",sep=""))
write.csv(oncoscan_summary, csvFileName)

#################
# segment visualization
#################

#gir2 = gain_GR_smooth[seqnames(gain_GR_smooth) == '3q']
#plotRanges(gir2,xlim=c(93517443,197852564))

##########################

}
=======
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

###########################################################
#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.
###########################################################

chr_alt_summary <- function(in_file,segment_thr,min_gap){
  
###########################################################
# Read User input file
###########################################################
#args <- commandArgs(trailingOnly = TRUE)

# Input file name.txt
#in_file<-args[1]


# Filter out not broad enough segments defined by threshold segment_thr, mega base pair range
#segment_thr<-args[2]

# minimum gap allowed between adjcent segments
#min_gap<-args[3]

##############
#Temporary inputs
#sample_id<-"test"
#in_file<-"test.txt"
#segment_thr<-1000
#min_gap<-2000
###############
sample_id<-gsub("_gene_list_full_location.txt","",basename(in_file))
###############
# getSegments
# list of GenomicRanges with one extra column for the copy number and another for the LOH.
###############

segments<-getSegments(in_file, chr_table)

###############
# segments_alt()
# Segment list based on Alteration
###############

gain_GR<-c(segments_alt(segments, "Gain"), segments_alt(segments, "Ampli"))
ampli_GR<-segments_alt(segments, "Ampli")
loss_GR<-segments_alt(segments, "Loss")
loh_GR<-segments_alt(segments, "LOH")

###############
# sum(segs)
# takes as input a list of segments and returns the sum of the length of all segments. 
###############

gain_GR_sum<-sum_seg(gain_GR)
ampli_GR_sum<-sum_seg(ampli_GR)
loss_GR_sum<-sum_seg(loss_GR)
loh_GR_sum<-sum_seg(loh_GR)

##################
#longest(segs)
# takes as input a list of segments and returns the sum of the longest segment 
##################

gain_GR_sum_long<-sum_long_seg(gain_GR, segment_thr, min_gap)
ampli_GR_sum_long<-sum_long_seg(ampli_GR, segment_thr, min_gap)
loss_GR_sum_long<-sum_long_seg(loss_GR, segment_thr, min_gap)
loh_GR_sum_long<-sum_long_seg(loh_GR, segment_thr, min_gap) 

#################
# % Gain, Ampli, Loss and LOH calculations based on sum of segments altered
#################

gain_per<-percent_alt(gain_GR_sum,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_sum,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_sum,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_sum,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

#print (oncoscan_summary)

# Write output file
csvFileName <- file.path("./chr_alt_out/",paste0(sample_id,"_","sum_seg",".csv",sep=""))
write.csv(oncoscan_summary, csvFileName)


#################
# % Gain, Ampli, Loss and LOH calculations based on longest segment
#################

gain_per<-percent_alt(gain_GR_sum_long,'GAIN', chr_table)
ampli_per<-percent_alt(ampli_GR_sum_long,'AMPLI', chr_table)
loss_per<-percent_alt(loss_GR_sum_long,'LOSS', chr_table)
loh_per<-percent_alt(loh_GR_sum_long,'LOH', chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

#print (oncoscan_summary)
# Write output file
csvFileName <- file.path("./chr_alt_out/",paste0(sample_id,"_","long_seg",".csv",sep=""))
write.csv(oncoscan_summary, csvFileName)

#################
# segment visualization
#################

#gir2 = gain_GR_smooth[seqnames(gain_GR_smooth) == '3q']
#plotRanges(gir2,xlim=c(93517443,197852564))

##########################

}
>>>>>>> 983abdba5eb490702df57fed88c808299e78c5ff
