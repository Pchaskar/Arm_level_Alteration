###########################################################
#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.
#Output arm level alterations
###########################################################
# Libraries
source("libraries.R")

# Functions
source("functions.R")

###########################
#getChrTable() 
###########################
#Input reference data

chr_table <- getChrTable()

#########################
# Read User input file
#########################
args <- commandArgs(trailingOnly = TRUE)

# Input file name.txt
in_file<-args[1]

# Get the sample ID from the user input file  
sample_id<-gsub("_gene_list_full_location.txt","",basename(in_file))

#######################
# Default inputs for filtering and smoothing
#######################

segment_thr<-300000 # 300k genomide resolution of oncoscan assay
min_gap<-0 # Join segments if within X distance of each other 

filter<-as.numeric(segment_thr)
smooth<-as.numeric(min_gap)

#######################
# getSegments
#######################
# list of GenomicRanges with one extra column for the copy number and another for the LOH.

segments<-getSegments(in_file, chr_table)

#######################
# segments_alt()
#######################
# Segment list based on Alteration

gain_GR<-c(segments_alt(segments, GAIN), segments_alt(segments, AMP))
ampli_GR<-segments_alt(segments, AMP)
loss_GR<-segments_alt(segments, LOSS)
loh_GR<-segments_alt(segments, LOH)

#######################
#trim(segs, x)
#######################
# takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.

gain_trim <- trim (gain_GR, segment_thr)
ampli_trim <- trim (ampli_GR, segment_thr)
loss_trim <- trim (loss_GR, segment_thr)
loh_trim <- trim (loh_GR, segment_thr)

#######################
# Function: smoothing(segs, x)
#######################
# join segements separated by min_gap after trimming

gain_smooth <- smoothing (gain_trim, min_gap)
ampli_smooth <- smoothing (ampli_trim, min_gap)
loss_smooth <- smoothing (loss_trim, min_gap)
loh_smooth <- smoothing (loh_trim, min_gap)

#######################
# sum(segs)
#######################
# takes as input a list of segments and returns the sum of the length of all segments. 

gain_GR_sum<-sum_seg(gain_smooth)
ampli_GR_sum<-sum_seg(ampli_smooth)
loss_GR_sum<-sum_seg(loss_smooth)
loh_GR_sum<-sum_seg(loh_smooth)

#######################
# percent_alt()
#######################
# % Gain, Ampli, Loss and LOH calculations based on sum of segments altered

gain_per<-percent_alt(gain_GR_sum, GAIN, chr_table)
ampli_per<-percent_alt(ampli_GR_sum, AMP, chr_table)
loss_per<-percent_alt(loss_GR_sum, LOSS, chr_table)
loh_per<-percent_alt(loh_GR_sum, LOH, chr_table)

# Final results as percent table of GAIN, LOSS adn LOH

oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)

oncoscan_summary<-round(oncoscan_summary, 2)

#print (oncoscan_summary)

#######################
#arm_level(df)
#######################
# Report Arm Level alterations based on oncoscan_summary

arm_level_scna<-arm_level(oncoscan_summary)

# Write output file with Arm Level alterations

csvFileName <- file.path(".",paste0(sample_id,"_", "arm_level","_", "scna",".csv",sep=""))
write.csv(arm_level_scna, csvFileName, row.names=FALSE)

