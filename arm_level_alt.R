###########################################################
#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.
###########################################################
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
source("functions.R")

###########################################################
#getChrTable() 
#Reads in the Oncoscan.na33.r1.chromStats.tsv file
###########################################################

chr_table <- getChrTable()

###########################################################
#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.
###########################################################

chr_alt_summary <- function(in_file,segment_thr,min_gap)
{
  
  ###########################################################
  # Get the sample ID from the user input file
  
  #Temporary inputs
  #sample_id<-"test"
  #in_file<-"test.txt"
  #segment_thr<-1000
  #min_gap<-2000
  
  sample_id<-gsub("_gene_list_full_location.txt","",basename(in_file))
  filter<-as.numeric(segment_thr)
  smooth<-as.numeric(min_gap)
  ###############
  # getSegments
  # list of GenomicRanges with one extra column for the copy number and another for the LOH.
  ###############
  
  segments<-getSegments(in_file, chr_table)
  
  ###############
  # segments_alt()
  # Segment list based on Alteration
  ###############
  
  gain_GR<-c(segments_alt(segments, GAIN), segments_alt(segments, AMP))
  ampli_GR<-segments_alt(segments, AMP)
  loss_GR<-segments_alt(segments, LOSS)
  loh_GR<-segments_alt(segments, LOH)
  
  #################
  #trim(segs, x)
  #################
  # takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.
  
  gain_trim <- trim (gain_GR, segment_thr)
  ampli_trim <- trim (ampli_GR, segment_thr)
  loss_trim <- trim (loss_GR, segment_thr)
  loh_trim <- trim (loh_GR, segment_thr)
  
  #################
  # Function: smoothing(segs, x)
  #################
  # join segements separated by min_gap after trimming
  
  gain_smooth <- smoothing (gain_trim, min_gap)
  ampli_smooth <- smoothing (ampli_trim, min_gap)
  loss_smooth <- smoothing (loss_trim, min_gap)
  loh_smooth <- smoothing (loh_trim, min_gap)
  
  ##################
  # Function: longest(segs)
  ##################
  # find the longest segment after smoothing
  
  gain_longest <-longest(gain_smooth)
  ampli_longest <- longest(ampli_smooth)
  loss_longest <- longest(loss_smooth)
  loh_longest <- longest(loh_smooth)
  
  ###############
  # sum(segs)
  ###############
  # takes as input a list of segments and returns the sum of the length of all segments. 
  
  gain_GR_sum<-sum_seg(gain_smooth)
  ampli_GR_sum<-sum_seg(ampli_smooth)
  loss_GR_sum<-sum_seg(loss_smooth)
  loh_GR_sum<-sum_seg(loh_smooth)
  
  ##################
  # takes as input a list of segments and returns the sum of the longest segment 
  
  gain_GR_sum_long<-sum_seg(gain_longest)
  ampli_GR_sum_long<-sum_seg(ampli_longest)
  loss_GR_sum_long<-sum_seg(loss_longest)
  loh_GR_sum_long<-sum_seg(loh_longest)
  
  #################
  # percent_alt()
  #################
  # % Gain, Ampli, Loss and LOH calculations based on sum of segments altered
  
  
  gain_per<-percent_alt(gain_GR_sum, GAIN, chr_table)
  ampli_per<-percent_alt(ampli_GR_sum, AMP, chr_table)
  loss_per<-percent_alt(loss_GR_sum, LOSS, chr_table)
  loh_per<-percent_alt(loh_GR_sum, LOH, chr_table)
  
  # Final results as percent table of GAIN, LOSS adn LOH
  
  oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)
  
  oncoscan_summary<-round(oncoscan_summary, 2)
  
  #print (oncoscan_summary)
  
  # Write output file
  csvFileName <- file.path("./chr_alt_out",paste0(sample_id,"_", "sum_seg","_", filter,"_", smooth,".csv",sep=""))
  write.csv(oncoscan_summary, csvFileName)
  
  #################
  # % Gain, Ampli, Loss and LOH calculations based on longest segment
  #################
  
  gain_per<-percent_alt(gain_GR_sum_long, GAIN, chr_table)
  ampli_per<-percent_alt(ampli_GR_sum_long, AMP, chr_table)
  loss_per<-percent_alt(loss_GR_sum_long, LOSS, chr_table)
  loh_per<-percent_alt(loh_GR_sum_long, LOH, chr_table)
  
  # Final results as percent table of GAIN, LOSS adn LOH
  
  oncoscan_summary<-cbind(gain_per, ampli_per, loss_per, loh_per)
  
  oncoscan_summary<-round(oncoscan_summary, 2)
  
  #print (oncoscan_summary)
  # Write output file
  csvFileName <- file.path("./chr_alt_out",paste0(sample_id,"_", "long_seg","_", filter,"_",smooth,".csv",sep=""))
  write.csv(oncoscan_summary, csvFileName)
  
  #################
  # segment visualization
  #################
  
  #gir2 = gain_GR_smooth[seqnames(gain_GR_smooth) == '3q']
  #plotRanges(gir2,xlim=c(93517443,197852564))
  
  ##########################
  
  seg_gain <- segment_no (gain_GR, GAIN, chr_table)
  seg_ampli <- segment_no (ampli_GR, AMP, chr_table)
  seg_loss <- segment_no (loss_GR, LOSS, chr_table)
  seg_loh <- segment_no (loh_GR, LOH, chr_table)
  
  # Final results as percent table of GAIN, LOSS adn LOH
  
  seg_summary<-cbind(seg_gain, seg_ampli, seg_loss, seg_loh)
  
  #print (oncoscan_summary)
  
  # Write output file
  csvFileName <- file.path("./chr_alt_out",paste0(sample_id,"_", "seg_num","_", filter,"_", smooth,".csv",sep=""))
  write.csv(seg_summary, csvFileName)
  
}