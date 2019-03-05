#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.

# Input text file of chromosome size
# Chromosome size, start and end based on Oncoscan coverage

# Packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

library(GenomicRanges)
suppressMessages(library(GenomicRanges))

#Reads in the Oncoscan.na33.r1.chromStats.tsv file
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

in_file<-args[1]

sample_id<-gsub(".txt","",basename(in_file))

oncoscan <- read.table(in_file, 
               sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"), stringsAsFactors = FALSE)

# To be removed
segment_info<-subset(oncoscan, select=c("Chromosome", "Full Location", "CN State", "Type"))

###############
#getSegments (Data, Alt)
###############


getSegments  <- function(Data, Alt) {
  seg_alt <- subset(Data , Data$Type == Alt)
  seg_alt <- subset(seg_alt, select=c("Full Location"))
  return(seg_alt$`Full Location`)
}

seg_gain<-getSegments(oncoscan, "Gain")
seg_loss<-getSegments(oncoscan, "Loss")
seg_loh<-getSegments(oncoscan, "LOH")

###############
# segment tuple
###############
tuple  <- function(List) {
  
  chr<-c()
  seg_start<-c()
  seg_end<-c()
  seg_length<-c()
  
for (line in 1:length(List))
{
  loc_cord<-List[line]
  loc_cord_list<-strsplit(loc_cord,split=':', fixed=TRUE)[[1]]
  chr[line]<-(loc_cord_list[1])
  chr[line]<-gsub("chr.*:","",chr[line])
  
  coord<-(loc_cord_list[2])
  
  coord_list<-strsplit(coord,split='-', fixed=TRUE)[[1]]
  
  seg_start[line]<-as.numeric(coord_list[1])
  seg_end[line]<-as.numeric(coord_list[2])
  seg_length[line]<-seg_end[line]-seg_start[line]+1
}
  df_tupule <- data.frame(
    chrom = c(chr),
    start = c(seg_start),
    stop = c(seg_end),
    length = c(seg_length)
  )
  
  return (df_tupule)
}

tup_seg_gain<-tuple(seg_gain)
tup_seg_loss<-tuple(seg_loss)
tup_seg_loh<-tuple(seg_loh)

# make GR ranges from data frame

gain_GR<-makeGRangesFromDataFrame(tup_seg_gain, keep.extra.columns = TRUE)
loss_GR<-makeGRangesFromDataFrame(tup_seg_loss, keep.extra.columns = TRUE)
loh_GR<-makeGRangesFromDataFrame(tup_seg_loh, keep.extra.columns = TRUE)

###############
# hasOverlaps (segs)
# reduce creates the smallest merged set of unique, non-overlapping pieces that cover all the bases the original set does
###############
gain_GR_reduced<-reduce(gain_GR, ignore.strand=T)
loss_GR_reduced<-reduce(loss_GR, ignore.strand=T)
loh_GR_reduced<-reduce(loh_GR, ignore.strand=T)

#################
# Get width -length or merge fragments
################
gain_GR_reduced$width<-width(gain_GR_reduced)
loss_GR_reduced$width<-width(loss_GR_reduced)
loh_GR_reduced$width<-width(loh_GR_reduced)

#################
#Filter segments smaller than x kb
#################
# Segment defination 
segment_thr <- args[2]

gain_GR_reduced_filt<-gain_GR_reduced[gain_GR_reduced$width >=segment_thr]
loss_GR_reduced_filt<-loss_GR_reduced[loss_GR_reduced$width >=segment_thr]
loh_GR_reduced_filt<-loh_GR_reduced[loh_GR_reduced$width >=segment_thr]


# Add gap lenth information as a metadata
gain_GR_reduced_filt$gaps<-width(gaps(gain_GR_reduced_filt))
loss_GR_reduced_filt$gaps<-width(gaps(loss_GR_reduced_filt))
loh_GR_reduced_filt$gaps<-width(gaps(loh_GR_reduced_filt))

print (gain_GR_reduced_filt)

################
# Smoothing Function
################
  

##########################################################
# Data extraction from Oncoscan file and % Gain, Loss and LOH calculations
##########################################################  
#Init data table
cnv <- data.frame(row.names = rownames(chr_table), 
                  LOH = rep(0, dim(chr_table)[1]), 
                  GAIN = rep(0, dim(chr_table)[1]), 
                  LOSS = rep(0, dim(chr_table)[1]))


for (arm in rownames(chr_table))
{
  chr_name <- chr_table[arm, 'Chromosome']
  #chr_name<-as.vector(chr_name)
  
  chr_arm <- chr_table[arm, 'Arm']
  #chr_arm<-as.vector(chr_arm)
  
  arm_start <- chr_table[arm, 'Arm_str']
  #arm_start<-as.numeric(arm_start)
  arm_start <- arm_start-2               #Buffer start site 2bp 
  
  arm_end <- chr_table[arm, 'Arm_end']
  #arm_end<-as.numeric(arm_end)
  arm_end <- arm_end+2                #Buffer end site 2bp 

  #arm_length<-arm_end-arm_start #Use the length column from chr_table
  arm_length <- chr_table[arm, 'Length']+4
  
  ####################################################################
  # Subset input file chromosome wise
  df_chr <- segment_info[segment_info$Chromosome == chr_name,]
  
  # Total number of predictions (rows) from oncoscan analysis for each chromosome
  total <- length(rownames(df_chr))
  
  # Defination and Initialization of percent loh, gain and loss vectors
  loh=0
  gain=0
  loss=0
  loh_segs <- 0
  gain_segs <- 0
  loss_segs <- 0
  
  if (total > 0)  # calculate % alteration if present in the oncoscan file
  {
    for (line in 1:total) # start of for loop over each chromosome
    {
      # For each alteration line extract location of the segment altered (segment start and end) 
      loc<-df_chr$`Full Location`[line]
      loc_cord<-gsub("chr.*:","",loc)
      loc_cord_list<-strsplit(loc_cord,split='-', fixed=TRUE)[[1]]
      
      seg_start<-as.numeric(loc_cord_list[1])
      seg_end<-as.numeric(loc_cord_list[2])


      # Define the length of the altered segment based on its location with respect to the chromosome arms
      if (chr_arm == "p")
      {
        if (seg_end <= arm_end)
        {
          size_seq=seg_end-seg_start+1
        }
        else if (seg_end > arm_end && seg_start < arm_end)
        {
          size_seq=arm_end-seg_start+1
        }
        else {
          size_seq=0
        }
      }
      else if (chr_arm == "q")
      {
        if (seg_start >= arm_start)
        {
          size_seq=seg_end-seg_start
        }
        else if (seg_start < arm_start && seg_end > arm_start)
        {
          size_seq=seg_end-arm_start
        }
        else {
          size_seq=0
        }
      }
      
      # Sum up the alterations for each chromosome
      if(df_chr$Type[line] == "LOH" && size_seq >= segment_thr) # if matched extract segments 
      {
       # loh=loh+size_seq
	loh_segs <-c(loh_segs, size_seq)
      }
      else if(df_chr$Type[line] == "Gain" && size_seq >= segment_thr) # if matched extract segments
      {
       # gain=gain+size_seq
	gain_segs <-c(gain_segs, size_seq)
      }
      else if(df_chr$Type[line] == "Loss" && size_seq >= segment_thr) # if matched extract segments
      {
       # loss=loss+size_seq
	loss_segs <-c(loss_segs, size_seq)
      }
    } # End of for loop over each chromosome

    # Sort the filtered segments
      loh_segs <- sort (loh_segs)
      gain_segs <- sort (gain_segs)
      loss_segs <- sort (loss_segs)

    # Sum the filtered segements based on segmemnt length cutoff
      loh <- sum (loh_segs)
      gain <- sum (gain_segs)
      loss <- sum (loss_segs)




    # percent loh, gain and loss calculation for each arm
    cnv[arm, 'GAIN'] <- gain*100/arm_length
    cnv[arm, 'LOSS'] <- loss*100/arm_length
    cnv[arm, 'LOH'] <- loh*100/arm_length
    
  } # End of if alteration  present condition
  
  
}

# Final results as percent table of GAIN, LOSS adn LOH
oncoscan_summary<-data.frame(sample_id, cnv)

#print (oncoscan_summary)
  
