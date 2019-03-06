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
suppressMessages(library(GenomicRanges))
library(Repitools)
library(ggplot2)

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
# Temporary hard coded values to be removed 
oncoscan <- read.table("test.txt", 
                       sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"), stringsAsFactors = FALSE)

sample_id<-"test"
segment_thr<-1000000
min_gap<-2000

##############################################################################################################
# Functions 
############################################################

###############
# getSegments 
###############

getSegments  <- function(Data, Alt) {
  seg_alt <- subset(Data , Data$Type == Alt)
  seg_alt <- subset(seg_alt, select=c("Full Location"))
  return(seg_alt$`Full Location`)
}

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
  chr[line]<-gsub("chr","",chr[line])

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

################
# Smoothing 
################
smoothing  <- function(GR, gap) {
 
# Add buffer region upstream and downstream of GRranges
start(GR) <- start(GR) - gap
end(GR) <- end(GR) + gap

# Smoothing: Reduce to join the buffered gragments
GR_smooth<-reduce(GR, ignore.strand=T)

# Remove the added buffer after smoothing

start(GR_smooth) <- start(GR_smooth) + gap
end(GR_smooth) <- end(GR_smooth) - gap

return(GR_smooth)
}

###############
# Plot ranges 
###############

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) 
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

##################
# % Alteration
##################

percent_alt <- function(DF, ALT) {
  
  alt_type<-ALT
  
  #Init data table
  if (alt_type %in% "GAIN")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      GAIN = rep(0, dim(chr_table)[1]))
  }
  else if (alt_type %in% "LOSS")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      LOSS = rep(0, dim(chr_table)[1]))
  }
  else if (alt_type %in% "LOH")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      LOH = rep(0, dim(chr_table)[1]))
  }  
  
  
  for (arm in rownames(chr_table))
  {
    chr_name <- chr_table[arm, 'Chromosome']
    
    chr_arm <- chr_table[arm, 'Arm']
    
    arm_start <- chr_table[arm, 'Arm_str']
    arm_start <- arm_start-2               #Buffer start site 2bp 
    
    arm_end <- chr_table[arm, 'Arm_end']
    arm_end <- arm_end+2                #Buffer end site 2bp 
    
    #arm_length<-arm_end-arm_start #Use the length column from chr_table
    arm_length <- chr_table[arm, 'Length']+4
    
    ####################################################################
    # Subset input file chromosome wise
    
    df_chr <- DF[DF$chr %in% chr_name,]
    
    # Total number of predictions (rows) from oncoscan analysis for each chromosome
    total <- length(rownames(df_chr))
    
    # Defination and Initialization of percent loh, gain and loss vectors
    loh=0
    gain=0
    loss=0
    alt_segs <- 0
    
    
    if (total > 0)  # calculate % alteration if present in the oncoscan file
    {
      for (line in 1:total) # start of for loop over each chromosome
      {
        # For each alteration line extract location of the segment altered (segment start and end) 
        
        seg_start<-as.numeric(df_chr$start[line])
        seg_end<-as.numeric(df_chr$end[line])
        
        
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
        alt_segs <-c(alt_segs, size_seq)
        
      } # End of for loop over each chromosome
      
      # Sum all the alterations
      totla_alt <- sum (alt_segs)
      
      # percent alteration calculation for each arm
      cnv[arm, alt_type] <- totla_alt*100/arm_length
      
    } # End of if alteration  present condition
    
  }
  return(cnv)
  
} # End of function percent alteration

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


