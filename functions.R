############################################################
# Functions Oncoscan Arm_level_Alteration
############################################################

#Definitions for alteration names
LOSS <- 'Loss'
GAIN <- 'Gain'
AMP <- 'Amp'
LOH <- 'LOH'


###############
# Function: getChrTable
# Reads in the Oncoscan.na33.r1.chromStats.tsv file
###############
getChrTable <- function(){
  
  chromstats <- read.table("OncoScan.na33.r1.chromStats.tsv", header = TRUE, na.strings = 'None', stringsAsFactors = FALSE)
  chr_table <- data.frame(Chromosome = c(chromstats$Chrom, chromstats$Chrom), 
                          Arm = c(rep('p', dim(chromstats)[1]), rep('q', dim(chromstats)[1])),
                          Arm_str = as.numeric(c(chromstats$P_Start, chromstats$Q_Start)),
                          Arm_end = as.numeric(c(chromstats$P_End, chromstats$Q_End)))
  rownames(chr_table) <- paste0(chr_table$Chromosome, chr_table$Arm) #rownames are the full arm name (e.g. '12p')
  chr_table$Names<-paste0(chr_table$Chromosome, chr_table$Arm) #rownames are the full arm name (e.g. '12p')
  
  #Remove arms not covered in Oncoscan
  chr_table <- chr_table[!is.na(chr_table$Arm_str),]
  return(chr_table)
}

###############
# Function: getSegments
# Read User input file
# list of GenomicRanges with one extra-
# column for the copy number and another for the LOH.
###############

getSegments  <- function(in_filename, chr_table) {
  #Load the oncoscan file
  oncoscan_table <- read.table(in_filename, 
                               sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"), stringsAsFactors = FALSE)
  segments_list <- vector(mode = "list", length = 2*dim(oncoscan_table)[1]) #At most we will end up with twice as much segments
  
  counter <- 0
  for (i in 1:dim(oncoscan_table)[1]){
    counter <- counter+1
    ######## 
    #Start: Extract chr no, start, end, copy number
    
    # Full location from oncoscan file
    loc_cord<-oncoscan_table[i, 'Full Location']
    loc_cord_list<-strsplit(loc_cord,split=':', fixed=TRUE)[[1]]
    
    # seg chr no based on oncoscan file
    seg_chr<- loc_cord_list[1]
    seg_chr<-gsub("chr","",seg_chr) #remove 'chr' if present
    
    # seg coordinates
    coord_list<-strsplit(loc_cord_list[2],split='-', fixed=TRUE)[[1]]
    seg_start<-as.numeric(coord_list[1])
    seg_end<-as.numeric(coord_list[2])
    
    # Define alteration:  copy gain (2<n<5), copy loss (n<2), LOH and amplification (n>=5)
    seg_cn <- oncoscan_table[i, 'CN State']
    seg_cntype <- oncoscan_table[i, 'Type']
    
    if (seg_cntype == GAIN && seg_cn >= 5){
      seg_cntype <- AMP
    }
    
    #Does the chromosome has a p arm in Oncoscan
    if (paste0(seg_chr, 'p') %in% rownames(chr_table)){
      p_arm <- chr_table[chr_table$Chromosome == seg_chr & chr_table$Arm == 'p', c('Arm_str','Arm_end')]
      q_arm <- chr_table[chr_table$Chromosome == seg_chr & chr_table$Arm == 'q', c('Arm_str','Arm_end')]
      
      if (seg_start < p_arm[1]-10 | seg_end > q_arm[2]+10){
        stop('Some segments go beyond the oncoscan coverage!')
      }
      
      if (seg_start >= p_arm[1]-10 & seg_start <= p_arm[2]){ #Segment start is in the p arm
        if (seg_end >= p_arm[1]-10 & seg_end < q_arm[1]){ #Segment end is in the p arm or in the centromere
          seg <- GRanges(seqnames = factor(paste0(seg_chr, 'p'), levels = rownames(chr_table)), 
                         ranges = IRanges(start = seg_start, end = seg_end),
                         cntype = seg_cntype,
                         cn = seg_cn)
          segments_list[[counter]] <- seg
          #append(segments, seg)
        } 
        else { #Segment end is in the q arm
          seg_p <- GRanges(seqnames = factor(paste0(seg_chr, 'p'), levels = rownames(chr_table)), 
                           ranges = IRanges(start = seg_start, end = as.numeric(p_arm[2])),
                           cntype = seg_cntype,
                           cn = seg_cn)
          seg_q <- GRanges(seqnames = factor(paste0(seg_chr, 'q'), levels = rownames(chr_table)), 
                           ranges = IRanges(start = as.numeric(q_arm[1]), end = seg_end),
                           cntype = seg_cntype,
                           cn = seg_cn)
          segments_list[[counter]] <- seg_p
          counter <- counter+1
          segments_list[[counter]] <- seg_q
          #append(segments, c(seg_p, seg_q))
        }
      }
      else { #The segment start is in the centromere or the q arm
        seg <- GRanges(seqnames = factor(paste0(seg_chr, 'q'), levels = rownames(chr_table)), 
                       ranges = IRanges(start = seg_start, end = seg_end),
                       cntype = seg_cntype,
                       cn = seg_cn)
        segments_list[[counter]] <- seg
        #append(segments, seg)
      }
    }
    else { #The chromosome has no p arm so it has to be the q arm
      seg <- GRanges(seqnames = factor(paste0(seg_chr, 'q'), levels = rownames(chr_table)), 
                     ranges = IRanges(start = seg_start, end = seg_end),
                     cntype = seg_cntype,
                     cn = seg_cn)
      segments_list[[counter]] <- seg
      #append(segments, seg)
    }
  }
  
  return(do.call("c", segments_list))
}

###############
# Function: segments_alt()
###############
# Segment list based on Alteration

segments_alt  <- function(GR, Alt) {
  if (Alt %in% c(GAIN, AMP, LOSS, LOH)){
    return(GR[GR$cntype == Alt])
  }
  else {
    stop(paste(Alt, 'is not a valid alteration.'))
  }
}

#################
#trim(segs, x)
#################
# takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.

trim <- function (GR,X)
  
{
  # calculate distance
  GR$dist <- width(GR)
  
  trimmed_seg <- GR[GR$dist > X]
  
  return (trimmed_seg)
}

#################
# Function: smoothing(segs, x)
#################

smoothing  <- function(GR, gap) {
  
  # Add buffer region upstream and downstream of GRranges
  start(GR) <- start(GR) - gap
  end(GR) <- end(GR) + gap
  
  # Smoothing: Reduce to join the buffered fragments
  GR_smooth<-reduce(GR, ignore.strand=T)
  
  # Remove the added buffer after smoothing
  
  start(GR_smooth) <- start(GR_smooth) + gap
  end(GR_smooth) <- end(GR_smooth) - gap
  
  return(GR_smooth)
}

##################
# Function: longest(segs)
# Warning: If there are overlapping segments, it will not merge the segments!
##################

longest <- function (GR) {
  
  # calculate distance
  GR$dist<-width(GR)
  
  # Sort based on distance
  GR<-sort(GR,  decreasing=TRUE, by = ~dist)
  
  # total sequences
  total<-unique(as.character(seqnames(GR)))
  
  # Define empty Grange
  Gr_long <- GRanges()
  
  for (i in 1:length(total))
  {
    sub_gr<-GR[seqnames(GR) == total[i]]
    sub_long<-head(sub_gr, n=1)
    Gr_long <- append(Gr_long, sub_long)
  }
  
  return(Gr_long)
}

###############
# Function: sum_seg(segs)
###############
# takes as input a list of segments and returns the sum of the length of all segments for all arms. 
# Overlapping segments are merged.

sum_seg <- function(GR) {
  gr_merged <- smoothing(GR, 0)
  gr_merged$dist<-width(gr_merged)
  
  lengths <- sapply(levels(seqnames(gr_merged)), function(arm){
    ifelse(arm %in% seqnames(gr_merged), sum(gr_merged[seqnames(gr_merged) == arm]$dist), 0)
  })
  return(lengths)
}

###############
# Function: Plot ranges 
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
# Function: Percent Alteration
##################

percent_alt <- function(sum_alt, alt, chr_table) {
  if (!(alt %in% c(GAIN, AMP, LOSS, LOH))){
    stop(paste(alt, 'is not a valid alteration.'))
  }

  #Init data table
  cnv <- data.frame()
  
  # Check if GR is not empty
  
  # Loop over each chromosome Arm
  for (i in 1:nrow(chr_table))
  {
    # chrm name, arm , start and end information based on chr_table
    chr_name<-as.character(chr_table$Names[i])
    chr_length<-as.numeric(chr_table$Arm_end[i]-chr_table$Arm_str[i])
    
    # Percentage alteration
    cnv[chr_name, alt]<-sum_alt[i]/chr_length*100
  }
  
  return(cnv)
  
} # End of function percent alteration


##################
# Function: # of segments altered per chromosome arm
##################

segment_no <- function(GR, alt, chr_table) {
  if (!(alt %in% c(GAIN, AMP, LOSS, LOH))){
    stop(paste(alt, 'is not a valid alteration.'))
  }

  #Init data table
  cnv <- data.frame()
  
  # Check if GR is not empty
  
  # Loop over each chromosome Arm
  for (i in 1:nrow(chr_table))
  {
    
    # chrm name, arm , start and end information based on chr_table
    chr_name<-as.character(chr_table$Names[i])
    chr_length<-as.numeric(chr_table$Arm_end[i]-chr_table$Arm_str[i])
    
    # of altered segments
    cnv[chr_name, alt]<-length( GR[seqnames(GR) == chr_name])
  }
  
  return(cnv)
  
} # End of function percent alteration



