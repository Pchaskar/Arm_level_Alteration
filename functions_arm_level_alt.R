<<<<<<< HEAD
############################################################
# Functions Oncoscan Arm_level_Alterations
############################################################

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
    seg_cntype <- NULL # LOH, Gain, Loss or Amp
    
    if (is.na(seg_cn)){
      seg_cntype <- "LOH"
    }
    else if (seg_chr %in% c('X','Y')){
      if (seg_cn >1 && seg_cn <5){
        seg_cntype <- "Gain"
      }
      else if (seg_cn >=5){
        seg_cntype <- "Amp"
      }
      else if (seg_cn <=1){
        seg_cntype <- "Loss"
      }
    }
    else {
      if (seg_cn >2 && seg_cn <5){
        seg_cntype <- "Gain"
      }
      else if (seg_cn >=5){
        seg_cntype <- "Amp"
      }
      else if (seg_cn <2){
        seg_cntype <- "Loss"
      }
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
  
  return(do.call(c, segments_list))
}

###############
# Function: segments_alt()
###############
# Segment list based on Alteration

segments_alt  <- function(GR, Alt) {
  return(GR[GR$cntype == Alt])
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
# Function: smooth(segs, x)
#################
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


#################
# sum_seg_all
# takes as input a list of all the segments and returns the sum of their lengths. 
#################

sum_seg_all <- function(GR) {
  gr_merged <- smoothing(GR, 0)
  gr_merged$dist<-width(gr_merged)
  
  lengths <- sapply(levels(seqnames(gr_merged)), function(arm){
    ifelse(arm %in% seqnames(gr_merged), sum(gr_merged[seqnames(gr_merged) == arm]$dist), 0)
  })
  return(lengths)
}

###############
# sum_seg_longest
# takes as input a list of all the segments and returns the length of longest segment. 
###############

sum_seg_longest <- function(GR,  segment_thr, min_gap) {
  
  # trim(segs, x)
  GR_filt<-trim(GR,segment_thr)
  
  # smooth(segs, x)
  GR_smooth<-smoothing(GR_filt, min_gap)
  
  # calculate longest segement
  Gr_long<-longest(GR_smooth)

  # calculate length of the longest segment
  Gr_long$dist<-width(Gr_long)
  
  lengths <- sapply(levels(seqnames(Gr_long)), function(arm){
    ifelse(arm %in% seqnames(Gr_long), sum(Gr_long[seqnames(Gr_long) == arm]$dist), 0)
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

percent_alt <- function(alt_sum, alt, chr_table) {
  
  alt_type<-alt
  
  cnv<-data.frame()
  
  #calculate % alteration if present in the oncoscan file
  # Loop over each chromosome Arm
  
  for (i in 1:nrow(chr_table))
  {

    # chrm name, arm , start and end information based on chr_table
    chr_name<-as.character(chr_table$Names[i])
    chr_length<-as.numeric(chr_table$Arm_end[i]-chr_table$Arm_str[i])
    
    # Percentage alteration
    cnv[chr_name, alt_type]<-alt_sum[i]/chr_length*100
  }
  
  return(cnv)
  
} # End of function percent alteration

#####################################################################################################################
=======
############################################################
# Functions Oncoscan Arm_level_Alteration
############################################################

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
    seg_cntype <- NULL # LOH, Gain, Loss or Amp
    
    if (is.na(seg_cn)){
      seg_cntype <- "LOH"
    }
    else if (seg_chr %in% c('X','Y')){
      if (seg_cn >1 && seg_cn <5){
        seg_cntype <- "Gain"
      }
      else if (seg_cn >=5){
        seg_cntype <- "Amp"
      }
      else if (seg_cn <=1){
        seg_cntype <- "Loss"
      }
    }
    else {
      if (seg_cn >2 && seg_cn <5){
        seg_cntype <- "Gain"
      }
      else if (seg_cn >=5){
        seg_cntype <- "Amp"
      }
      else if (seg_cn <2){
        seg_cntype <- "Loss"
      }
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
  
  return(do.call(c, segments_list))
}

###############
# Function: segments_alt()
###############
# Segment list based on Alteration

segments_alt  <- function(GR, Alt) {
  return(GR[GR$cntype == Alt])
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
  
  # Smoothing: Reduce to join the buffered gragments
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
# Function: sum_long_seg
###############
# takes as input a list of segments and returns the sum of the logest segement only. 
# Overlapping segments are merged.

sum_long_seg <- function(GR, segment_thr, min_gap) {
  
  # takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.
  gr_filtered <- trim (GR, segment_thr)
  
  # join segements separated by min_gap
  gr_merged <- smoothing(gr_filtered, min_gap)
  
  # find the longest segment after smoothing
  gr_longest <- longest(gr_merged)
  
  # calculate segment length of longest fragment
  gr_longest$dist<-width(gr_longest)
  
  lengths <- sapply(levels(seqnames(gr_longest)), function(arm){
    ifelse(arm %in% seqnames(gr_longest), sum(gr_longest[seqnames(gr_longest) == arm]$dist), 0)
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
  
  alt_type<-alt
  
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
    cnv[chr_name, alt_type]<-sum_alt[i]/chr_length*100
  }
  
  return(cnv)
  
} # End of function percent alteration

#####################################################################################################################
>>>>>>> 983abdba5eb490702df57fed88c808299e78c5ff
