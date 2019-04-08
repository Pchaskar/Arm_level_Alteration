############################################################
# Functions Oncoscan Arm_level_Alteration
############################################################

###############
# 1: getSegments
###############
# list of GenomicRanges with one extra column for the copy number and another for the LOH.


# Format data, obtain chromosome arm names e.g. p or q

getSegments  <- function(List) {
  
  # Initialize necessary variables
  seg_chr<-c()
  seg_start<-c()
  seg_end<-c()
  seg_length<-c()
  seg_cn<-c()
  seg_alt<-c()
  seg_name<-c()
  seg_arm<-c()
  mod_seg_start<-c()
  mod_seg_end<-c()
  seg_name<-c()
  
  
  for (line in 1:nrow(List))
  {
    ######## 
    #Start: Extract chr no, start, end, copy number
    
    # Full location from oncoscan file
    loc_cord<-List$`Full Location`[line]
    loc_cord_list<-strsplit(loc_cord,split=':', fixed=TRUE)[[1]]
    
    # seg chr no based on oncoscan file
    chr<-(loc_cord_list[1])
    chr<-gsub("chr","",chr)
    seg_chr[line]<-as.character(chr)
    
    # seg coordinates
    coord<-(loc_cord_list[2])
    coord_list<-strsplit(coord,split='-', fixed=TRUE)[[1]]
    
    # seg start and end based on oncoscan
    seg_start[line]<-as.numeric(coord_list[1])
    seg_end[line]<-as.numeric(coord_list[2])
    
    # seg copy number and type of alteration
    # Define alteration:  copy gain (2<n<5), copy loss (n<2), LOH and amplification (n>=5)
    
    seg_cn[line]<-List$`CN State`[line]
    if (is.na(seg_cn[line]))
    {
      seg_alt[line]<-"LOH"
    }
    else if (seg_cn[line] >2 && seg_cn[line] <5)
    {
      seg_alt[line]<-"Gain"
    }
    else if (seg_cn[line] >5)
    {
      seg_alt[line]<-"Ampli"
    }
    else if (seg_cn[line] <2)
    {
      seg_alt[line]<-"Loss"
    }
    
    ########
    # End: Extract chr no, start, end, copy number
    
    # Subset chromosome table based on chromosome number
    chr_tab_sub <- chr_table[chr_table$Chromosome %in% seg_chr[line],]
    
    
    for (i in 1:nrow(chr_tab_sub))
    {
      start<-c()
      end<-c()
      arm<-c()
      length<-c()
      
      # chrm name, arm , start and end information based on chr_table
      chr_arm<-as.character(chr_tab_sub$Arm[i])
      chr_name<-as.character(chr_tab_sub$Chromosome[i])
      arm_start<-chr_tab_sub$Arm_str[i]
      arm_end<-chr_tab_sub$Arm_end[i]
      
      if (chr_arm == "p" && chr_name %in% seg_chr[line])
      {
        if (seg_end[line] <= arm_end)
        {
          length[i]<-seg_end[line]-seg_start[line]+1
          start[i]<-seg_start[line]
          end[i]<-seg_end[line]
          arm[i]<-"p"
        }
        else if (seg_end[line] > arm_end && seg_start[line] < arm_end)
        {
          length[i]<-arm_end-seg_start[line]+1
          start[i]<-seg_start[line]
          end[i]<-arm_end
          arm[i]<-"p"
        }
        else {
          length[i]<-"0"
          start[i]<-"0"
          end[i]<-"0"
          arm[i]<-"p"
        }
        
        if (length[i] > 0)
        {
          mod_seg_start[line]<-start[i]
          mod_seg_end[line]<-end[i]
          seg_length[line]<-length[i]
          seg_arm[line]<-arm[i]
        }
        
      }
      
      else if (chr_arm == "q" && chr_name %in% seg_chr[line])
      {
        if (seg_start[line] >= arm_start)
        {
          length[i]<-seg_end[line]-seg_start[line]
          start[i]<-seg_start[line]
          end[i]<-seg_end[line]
          arm[i]<-"q"
        }
        else if (seg_start[line] < arm_start && seg_end[line] > arm_start)
        {
          length[i]<-seg_end[line]-arm_start
          start[i]<-arm_start
          end[i]<-seg_end[line]
          arm[i]<-"q"
        }
        else {
          length[i]<-"0"
          start[i]<-"0"
          end[i]<-"0"
          arm[i]<-"q"
        }
      }  
      
      if (length[i] > 0)
      {
        mod_seg_start[line]<-start[i]
        mod_seg_end[line]<-end[i]
        seg_length[line]<-length[i]
        seg_arm[line]<-arm[i]
        seg_name[line]<-paste0(seg_chr[line], seg_arm[line]) #rownames are the full arm name (e.g. '12p')
      }
      
    }
    
  }
  df_tupule <- data.frame(
    seqnames=c(seg_name),
    start = c(mod_seg_start),
    stop = c(mod_seg_end),
    type=c(seg_alt),
    cn=c(seg_cn)
  )
  
  
  return (df_tupule)
}


###############
# segments_alt()
###############
# Segment list based on Alteration

segments_alt  <- function(Data, Alt) {
  seg_alt <- subset(Data , Data$type == Alt)
  return(seg_alt)
}

###############
# hasOverlaps (segs)
###############
# TRUE: No overlap
# False: Has overlap

hasOverlaps <- function(GR) {

  #findoverlaps
  hits <- findOverlaps(GR,GR, type="any")
  
  # remove the ones where an interval is being compared to itself
  hits <- hits[queryHits(hits) != subjectHits(hits)]

  ov<-lengths(hits)
  
  keep_ov<-ov==0
  
  test_ol<-all(keep_ov)
  
  return(test_ol)
}

###############
# sum_seg(segs)
###############
# takes as input a list of segments and returns the sum of the length of all segments. 
# Segments are expected to be non-overlapping

sum_seg <- function(GR, olaps) {
  
  if (olaps=="TRUE")
  {
    sum_seq<-reduce(GR, ignore.strand=T)
    sum_seq$dist<-width(sum_seq)
  }
  return(sum_seq)
 
}

#################
#trim(segs, x)
#################
# takes as input a list of segments 'segs' and returns a list of all segments larger than 'x' Kb. x as to be >=0.


trim <- function (GR,X)
  
{
  trimmed_seg<-GR[GR$dist >=X]
  
  return (trimmed_seg)
}


#################
#smooth(segs, x)
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
#longest(segs)
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

percent_alt <- function(GR, ALT) {
  
  alt_type<-ALT
  
  #Init data table
  if (alt_type %in% "GAIN")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      GAIN = rep(0, dim(chr_table)[1]))
  }
  else if (alt_type %in% "AMPLI")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      AMPLI = rep(0, dim(chr_table)[1]))
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
  
  
  # Check if GR is not empty
  
  total <-length(GR)

  if (total > 0)  # calculate % alteration if present in the oncoscan file
  {
  # Loop over each chromosome Arm
  
  for (i in 1:nrow(chr_table))
  {
    sub_gr<-c()
    sum_gr<-c()
    
    # chrm name, arm , start and end information based on chr_table
    chr_name<-as.character(chr_table$Names[i])
    chr_length<-as.numeric(chr_table$Length[i])
    
    # subset Grange
    sub_gr<-GR[seqnames(GR) %in% chr_name]
    sum_gr<-sum(width(reduce(sub_gr, ignore.strand=T)))
    
    # Percentage alteration
    cnv[chr_name, alt_type]<-sum_gr/chr_length*100
  }
  
  }
  return(cnv)
  
} # End of function percent alteration

#####################################################################################################################
