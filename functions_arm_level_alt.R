############################################################
# Functions Oncoscan Arm_level_Alteration
############################################################

###############
# Format data, obtain chromosome arm names e.g. p or q
###############
getSegments  <- function(List) {
  
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
    # Extract chr no, start and end
    
    # Full location from oncoscan file
    loc_cord<-List$`Full Location`[line]
    loc_cord_list<-strsplit(loc_cord,split=':', fixed=TRUE)[[1]]
    
    # seg chr no
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
    #seg_alt[line]<-List$Type[line]
    
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
    chrom = c(seg_chr),
    seg=c(seg_name),
    start = c(mod_seg_start),
    stop = c(mod_seg_end),
#    length = c(seg_length),
#    chrom = c(seg_chr),
#    arm = c(seg_arm),
    type=c(seg_alt),
    cn=c(seg_cn)
  )
  
  
  return (df_tupule)
}

###############
# getSegments 
###############

segments_alt  <- function(Data, Alt) {
  seg_alt <- subset(Data , Data$type == Alt)
  return(seg_alt)
}

###############
# hasOverlaps (segs)
###############

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
# sum(segs)
###############
sum <- function(GR, olaps) {
  
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
  else if (alt_type %in% "Ampli")
  {
    cnv <- data.frame(row.names = rownames(chr_table), 
                      Ampli = rep(0, dim(chr_table)[1]))
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
    ampli=0
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

#####################################################################################################################
