############################################################
# Functions Oncoscan Arm_level_Alteration
############################################################

###############
# getSegments 
###############

getSegments  <- function(Data, Alt) {
  seg_alt <- subset(Data , Data$Type == Alt)
  seg_alt <- subset(seg_alt, select=c("Full Location", "CN State"))
  return(seg_alt)
}

###############
# segment tuple
###############

tuple  <- function(List) {
  
  chr<-c()
  seg_start<-c()
  seg_end<-c()
  seg_length<-c()
  seg_cn<-c()
  
  for (line in 1:nrow(List))
  {
    loc_cord<-List$`Full Location`[line]
    loc_cord_list<-strsplit(loc_cord,split=':', fixed=TRUE)[[1]]
    chr[line]<-(loc_cord_list[1])
    chr[line]<-gsub("chr","",chr[line])
    
    coord<-(loc_cord_list[2])
    
    coord_list<-strsplit(coord,split='-', fixed=TRUE)[[1]]
    
    seg_start[line]<-as.numeric(coord_list[1])
    seg_end[line]<-as.numeric(coord_list[2])
    seg_length[line]<-seg_end[line]-seg_start[line]+1
    seg_cn[line]<-List$`CN State`[line]
  }
  df_tupule <- data.frame(
    chrom = c(chr),
    start = c(seg_start),
    stop = c(seg_end),
    length = c(seg_length),
    cn=c(seg_cn)
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

#####################################################################################################################
