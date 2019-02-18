#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.

# Input text file of chromosome size
# Chromosome size, start and end based on Oncoscan coverage

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
  df_chr <- oncoscan[oncoscan$Chromosome == chr_name,]
  
  # Total number of predictions (rows) from oncoscan analysis for each chromosome
  total <- length(rownames(df_chr))
  
  # Defination and Initialization of percent loh, gain and loss vectors
  loh=0
  gain=0
  loss=0
  
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
      if(df_chr$Type[line] == "LOH") # if matched extract and sum loh 
      {
        loh=loh+size_seq
      }
      else if(df_chr$Type[line] == "Gain") # if matched extract and sum gain
      {
        gain=gain+size_seq
      }
      else if(df_chr$Type[line] == "Loss") # if matched extract and sum loss
      {
        loss=loss+size_seq
      }
    } # End of for loop over each chromosome
    # percent loh, gain and loss calculation for each arm
    cnv[arm, 'GAIN'] <- gain*100/arm_length
    cnv[arm, 'LOSS'] <- loss*100/arm_length
    cnv[arm, 'LOH'] <- loh*100/arm_length
    
  } # End of if alteration  present condition
  
  
}

# Final results as percent table of GAIN, LOSS adn LOH
oncoscan_summary<-data.frame(sample_id, cnv)

print (oncoscan_summary)
  
