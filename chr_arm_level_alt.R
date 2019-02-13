#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.

# Input text file of chromosome size
# Chromosome size, start and end based on Oncoscan coverage

#Do not hard code! Use the text file instead (see below)
# Chromosome<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "16", "17", "18", "19", "20", "21", "Y", "X", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y", "X")
# Arm<-c("p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q")
# Length<-c("120681421", "92250747", "90410210", "49023050", "46389046", "58565593", "57984002", "43765384", "46992414", "39021547", "51383187", "34767256", "35187838", "21841396", "15380695", "24297088", "26240161", "1583738", "7381837", "58346321", "106082854", "147623134", "104335121", "138343687", "131256346", "109026658", "98053925", "99395762", "72896070", "93031003", "80143490", "95915127", "96018327", "88016877", "82409721", "43696696", "54936486", "59460873", "31348601", "33457921", "33753073", "35159113", "15665404", "93487145")
# Length<-as.numeric(Length)
# Arm_str<-c("754192", "21494", "63411", "69404", "12225", "204909", "41421", "23417", "204738", "126070", "192764", "83848", "83887", "400959", "12842", "247232", "69094", "9412603", "2655180", "177942", "143130024", "95429197", "93517443", "52684891", "49441966", "61886393", "61064518", "46896972", "68170421", "42403300", "54795357", "37902988", "19084823", "19265147", "20019332", "46461309", "25326941", "18554307", "27744638", "29454542", "14344537", "16054713", "13134531", "61732219")
# Arm_str<-as.numeric(Arm_str)
# Arm_end<-c("121435613", "92272241", "90473621", "49092454", "46401271", "58770502", "58025423", "43788801", "47197152", "39147617", "51575951", "34851104", "35271725", "22242355", "15393537", "24544320", "26309255", "10996341", "10037017", "58524263", "249212878", "243052331", "197852564", "191028578", "180698312", "170913051", "159118443", "146292734", "141066491", "135434303", "134938847", "133818115", "115103150", "107282024", "102429053", "90158005", "80263427", "78015180", "59093239", "62912463", "48097610", "51213826", "28799935", "155219364")
# Arm_end<-as.numeric(Arm_end)
# chr_table<-data.frame(Chromosome, Arm, Length, Arm_str, Arm_end)

#Reads in the Oncoscan.na33.r1.chromStats.tsv file
chromstats <- read.table('Q:/1.BIOINFORMATIQUE/GitCentralRepo/CNV-tools/CNV-tools/data/Oncoscan.na33.r1.chromStats.tsv', header = TRUE, na.strings = 'None')
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
    cnv[arm, 'GAIN'] <- gain/arm_length
    cnv[arm, 'LOSS'] <- loss/arm_length
    cnv[arm, 'LOH'] <- loh/arm_length
    
  } # End of if alteration  present condition
  
  
}

# Final results as percent table of GAIN, LOSS adn LOH
oncoscan_summary<-data.frame(sample_id, chr_table, cnv)

print (oncoscan_summary)
  
