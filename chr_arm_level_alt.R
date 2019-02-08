#Calculate percentage chromosome copy number alteration (Gain, Loss, LOH) based on Affymetrix Oncoscan Assay 
#Percentage alteration defined based on the chromosome length covered by oncoscan.

# Input text file of chromosome size
# Chromosome size, start and end based on Oncoscan coverage

Chromosome<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "16", "17", "18", "19", "20", "21", "Y", "X", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y", "X")

Arm<-c("p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "p", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q", "q")

Length<-c("120681421", "92250747", "90410210", "49023050", "46389046", "58565593", "57984002", "43765384", "46992414", "39021547", "51383187", "34767256", "35187838", "21841396", "15380695", "24297088", "26240161", "1583738", "7381837", "58346321", "106082854", "147623134", "104335121", "138343687", "131256346", "109026658", "98053925", "99395762", "72896070", "93031003", "80143490", "95915127", "96018327", "88016877", "82409721", "43696696", "54936486", "59460873", "31348601", "33457921", "33753073", "35159113", "15665404", "93487145")
Length<-as.numeric(Length)

Arm_str<-c("754192", "21494", "63411", "69404", "12225", "204909", "41421", "23417", "204738", "126070", "192764", "83848", "83887", "400959", "12842", "247232", "69094", "9412603", "2655180", "177942", "143130024", "95429197", "93517443", "52684891", "49441966", "61886393", "61064518", "46896972", "68170421", "42403300", "54795357", "37902988", "19084823", "19265147", "20019332", "46461309", "25326941", "18554307", "27744638", "29454542", "14344537", "16054713", "13134531", "61732219")
Arm_str<-as.numeric(Arm_str)

Arm_end<-c("121435613", "92272241", "90473621", "49092454", "46401271", "58770502", "58025423", "43788801", "47197152", "39147617", "51575951", "34851104", "35271725", "22242355", "15393537", "24544320", "26309255", "10996341", "10037017", "58524263", "249212878", "243052331", "197852564", "191028578", "180698312", "170913051", "159118443", "146292734", "141066491", "135434303", "134938847", "133818115", "115103150", "107282024", "102429053", "90158005", "80263427", "78015180", "59093239", "62912463", "48097610", "51213826", "28799935", "155219364")
Arm_end<-as.numeric(Arm_end)

chr_table<-data.frame(Chromosome, Arm, Length, Arm_str, Arm_end)

##############################################################################################################
# Read User input file
###########################################################
args <- commandArgs(trailingOnly = TRUE)

in_file<-args[1]

sample_id<-gsub(".txt","",in_file)

oncoscan <- read.table(in_file, 
               sep = "\t", header = TRUE, check.names = FALSE, na.strings=c("","NA"))
  
##########################################################
# Data exrtraction from Oncoscan file and % Gain, Loss and LOH calculations
##########################################################  
    # Variable initialization
    LOH<-NULL
    GAIN<-NULL
    LOSS<-NULL
    CHROMOSOME<-NULL
    Arm<-NULL
    start<-NULL
    end<-NULL
    
    
    for (cn in 1:nrow(chr_table))
    {
      chr_name<-chr_table$Chromosome[cn]
      chr_name<-as.vector(chr_name)
      
      chr_arm<-chr_table$Arm[cn]
      chr_arm<-as.vector(chr_arm)
      
      arm_start<-chr_table$Arm_str[cn]
      arm_start<-as.numeric(arm_start)
      arm_start<-arm_start-2               #Buffer start site 2bp 
      
      arm_end<-chr_table$Arm_end[cn]
      arm_end<-as.numeric(arm_end)
      arm_end<-arm_end+2                #Buffer end site 2bp 

      arm_length<-arm_end-arm_start
      
      ####################################################################
      # Subset input file chromosome wise
      df_chr <- oncoscan[oncoscan$Chromosome %in% chr_table$Chromosome[cn],]
      
      # Total number of predictions (rows) from oncoscan analysis for each chromosome
      total<-length(rownames(df_chr))
      
      # Defination and Initializatioin of percent loh, gain and loss vectors
      loh=0
      gain=0
      loss=0
      
      if (total >0)  # calculate % alteration if present in the oncoscan file
      {
        
      for (line in 1:total) # start of for loop over each chromosome
      {
        # For each alteration line extract location of the segment altered (segment start and end) 
        
        loc<-df_chr$`Full Location`[line]
        loc_cord<-gsub("chr.*:","",loc)
        loc_cord_list<-strsplit(loc_cord,split='-', fixed=TRUE)[[1]]
        seg_start<-loc_cord_list[1]
        seg_start<-as.numeric(seg_start)
        seg_end<-loc_cord_list[2]
        seg_end<-as.numeric(seg_end)

        # Define the length of the altered segment based on its location with respect to the chromosome arms
        
       if (chr_arm %in% "p" && seg_end <= arm_end)
       {
         size_seq=seg_end-seg_start
       }
       else if (chr_arm %in% "p" && seg_end > arm_end && seg_start < arm_end)
        {
          size_seq=arm_end-seg_start
       }
        else {size_seq=0}
        
        if (chr_arm %in% "q" && seg_start >= arm_start)
        {
          size_seq=seg_end-seg_start
        }
       else if (chr_arm %in% "q" && seg_start < arm_start && seg_end > arm_start)
        {
          size_seq=seg_end-arm_start
       }
        
        # Sum up the alterations for each chromosome
          
          if(df_chr$Type[line] == "LOH") # if matched extract and sum loh 
          {
            loh=loh+size_seq
          }
          if(df_chr$Type[line] == "Gain") # if matched extract and sum gain
          {
            gain=gain+size_seq
          }
          if(df_chr$Type[line] == "Loss") # if matched extract and sum loss
          {
            loss=loss+size_seq
          }
      } # End of for loop over each chromosome
        
      } # End of if alteration  present condition
      
      else # if no alteration is present in the oncoscan file, set them to zero
        {        
        loh=0
        gain=0
        loss=0
      }
      
      total=0 # Flush the number of alterations per chromosome
      
      
      # percent loh, gain and loss calculation
      LOH[cn]=as.integer(((loh)/arm_length)*100)

      GAIN[cn]=as.integer(((gain)/arm_length)*100)

      LOSS[cn]=as.integer(((loss)/arm_length)*100)
      
      CHROMOSOME[cn]=chr_name
      Arm[cn]=chr_arm
    }
    
    # Final results as percent table of GAIN, LOSS adn LOH
    oncoscan_summary<-data.frame(sample_id, CHROMOSOME,Arm, GAIN, LOSS, LOH)
    
    print (oncoscan_summary)
  
