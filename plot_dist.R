# test code to validate predictions from chr_arm_level_alt.R
source("arm_level_alt.R")

cutoff <- read.table("cutoffs.dat", 
                             sep = " ", header = TRUE)

for (i in 1:nrow(cutoff))
{
segment_thr<-as.numeric(cutoff$filt[i])
min_gap<-as.numeric(cutoff$smooth[i])

# List all files in the oncoscan_files directory

file_list<-list.files(path = "onscan_files/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


#get the chr_alt_summary for all samples based on the specified cut offs
lapply(file_list, chr_alt_summary, segment_thr=segment_thr,min_gap=min_gap)
}
#####################
# Plot distribution

library(ggplot2)
library(fitdistrplus)

#################################
# List all files in the chr_alt_out directory

chr_alt_out<-list.files(path = "./results", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#################################
# readFile
#################################

readFile <- function(filename) {
  df <- read.csv(filename, header = FALSE, sep = " ")
  df <-df[,c(3,4,5,6)]
  
  colnames(df)<-c("GAIN","AMPLI","LOSS","LOH")
  
  # remove first row
  df = df[-1,]
  
  # set numeric
  df$GAIN <- as.numeric(as.character(df$GAIN))
  df$LOSS <- as.numeric(as.character(df$LOSS))
  
  # set all zero to NA
  df[df == 0] <- NA
  
#  df$CHR<-gsub("p|q", "", df$CHR)
  
return(df)
}
################################
#Q1: 
#Sum

df1<-readFile(chr_alt_out[23])
df2<-readFile(chr_alt_out[24])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='0k_0k')) + 
  geom_density(data=df2, aes(x=GAIN, col='10k_0k'))
  
ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='0k_0k')) + 
  geom_density(data=df2, aes(x=LOSS, col='10k_0k'))  
  
##################################
#Q2:
#Sum

df1<-readFile(chr_alt_out[23])
df2<-readFile(chr_alt_out[24])
df3<-readFile(chr_alt_out[25])
df4<-readFile(chr_alt_out[26])
df5<-readFile(chr_alt_out[27])
df6<-readFile(chr_alt_out[28])
df7<-readFile(chr_alt_out[29])
df8<-readFile(chr_alt_out[30])
df9<-readFile(chr_alt_out[31])
df10<-readFile(chr_alt_out[32])
df11<-readFile(chr_alt_out[33])

#GAIN
ggplot() + 
  geom_density(data=df8, aes(x=GAIN, col='300k_300k'))+
  geom_density(data=df6, aes(x=GAIN, col='300k_1M'))+
  geom_density(data=df7, aes(x=GAIN, col='300k_10M'))


#LOSS
ggplot() + 
  geom_density(data=df8, aes(x=LOSS, col='300k_300k'))+
  geom_density(data=df6, aes(x=LOSS, col='300k_1M'))+
  geom_density(data=df7, aes(x=LOSS, col='300k_10M'))
  
##############################################  
# (graph/statistic?) stating that there is no arm-bias. distribution (gains) for each arm
seg_0k_0K<-read.csv(chr_alt_out[17], header = FALSE, sep = " ")

# remove first row and first column
seg_0k_0K<-seg_0k_0K[-1,-1]

# change columns to numeric
cols.num <- c("V3", "V4", "V5", "V6")
seg_0k_0K[cols.num] <- sapply(seg_0k_0K[cols.num],as.character)
seg_0k_0K[cols.num] <- sapply(seg_0k_0K[cols.num],as.numeric)

# Remove amplifications from gains
seg_0k_0K$V3<-seg_0k_0K$V3-seg_0k_0K$V4

colnames(seg_0k_0K)<-c("Chr", "Gain", "Ampli", "Loss", "LOH")

# Change order of chromosomes
chrOrder<-c(paste(1:12,"p",sep=""), paste(16:21,"p",sep=""),"Xp","Yp",paste(1:22,"q",sep=""),"Xq","Yq")
chrOrder

seg_0k_0K$Chr<-factor(seg_0k_0K$Chr, levels=chrOrder)


# Bar plot
library(ggplot2)
library(reshape2)

# Met data frame for barchart
seg_0k_0K_melt<-melt(seg_0k_0K, id.vars="Chr")

# Plot

# legend , labels theme
My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 6),
  axis.title.y = element_text(size = 12))

#Rename columns
colnames(seg_0k_0K_melt)<-c("Chromosomes", "Alterations", "Segments")

bar1 <- ggplot(data=seg_0k_0K_melt, aes(x=Chromosomes, y=Segments, fill=Alterations))

bar1 + geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("orange", "blue", "darkgreen", "purple"),
                    name="Alterations",
                    breaks=c("Gain", "Ampli", "Loss", "LOH"),
                    labels=c("Gain", "Ampli", "Loss", "LOH")) + My_Theme

ggsave("Arm_level_dist_alterations.pdf", width = 20, height = 15, units = "cm", dpi = 300)

