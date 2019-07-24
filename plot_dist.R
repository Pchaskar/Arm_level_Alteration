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

chr_alt_out<-list.files(path = "./results/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#################################
# readFile
#################################

readFile <- function(filename) {
  df <- read.csv(filename, header = FALSE, sep = " ")
  df <-df[,c(2,3,4,5,6,1)]
  
  colnames(df)<-c("CHR","GAIN","AMPLI","LOSS","LOH", "Filename")
  
  df[df == 0] <- NA
  
#  df$CHR<-gsub("p|q", "", df$CHR)
  
return(df)
}
################################
#lONG

df1<-readFile(chr_alt_out[1])
df2<-readFile(chr_alt_out[2])
df3<-readFile(chr_alt_out[3])
df4<-readFile(chr_alt_out[4])
df5<-readFile(chr_alt_out[5])
df6<-readFile(chr_alt_out[6])
df7<-readFile(chr_alt_out[7])
df8<-readFile(chr_alt_out[8])
df9<-readFile(chr_alt_out[9])
df10<-readFile(chr_alt_out[10])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=GAIN, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=GAIN, col='10k_1K')) +
  geom_density(data=df4, aes(x=GAIN, col='10k_300k')) +
  geom_density(data=df5, aes(x=GAIN, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=GAIN, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=GAIN, col='300k_300K')) +
  geom_density(data=df8, aes(x=GAIN, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=GAIN, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=GAIN, col='3Mb_300K'))
  
  
  

ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=LOSS, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=LOSS, col='10k_1K')) +
  geom_density(data=df4, aes(x=LOSS, col='10k_300k')) +
  geom_density(data=df5, aes(x=LOSS, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=LOSS, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=LOSS, col='300k_300K')) +
  geom_density(data=df8, aes(x=LOSS, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=LOSS, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=LOSS, col='3Mb_300K'))


################################
#SEG
#lONG

df1<-readFile(chr_alt_out[11])
df2<-readFile(chr_alt_out[12])
df3<-readFile(chr_alt_out[13])
df4<-readFile(chr_alt_out[14])
df5<-readFile(chr_alt_out[15])
df6<-readFile(chr_alt_out[16])
df7<-readFile(chr_alt_out[17])
df8<-readFile(chr_alt_out[18])
df9<-readFile(chr_alt_out[19])
df10<-readFile(chr_alt_out[20])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=GAIN, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=GAIN, col='10k_1K')) +
  geom_density(data=df4, aes(x=GAIN, col='10k_300k')) +
  geom_density(data=df5, aes(x=GAIN, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=GAIN, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=GAIN, col='300k_300K')) +
  geom_density(data=df8, aes(x=GAIN, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=GAIN, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=GAIN, col='3Mb_300K'))




ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=LOSS, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=LOSS, col='10k_1K')) +
  geom_density(data=df4, aes(x=LOSS, col='10k_300k')) +
  geom_density(data=df5, aes(x=LOSS, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=LOSS, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=LOSS, col='300k_300K')) +
  geom_density(data=df8, aes(x=LOSS, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=LOSS, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=LOSS, col='3Mb_300K'))

##########################
#SUM
#lONG

df1<-readFile(chr_alt_out[21])
df2<-readFile(chr_alt_out[22])
df3<-readFile(chr_alt_out[23])
df4<-readFile(chr_alt_out[24])
df5<-readFile(chr_alt_out[25])
df6<-readFile(chr_alt_out[26])
df7<-readFile(chr_alt_out[27])
df8<-readFile(chr_alt_out[28])
df9<-readFile(chr_alt_out[29])
df10<-readFile(chr_alt_out[30])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=GAIN, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=GAIN, col='10k_1K')) +
  geom_density(data=df4, aes(x=GAIN, col='10k_300k')) +
  geom_density(data=df5, aes(x=GAIN, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=GAIN, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=GAIN, col='300k_300K')) +
  geom_density(data=df8, aes(x=GAIN, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=GAIN, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=GAIN, col='3Mb_300K'))




ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='10k_1MB')) + 
  geom_density(data=df2, aes(x=LOSS, col='10k_10Mb')) +
  geom_density(data=df3, aes(x=LOSS, col='10k_1K')) +
  geom_density(data=df4, aes(x=LOSS, col='10k_300k')) +
  geom_density(data=df5, aes(x=LOSS, col='300k_1Mb')) +
  geom_density(data=df6, aes(x=LOSS, col='300k_10Mb'))+
  geom_density(data=df7, aes(x=LOSS, col='300k_300K')) +
  geom_density(data=df8, aes(x=LOSS, col='3Mb_1Mb')) +
  geom_density(data=df9, aes(x=LOSS, col='3Mb_10Mb')) +
  geom_density(data=df10, aes(x=LOSS, col='3Mb_300K'))




#f <- fitdist(data$GAIN/100, 'beta', method='mme', start=list(shape1=0.001, shape2=0.001))
#f$estimate
#ratio=log(f$estimate[1]/f$estimate[2])
#ratio

############################################################
#################################
# List all files in the chr_alt_out directory

results_10k_300k<-list.files(path = "./tmp", pattern = NULL, all.files = FALSE,
                        full.names = TRUE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

long<-readFile(chr_alt_out[4])
sum<-readFile(chr_alt_out[24])
seg<-readFile(chr_alt_out[14])

chr_table$length<-chr_table$Arm_end-chr_table$Arm_str


# Plot number of segmets vs %alteration based on sum
boxplot(split(sum$GAIN,seg$GAIN), main="Gain: trim 10k, smooth 300k ",
        xlab="Number of segments gained", ylab="Percentage Gain")

plot(sum$GAIN,seg$GAIN,xlab="Percentage Gain", ylab="Number of segments gained")

boxplot(split(sum$LOSS,seg$LOSS), main="LOSS: trim 10k, smooth 300k ",
        xlab="Number of segments lost", ylab="Percentage Alteration")

plot(sum$LOSS,seg$LOSS,xlab="Percentage LOSS", ylab="Number of segments Lost")


# Normalize the seg lengths by arm length
for(i in 1:nrow(chr_table))
{
  for(j in 1:nrow(seg))
  {
    if (chr_table$Names[i]==seg$CHR[j])
    {
      seg$GAIN_norm[j]=seg$GAIN[j]/chr_table$length[i]
      seg$LOSS_norm[j]=seg$LOSS[j]/chr_table$length[i]
    }
  }
}

boxplot(split(sum$GAIN,seg$GAIN_norm), main="Gain: trim 10k, smooth 300k ",
        xlab="Number of segments", ylab="Percentage Alteration")

boxplot(split(sum$LOSS,seg$LOSS_norm), main="LOSS: trim 10k, smooth 300k ",
        xlab="Number of segments", ylab="Percentage Alteration")

plot(sum$GAIN,seg$GAIN_norm,xlab="Percentage Alteration", ylab="Number of segments")

plot(sum$LOSS,seg$LOSS_norm,xlab="Percentage Alteration", ylab="Number of segments")

#######################3
# Correlation sum vs long
x_g<-sum$GAIN
y_g<-long$GAIN

cor(x_g, y_g, use = "complete.obs")

plot(sum$GAIN,long$GAIN,xlab="Percentage GAIN (SUM)", ylab="Percentage GAIN (Longest)")

x_l<-sum$LOSS
y_l<-long$LOSS

cor(x_l, y_l, use = "complete.obs")

plot(sum$LOSS,long$LOSS,xlab="Percentage LOSS (SUM)", ylab="Percentage LOSS (Longest)")
