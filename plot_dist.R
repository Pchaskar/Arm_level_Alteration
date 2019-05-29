# test code to validate predictions from chr_arm_level_alt.R
source("chr_arm_level_alt.R")

segment_thr<-3000000
min_gap<-10000000

# List all files in the oncoscan_files directory

file_list<-list.files(path = "./oncoscan_files/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


#get the chr_alt_summary for all samples based on the specified cut offs
lapply(file_list, chr_alt_summary, segment_thr=segment_thr,min_gap=min_gap)

#####################
# Plot distribution

library(ggplot2)
library(fitdistrplus)

#################################
# List all files in the chr_alt_out directory

chr_alt_out<-list.files(path = "./chr_alt_out/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#################################
# readFile
#################################

readFile <- function(filename) {
  df <- read.csv(filename, header = FALSE, sep = ",")
  
  colnames(df)<-c("CHR","GAIN","AMPLI","LOSS","LOH")
  
  df[df == 0] <- NA
  
  df$CHR<-gsub("p|q", "", df$CHR)
  
return(df)
}
################################
#lONG

df1<-readFile(chr_alt_out[6])
df2<-readFile(chr_alt_out[7])
df3<-readFile(chr_alt_out[8])
df4<-readFile(chr_alt_out[9])
df5<-readFile(chr_alt_out[10])
df6<-readFile(chr_alt_out[11])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='300k_10MB')) + 
  geom_density(data=df2, aes(x=GAIN, col='300k_1Mb')) +
  geom_density(data=df3, aes(x=GAIN, col='300k_300K')) +
  geom_density(data=df4, aes(x=GAIN, col='3Mb_10Mb')) +
  geom_density(data=df5, aes(x=GAIN, col='3Mb_1Mb')) +
  geom_density(data=df6, aes(x=GAIN, col='3Mb_300K'))

ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='300k_10MB')) + 
  geom_density(data=df2, aes(x=LOSS, col='300k_1Mb')) +
  geom_density(data=df3, aes(x=LOSS, col='300k_300K')) +
  geom_density(data=df4, aes(x=LOSS, col='3Mb_10Mb')) +
  geom_density(data=df5, aes(x=LOSS, col='3Mb_1Mb')) +
  geom_density(data=df6, aes(x=LOSS, col='3Mb_300K'))


################################
#SUM
df1<-readFile(chr_alt_out[17])
df2<-readFile(chr_alt_out[18])
df3<-readFile(chr_alt_out[19])
df4<-readFile(chr_alt_out[20])
df5<-readFile(chr_alt_out[21])
df6<-readFile(chr_alt_out[22])


ggplot() + 
  geom_density(data=df1, aes(x=GAIN, col='300k_10MB')) + 
  geom_density(data=df2, aes(x=GAIN, col='300k_1Mb')) +
  geom_density(data=df3, aes(x=GAIN, col='300k_300K')) +
  geom_density(data=df4, aes(x=GAIN, col='3Mb_10Mb')) +
  geom_density(data=df5, aes(x=GAIN, col='3Mb_1Mb')) +
  geom_density(data=df6, aes(x=GAIN, col='3Mb_300K'))

ggplot() + 
  geom_density(data=df1, aes(x=LOSS, col='300k_10MB')) + 
  geom_density(data=df2, aes(x=LOSS, col='300k_1Mb')) +
  geom_density(data=df3, aes(x=LOSS, col='300k_300K')) +
  geom_density(data=df4, aes(x=LOSS, col='3Mb_10Mb')) +
  geom_density(data=df5, aes(x=LOSS, col='3Mb_1Mb')) +
  geom_density(data=df6, aes(x=LOSS, col='3Mb_300K'))

##########################
#f <- fitdist(data$GAIN/100, 'beta', method='mme', start=list(shape1=0.001, shape2=0.001))
#f$estimate
#ratio=log(f$estimate[1]/f$estimate[2])
#ratio
