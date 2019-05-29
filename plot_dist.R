# test code to validate predictions from chr_arm_level_alt.R
source("chr_arm_level_alt.R")

segment_thr<-3000000
min_gap<-10000000

# List all files in the oncoscan_files directory

file_list<-list.files(path = "./oncoscan_files/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

lapply(file_list, chr_alt_summary, segment_thr=segment_thr,min_gap=min_gap)

#####################
# Plot distribution


library(ggplot2)
library(fitdistrplus)

data<-read.csv("./chr_alt_out/chr_alt_summary_sum_seg_3Mb_10Mb.csv", header = FALSE, sep = ",")
colnames(data)<-c("CHR","GAIN","AMPLI","LOSS","LOH")

#chr arm
data$ARM<-data$CHR

#remove arm information from CHR column
data$CHR<-gsub("p|q", "", data$CHR)


q<-ggplot(data) + geom_density(aes(x=GAIN, col='GAIN')) + geom_density(aes(x=LOSS, col='LOSS')); 
plot(q)


f <- fitdist(data$GAIN/100, 'beta', method='mme', start=list(shape1=0.001, shape2=0.001))
f$estimate

ratio=log(f$estimate[1]/f$estimate[2])


