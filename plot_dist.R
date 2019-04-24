# test code to validate predictions from chr_arm_level_alt.R
source("chr_arm_level_alt.R")

segment_thr<-1000
min_gap<-2000

# List all files in the oncoscan_files directory

file_list<-list.files(path = "./oncoscan_files/", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

# select first 10 files
file_list<-head(file_list, 200)

lapply(file_list, chr_alt_summary, segment_thr=segment_thr,min_gap=min_gap)

#####################
# Plot distribution


library(ggplot2)
library(fitdistrplus)

data_sum_seg<-read.csv("./chr_alt_out/chr_alt_summary_sum_seg.csv", header = TRUE, sep = ",")

data_long_seg<-read.csv("./chr_alt_out/chr_alt_summary_long_seg.csv", header = TRUE, sep = ",")

dat<-data_sum_seg

dat$chr_arm <- dat$CHR


sapply(unique(dat$chr_arm), function(c) {
  q <- ggplot(data = dat[dat$chr_arm == c,]) + geom_density(aes(x=GAIN, col='gain')) + geom_density(aes(x=LOSS, col='loss')) + ggtitle(c); 
  plot(q)})

sapply(unique(dat$chr_arm), function(c) {
  q <- ggplot(data = dat[dat$chr_arm == c,]) + geom_density(aes(x=LOSS, col='loss')) + ggtitle(c); 
  plot(q)})

sapply(unique(dat$chr_arm), function(c) {
  q <- ggplot(data = dat[dat$chr_arm == c,]) + geom_density(aes(x=GAIN, col='gain')) + ggtitle(c); 
  plot(q)})

sapply(unique(dat$chr_arm), function(c) {
  f <- fitdist(dat$GAIN[dat$chr_arm==c]/100, 'beta', method='mme', start=list(shape1=0.001, shape2=0.001))
  return(c(f$estimate, ratio=log(f$estimate[1]/f$estimate[2])))
})




