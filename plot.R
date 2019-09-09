# test code to validate predictions from chr_arm_level_alt.R
source("script2.R")


file_list<-list.files(path = "test_onscan_files", pattern = NULL, all.files = FALSE,
                      full.names = TRUE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

#get the chr_alt_summary for all samples based on the specified cut offs
lapply(file_list, chr_alt_summary)

