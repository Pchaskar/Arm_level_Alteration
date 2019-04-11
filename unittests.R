# Unit tests for the file functions_arm_level_alt.R
# Author: Yann Christinat
# Date: 13/3/2019

library(GenomicRanges)
setwd('C:/Users/YCHR/PycharmProjects/Arm_level_Alteration/')

source("functions_arm_level_alt.R")

#
# Test getSegments
#
test_getSegments <- function(){
  chrtable <- getChrTable()
  
  fn1 <- 'testfiles/H19002057_gene_list_full_location.txt'
  fn2 <- 'testfiles/H19001711_gene_list_full_location.txt'
  fn3 <- 'testfiles/testmix_gene_list_full_location.txt'
  segs1 <- getSegments(fn1, chrtable)
  segs2 <- getSegments(fn2, chrtable)
  segs3 <- getSegments(fn3, chrtable)
}

#
# Test smoothing
#
test_smoothing <- function(){
  gr <- GRanges(
    seqnames = c('1p', '1p', '1p', '1q', '1p', '1q', '1q', '1q'),
    ranges = IRanges(start = c(1, 10, 100, 50, 30, 30, 35, 40),
                       end = c(9, 20, 150, 54, 40, 37, 42, 50)))

  s0 <- smoothing(gr, 0) #Overlapping and contiguous segments are merged 
  s0_expected <- GRanges(
    seqnames = c('1p', '1p', '1p', '1q'),
    ranges = IRanges(start = c(1, 100, 30, 30),
                       end = c(20, 150, 40, 54)))
  
  s4 <- smoothing(gr, 4) 
  #Same as s0
  
  s5 <- smoothing(gr, 5)
  s5_expected <- GRanges(
    seqnames = c('1p', '1p', '1q'),
    ranges = IRanges(start = c(1, 100, 30),
                     end = c(40, 150, 54)))
  
  s6 <- smoothing(gr, 6)
  s5_expected <- GRanges(
    seqnames = c('1p', '1p', '1q'),
    ranges = IRanges(start = c(1, 100, 30),
                     end = c(40, 150, 54)))
  
  
  s1000 <- smoothing(gr, 1000) #all segments are merged
  s1000_expected <- GRanges(
    seqnames = c('1p', '1q'),
    ranges = IRanges(start = c(1, 30),
                     end = c(150, 54)))
  
  return(c(sort(s0) == sort(s0_expected),
    sort(s4) == sort(s0_expected),
    sort(s5) == sort(s5_expected),
    sort(s6) == sort(s5_expected),
    sort(s1000) == sort(s1000_expected)))
}


#
# Test trimming
#
test_trim <- function(){
  gr <- GRanges(
    seqnames = c('1p', '1p', '1p', '1q', '1q'),
    ranges = IRanges(start = c(1, 10, 100, 30, 30),
                       end = c(9, 20, 150, 40, 50)))
  
  t0 <- trim(gr, 0) 
  #nothing should change
  
  t8 <- trim(gr, 8) 
  #nothing should change
  
  t9 <- trim(gr, 9) 
  t9_expected <- GRanges(
    seqnames = c('1p', '1p', '1q', '1q'),
    ranges = IRanges(start = c(10, 100, 30, 30),
                       end = c(20, 150, 40, 50)))
  
  t10 <- trim(gr, 10)
  #should be the same as t9
  
  
  t11 <- trim(gr, 11)
  t11_expected <- GRanges(
    seqnames = c('1p', '1q'),
    ranges = IRanges(start = c(100, 30),
                     end = c(150, 50)))
  
  t1000 <- trim(gr, 1000)
  #should eliminate all segments
  
  return(c(sort(t0) == sort(gr),
    sort(t8) == sort(gr),
    sort(t9) == sort(t9_expected),
    sort(t10) == sort(t9_expected),
    sort(t11) == sort(t11_expected),
    length(t1000) == 0))
}


test_smoothing()
test_trim()

