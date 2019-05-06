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
  fn4 <- 'testfiles/H18003069_gene_list_full_location.txt'

  segs1 <- getSegments(fn1, chrtable)
  segs2 <- getSegments(fn2, chrtable)
  segs3 <- getSegments(fn3, chrtable)
  segs4 <- getSegments(fn4, chrtable)
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

#
# Test the function sum_seg
#
test_sum <- function(){
  fn <- 'testfiles/testmix_gene_list_full_location.txt'
  chrtable <- getChrTable()
  
  #Get segments
  segs <- getSegments(fn, chrtable)
  
  gains <- segments_alt(segs, 'Gain')
  amps <- segments_alt(segs, 'Amp')
  losses <- segments_alt(segs, 'Loss')
  lohs <- segments_alt(segs, 'LOH')
  
  #Compute sum
  sum_gains <- sum_seg(gains)
  sum_amps <- sum_seg(amps)
  sum_losses <- sum_seg(losses)
  sum_lohs <- sum_seg(lohs)
  sum_all <- sum_seg(segs)
  
  #Define expected output
  arm_levels <- levels(seqnames(gains))
  
  expected_sum_gains <- rep(0, length(arm_levels))
  names(expected_sum_gains) <- arm_levels
  expected_sum_gains['7q'] <- 23201
  
  expected_sum_amps <- rep(0, length(arm_levels))
  names(expected_sum_amps) <- arm_levels
  expected_sum_amps['7p'] <- 2806224
  expected_sum_amps['7q'] <- 55897885
  
  expected_sum_losses <- rep(0, length(arm_levels))
  names(expected_sum_losses) <- arm_levels
  expected_sum_losses['1p'] <- 3816201
  expected_sum_losses['6p'] <- 50201
  expected_sum_losses['20p'] <- 38801
  expected_sum_losses['20q'] <- 314901
  expected_sum_losses['21q'] <- 330401
  
  expected_sum_lohs <- rep(0, length(arm_levels))
  names(expected_sum_lohs) <- arm_levels
  expected_sum_lohs['1p'] <- 3816201
  expected_sum_lohs['2q'] <- 2527001
  expected_sum_lohs['21q'] <- 330401
  expected_sum_lohs['7p'] <- 30001
  
  expected_sum_all <- rep(0, length(arm_levels))
  names(expected_sum_all) <- arm_levels
  expected_sum_all['1p'] <- 3816201
  expected_sum_all['2q'] <- 2527001
  expected_sum_all['6p'] <- 50201
  expected_sum_all['7p'] <- 2815424
  expected_sum_all['7q'] <- 55901084
  expected_sum_all['20p'] <- 38801
  expected_sum_all['20q'] <- 314901
  expected_sum_all['21q'] <- 330401
  
  
  return(c(
    amps=sum(expected_sum_amps == sum_amps) == length(expected_sum_amps),
    gains=sum(expected_sum_gains == sum_gains) == length(expected_sum_gains),
    lohs=sum(expected_sum_lohs == sum_lohs) == length(expected_sum_lohs),
    losses=sum(expected_sum_losses == sum_losses) == length(expected_sum_losses),
    all=sum(expected_sum_all == sum_all) == length(expected_sum_all)))

}


test_longest <- function(){
  fn <- 'testfiles/testmix_gene_list_full_location.txt'
  chrtable <- getChrTable()
  
  #Get segments
  segs <- getSegments(fn, chrtable)
  
  gains <- segments_alt(segs, 'Gain')
  amps <- segments_alt(segs, 'Amp')
  
  #Compute sum
  ampsgains_merged <- smoothing(c(amps, gains), 0)
  
  return(c(
  amps = sort(amps)[c(2,3)] == sort(longest(amps)),
  ampsgains = sort(amps)[c(2,3)] == sort(longest(c(amps, gains))),
  ampsgains_merged = sort(ampsgains_merged)[c(1,2)] == sort(longest(ampsgains_merged))
  ))
}

test_smoothing()
test_trim()
test_getSegments()
test_sum()
test_longest()
