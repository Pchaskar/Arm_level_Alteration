# Load all oncoscan data
source("functions.R")

chrtab <- getChrTable()
datadir <- 'Q:/labo biol mol/Oncoscan/genelists-wLocation'
files <- list.files(path = datadir, pattern = '*_gene_list_full_location.txt')

segs.perFile <- sapply(files, function(fn){
  getSegments(paste0(datadir, '/', fn), chrtab)
})

saveRDS(segs.perFile, file = 'segs_perFile.rds')
