# amplicon only reference 

library(data.table)
library(tidyverse)
library(glue)

library(valr)
####
chroms = c('2L', '2R', '3L', '3R', 'X')

targets = fread("resources/AgamDaoLoci.bed")
targets$V2 = targets$V2 - 1000
targets$V3 = targets$V3 + 1000

colnames(targets) = c("chrom", "start", "end")

targets_noverlap = bed_merge(targets)


targets_noverlap %>% 
  mutate("col" = paste0(chrom, ":", start, "-", end)) %>% 
  select(col) %>% 
  fwrite(., "resources/AgamDao_nonoverlaps.regions", 
         sep="\t", col.names = FALSE)


snps = fread("resources/AgamDaoLoci.bed")
colnames(snps) = c("chrom", "start", "end")

targets_noverlap %>% 
  mutate("chr" = paste0(chrom, ":", start, "-", end))

targets$start





for (chr in chroms){
    
  targ=targets_noverlap[targets_noverlap$chrom == chr,]
  
  for (i in 1:nrow(snps)){
    #print(snps$start[i])
    #print(targ[i,])
    bool_ = snps$start[i]-1000 == targ$start[i]
    
    a = snps$start[bool_] - targ$start[i]
    
    
    print(a)
  }
}
snps$start[snps$start-1000 == targ$start[3]] - targ$start[3]

