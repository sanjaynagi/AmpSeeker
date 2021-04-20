# amplicon only reference 

library(data.table)
library(tidyverse)
library(glue)
library(openxlsx)

library(valr)

#### chroms 
chroms = c('2L', '2R', '3L', '3R', 'X')

### read snp targets, and 1000 bases +-
targets = fread("resources/AgamDaoLoci.bed")
targets$V2 = targets$V2 - 1000
targets$V3 = targets$V3 + 1000

colnames(targets) = c("chrom", "start", "end")

# merge enlarged amplicons so there is no overlap (valr package)
targets_noverlap = bed_merge(targets)

# mutate to join chrom start and end 
targets_noverlap %>% 
  mutate("col" = paste0(chrom, ":", start, "-", end)) %>% 
  select(col) %>% 
  fwrite(., "resources/AgamDao_nonoverlaps.regions", 
         sep="\t", col.names = FALSE)

# mutate to join chrom start and end 
#targets_noverlap %>% 
#  mutate("chr" = paste0(chrom, ":", start, "-", end))

targets_noverlap = targets_noverlap %>% 
  mutate("amplicon_id" = glue("Agam_{row_number()-1}")) %>% as.data.table()


# read in SNP data 
snps = read.xlsx("../AmpSeq/IR_AGAMDAO_panel_info.xlsx")
snps = snps %>% select(TargetName, `Mutation.(blank.if.synonymous/unknown)`, Gene, Location)
colnames(snps) = c("TargetName", "Mutation", "Gene", "Location")

# correct SNP target data to format for foverlaps
snps = snps %>% mutate("chrom" = str_extract(Location, "2R|2L|3R|3L|X"),
                "start" = as.numeric(str_remove(Location, "Ag_2L:|Ag_2R:|Ag_3R:|Ag_3L:|Ag_X:"))) %>% 
  mutate("end" = start +1) %>% as.data.table()

### move bed columns to beginning
snps = snps %>% select(chrom, start,end, everything())

#set Key for  foverlaps
setkey(targets_noverlap)
#run foverlaps - merge based on SNP + amplicon ranges
snps_targets = foverlaps(snps, 
          targets_noverlap, 
          by.x = c("chrom", "start", "end"), 
          by.y = c("chrom", "start", "end"))

# find pos of target SNP in amplicon
snps_targets = snps_targets %>% 
  select(-i.end) %>% 
  dplyr::rename("snploc" = "i.start") %>% 
  mutate("pos" = snploc-start)

#Write BED file of targets 
snps_targets %>% 
  select(amplicon_id, pos) %>% 
  mutate(end = pos+1) %>% 
  arrange(amplicon_id) %>% 
  fwrite("resources/AgamDao_amplicon_snptargets.bed", sep="\t", col.names = FALSE)



# targets noverlap
targets_noverlap %>% 
  mutate("loc" = paste0(chrom, ":",  start, "-", end)) %>% 
  select(loc, amplicon_id) %>% 
  fwrite("resources/amplicon_names.txt", sep="\t", col.names = FALSE)
                            
sessionInfo()

                            