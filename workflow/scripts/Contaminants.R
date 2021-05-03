# contaminant 

library(tidyverse)
library(data.table)
library(glue)
library(Biostrings)

#####################

metadata = fread("config/samples.tsv")
#metadata = metadata %>% rename("i7_index" = "i5_index", "i5_index" = "i7_index")

#seq = "TCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
#pcrdimer = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
#pcrdimer2 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
#DNAStringSet(pcrdimer2)
#reverseComplement(DNAStringSet(pcrdimer))


### Read in contamination stats
bbduk = list()

matched = c()
pcr_dimer = c()
TruSeq = c()
pcr_primers = c()

for (sample in metadata$sampleID){
  bbduk = fread(glue("resources/reads/trimmed/stats/{sample}.txt"), sep="\t")
  
  pcr_dimer = c(pcr_dimer, as.character(bbduk[V1 == "pcr_dimer"][,3]))
  TruSeq = c(TruSeq, as.character(bbduk[V1 == "TruSeq_Universal_Adapter"][,3]))
  matched = c(matched, as.character(bbduk[V1 == "#Matched"][,3]))
  pcr_primers = c(pcr_primers, as.character(bbduk[V1 == "PCR_Primers"][,3]))
  
} 


pcr_dimer = as.numeric(str_remove(pcr_dimer, "%"))
pcr_primers = as.numeric(str_remove(pcr_primers, "%"))
matched = as.numeric(str_remove(matched, "%"))
TruSeq = as.numeric(str_remove(TruSeq, "%"))


tapply(pcr_dimer, metadata$external_plate, median)
tapply(matched, metadata$external_plate, median)
tapply(pcr_primers, metadata$external_plate, median)
tapply(TruSeq, metadata$external_plate, median)



### make DF of summarised medians data per plate
ContamStatsPlate = data.frame("Medianpcr_dimer" = tapply(pcr_dimer, metadata$external_plate, median),
                               "MediantotalMatches" = tapply(matched, metadata$external_plate, median),
                               "MedianPCRPrimers" = tapply(pcr_primers, metadata$external_plate, median),
                               "MedianTruSeq" = tapply(TruSeq, metadata$external_plate, median)) %>% 
  rownames_to_column("Plate")

# write to file 
fwrite(ContamStatsPlate, glue("results/ContamStats_perPlate.tsv"), sep="\t", col.names = TRUE)





