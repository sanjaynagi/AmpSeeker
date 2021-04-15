library(data.table)
library(tidyverse)
library(stringdist)
library(itertools)

samples = fread("resources/210326_M05658_0003_000000000-J54FT/SampleSheet.csv")
samples = samples[20:nrow(samples),]
colnames(samples) = unlist(samples[1,])
samples = samples[2:nrow(samples),] 

i5 = samples$index
i7 = samples$index2

i5_distmatrix = stringdistmatrix(unique(i5), unique(i5), method="hamming")
i7_distmatrix = stringdistmatrix(unique(i7), unique(i7), method="hamming")

indices = samples %>% select(index, index2, Sample_ID)
fwrite(indices, "resources/indices.txt", sep="\t", row.names = FALSE)
fwrite(list(indices$index), "i5.txt", row.names = FALSE)
fwrite(list(indices$index2), "i7.txt", row.names = FALSE)

#write fasta 
for (i in 1:length(unique(indices$index))){
  cat(paste0("> i5Tag-", i ,"\n", unique(indices$index)[i]), sep = "\n", file ="resources/i5.fasta", append = TRUE)
}

#write fasta 

for (i in 1:length(unique(indices$index2))){
  cat(paste0("> i7Tag-", i ,"\n", unique(indices$index2)[i]), sep = "\n", file ="resources/i7.fasta", append = TRUE)
}


#### Full adaptor sequences for bbduk to remove

i5 = read.xlsx("~/projects/AmpSeq/reagents/IDT_i5_adaptor_order_SN.xlsx")
i7 = read.xlsx("~/projects/AmpSeq/reagents/IDT_i7_adaptor_order_SN.xlsx")

i5$Sequence = str_remove(i5$Sequence, "\\*")
i7$Sequence = str_remove(i7$Sequence, "\\*")

#write fasta 
for (i in 1:length(unique(i5$Sequence))){
  cat(paste0("> i5Tag-", i ,"\n", unique(i5$Sequence)[i]), sep = "\n", file ="resources/i5.full.fasta", append = TRUE)
}

#write fasta 

for (i in 1:length(unique(i7$Sequence))){
  cat(paste0("> i7Tag-", i ,"\n", unique(i7$Sequence)[i]), sep = "\n", file ="resources/i7.full.fasta", append = TRUE)
}

#### Panel info 

library(data.table)
library(tidyverse)
library(openxlsx)

agamdao = read.xlsx("../AmpSeq/IR_AGAMDAO_panel_info.xlsx")
agamdao$Location = str_remove(agamdao$Location, "Ag_")


agamdao = agamdao %>% mutate("loc" = paste0(Location, "-", str_remove(agamdao$Location, "2L:|2R:|3R:|3L:|X:")))

fwrite(agamdao, "resources/agamdao.tsv", sep="\t", col.names = TRUE)

agamdao %>% select(loc) %>% fwrite(., "resources/AgamDaoLocs.txt", col.names = FALSE, sep="\t")



test = str_split(agamdao$loc, ":")
agamdao$chrom = sapply(test, "[[", 1)
agamdao$start = sapply(test, "[[", 2)
agamdao$start = as.numeric(sapply(two, "[[", 1)) -1
agamdao$end = as.numeric(agamdao$start) + 1

two = str_split(agamdao$start,"-")

agamdao %>% 
  select(chrom, start, end) %>% 
  arrange(chrom, start) %>% 
  fwrite(., "resources/AgamDaoLoci.bed", sep="\t", col.names = FALSE)



bed = fread("resources/AgamDao.bed")
bed$Chromosome = str_remove(bed$Chromosome, "Ag_")
bed %>% select(Chromosome, Start, Stop) %>% fwrite(., "resources/AgamDao.bed", col.names = FALSE, sep="\t")
bed = bed %>% mutate("hello" = paste0(Chromosome,":" ,Start-50,"-", Stop+50)) #%>% 
  #select(hello) %>%  
  #fwrite(., "resources/AgamDao.bed", col.names = FALSE, sep="\t")

bed$Start - bed$Stop            
            



####  
library(Biostrings)

#fwrite(metadata, "config/samples.tsv", sep="\t", col.names = TRUE)

metadata %>% 
  mutate("barcode" = paste0(i7_index,"+", i5_index)) %>% 
  select(barcode, sampleID) %>% 
  fwrite("resources/reads/temp-names.txt", sep="\t", col.names = FALSE)

metadata %>% 
  mutate("barcode" = paste0(i5_index,"+", i7_index)) %>% 
  select(barcode) %>% 
  fwrite("resources/reads/barcodes2.txt", sep="\t", col.names = FALSE)

library(Biostrings)



      