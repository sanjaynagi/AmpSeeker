### Prep primer seqs for bwa aln 

library(tidyverse)
library(data.table)
library(glue)

## Read metadata 
metadata = fread("config/samples.tsv", sep="\t")
present = fread("resources/present.samples.txt", sep="\t", header = FALSE)
metadata = metadata[metadata$sampleID %in% present$V1]
## targets
snps_targets = fread("resources/AgamDao_info.tsv", sep="\t")

## primer sequences
primers = fread("resources/AgamDaoPrimerSeqs.tsv", sep="\t")

## substring 
primers$fwd = substring(primers$fwd, 34, 100)
primers$rev = substring(primers$rev, 34, 100)

## m, d, y, u
primers = pivot_longer(primers, cols=c(fwd, rev), values_to = "seq") %>% 
  mutate("name" = paste0(Name,"_", name)) %>% 
  select(-Name)

## write Fasta
writeFasta = function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn = file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

## write fasta
writeFasta(primers, "resources/AgamDaoPrimerSeqs_nolinker.fa")

#######
# bwa aln -n 2 ~/ag1000g/data/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa $fasta |
# bwa samse ~/ag1000g/data/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa - $fasta
####

library(Rsamtools)

sam = fread("resources/AgamDaoPrimerSeqs_nolinker.fa.sam", skip = 9, fill=TRUE) %>% select(V1, V2, V3, V4, V10, V20)
head(sam)

sam = sam %>% separate(., V20, sep=",", into=c("chrom2", "pos2"))
head(sam)

sam$chrom2 = str_extract(sam$chrom2, "2L|2R|3R|3:L|X")
sam$pos2 = str_remove(sam$pos2, "\\+|\\-")



