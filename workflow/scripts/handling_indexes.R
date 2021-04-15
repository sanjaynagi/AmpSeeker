library(data.table)
library(tidyverse)
library(glue)



### metdata


metadata = fread("config/samples.tsv")
metadata = metadata %>% rename("i7_index" = "i5_index", "i5_index" = "i7_index")


present = fread("resources/present_samples.txt") %>% deframe()

present = present %>% 
  str_remove("resources/reads/") %>% 
  str_remove("_1_fq.gz") %>% 
  str_remove("_2_fq.gz") %>% 
  unique()

present %>% list() %>% fwrite("present.samples.txt", sep="\t")

#fwrite(metadata, "config/samples.tsv", sep="\t", col.names = TRUE)

metadata %>% 
  mutate("barcode" = paste0(i5_index,"+", i7_index)) %>% 
  select(barcode, sampleID) %>% 
  fwrite("resources/reads/temp-names.txt", sep="\t", col.names = FALSE)


library(Biostrings)
metadata$i5_index = as.vector(reverseComplement(DNAStringSet(metadata$i5_index)))

#### write out barcodes, it should be revcomp(i5) + i7
metadata %>% 
  mutate("barcode" = paste0(i5_index,"+", i7_index)) %>% 
  select(barcode) %>% 
  fwrite("resources/reads/barcodes2.txt", sep="\t", col.names = FALSE)


metadata %>% 
  fwrite("config/samples_i5rev.tsv", sep="\t")





#### AmpSeq stats ####

ampseq = fread("resources/AGAMDAO_LSTM2/AmpSeqStats_290321.csv")

ampseq = file.info(list.files("resources/AGAMDAO_LSTM2/", full.names = TRUE)) %>% 
  rownames_to_column("path") %>% 
  select(path, size)

ampseq$name = substring(ampseq$path, 26,100)
ampseq$sample = substring(ampseq$path, 26)
ampseq$sample = str_extract(ampseq$name, '^[^_]+')
ampseq = ampseq[5:nrow(ampseq), c(4,3,2)]


ampseq = ampseq %>%
  mutate("readNo" = case_when(str_detect(name, "R1") == TRUE ~ "R1",
                              str_detect(name, "R2") == TRUE ~ "R2")) %>% 
  pivot_wider(., names_from = readNo,
                       values_from = c(name, size)) %>% arrange(sample) %>% 
  mutate("ID" = str_extract(sample, "[0-9]+"))

ampseq = ampseq %>% arrange(as.numeric(ID))
ampseq$plate = c(rep(1:16, each=96), 0)

head(ampseq)

tapply(ampseq$size_R1,ampseq$plate, sum)/1000

### demultiplex stats 

sum(unlist(ampseq$size_R1) > 100000)
