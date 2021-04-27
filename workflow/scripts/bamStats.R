# Bam statistics
# mapped to reduced amplicon only reference

library(data.table)
library(tidyverse)
library(glue)

## Read metadata 
metadata = fread("config/samples.tsv", sep="\t")
present = fread("resources/present.samples.txt", sep="\t", header = FALSE)
metadata = metadata[metadata$sampleID %in% present$V1]

#### targets
snps_targets = fread("resources/AgamDao_info.tsv", sep="\t")

## Produce plots for each alignment method
for (ref in c("amplicon", "wholegenome")){
  
  totalReads = c()
  dupReads = c()
  mappedReads = c()
  pMapped = c()
  
  for (sample in metadata$sampleID){
    flagstat = fread(glue("resources/{ref}/alignments/bamStats/{sample}.flagstat"), fill=TRUE)
    
    totalReads = c(totalReads, flagstat$V1[1])
    dupReads = c(dupReads, flagstat$V1[4])
    mappedReads = c(mappedReads, flagstat$V1[5])
    pMapped = c(pMapped, ((flagstat$V1[5]/flagstat$V1[1]) * 100))
    
  }
  
  pdf(glue("results/{ref}/qc/totalReads.pdf"))
  plt_tot = ggplot(data.frame(totalReads), aes(x=totalReads)) + 
    geom_density() + ggtitle(glue("{ref}_totalReads")) + 
    theme_light()
  null = dev.off()
  
  pdf(glue("results/{ref}/qc/dupReads.pdf"))
  plt_dup = ggplot(data.frame(dupReads), aes(x=dupReads)) + 
    geom_density() +  ggtitle(glue("{ref}_dupReads")) + 
    theme_light()
  null = dev.off()
  
  pdf(glue("results/{ref}/qc/mappedReads.pdf"))
  plt_mapped = ggplot(data.frame(mappedReads), aes(x=mappedReads)) + 
    geom_density() + ggtitle(glue("{ref}_mappedReads")) + 
    theme_light()
  null = dev.off()
  
  pdf(glue("results/{ref}/qc/pMapped.pdf"))
  plt_pmapped = ggplot(data.frame(pMapped), aes(x=pMapped)) + 
    geom_density() + ggtitle(glue("{ref}_pMapped")) +  
    theme_light()
  null = dev.off()
  
  ### get avergaes per plate / well 
  pmap = tapply(pMapped, metadata$plate, mean)
  mapped = tapply(mappedReads, metadata$plate, mean)
  edge = tapply(mappedReads, metadata$well, sum)
  
  ### plot edge effects 
  edge = data.frame(edge) %>% 
    rownames_to_column("Well") %>% 
    separate(Well,sep="(?<=[A-Za-z])(?=[0-9])", 
             into = c("Row", "column"))
  edge$column = as.numeric(edge$column)
  
  # edge effects 
  pdf(glue("results/{ref}/qc/edge_effects.pdf"))
  ggplot(edge, aes(x=column, y=reorder(Row, desc(Row)), fill=edge)) + 
    ggtitle(glue("Visualising edge effects - Average mapped Reads - {ref} reference")) + 
    geom_tile() + 
    scale_x_continuous(breaks = seq(1,12)) +
    theme_light() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  null = dev.off()

}

# mapped reads
metadata$pMapped = pMapped
metadata$MappedReads = mappedReads 
metadata$duplicated = dupReads





############## Where are the reads being aligned to? ##########

cover = list()

# read in mosdepth coverage data and store coverage per 300bp
for (sample in metadata$sampleID){
  cov = fread(glue("results/amplicon/coverage/{sample}.regions.bed.gz"))
  cover[[sample]] = cov$V4
}






### Overall means ###
# get mean of all positions across samples 
cov$TotalReads = apply(simplify2array(cover), 1, sum, na.rm=T)

# calc overall cov
cov = cov %>% select(V1, TotalReads) %>% 
  group_by(V1) %>% 
  summarise(AmpliconReads = sum(TotalReads)) %>% 
  mutate("TotalReads" = sum(AmpliconReads)) %>% 
  mutate("Perc" = (AmpliconReads/TotalReads)*100)







sessionInfo()

