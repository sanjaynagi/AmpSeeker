# Bam statistics
# mapped to reduced amplicon only reference and whole genome 
library(data.table)
library(tidyverse)
library(glue)


## Read metadata 
metadata = fread("config/samples.tsv", sep="\t")
present = fread("resources/present.samples.txt", sep="\t", header = FALSE)
metadata = metadata[metadata$sampleID %in% present$V1]
#### targets
snps_targets = fread("resources/AgamDao_info.tsv", sep="\t")



flagstat = fread("resources/amplicon/alignments/bamStats/AgamDaoLSTM1_0058.flagstat", fill=TRUE)
flagstat




## Produce plots for each alignment method
for (ref in c("wholegenome","amplicon")){
  
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
  
  ### get avergaes per plate / well 
  pmap = tapply(pMapped, metadata$plate, median)
  mapped = tapply(mappedReads, metadata$plate, median)
  edge = tapply(mappedReads, metadata$well, median)
  
  ### make DF of summarised medians data per plate
  mappingStatsPlate = data.frame("MedianTotalReads" = tapply(totalReads, metadata$external_plate, median),
                            "MedianMappedReads" = tapply(mappedReads, metadata$external_plate, median),
                            "MedianPercentageMapped" = tapply(pMapped, metadata$external_plate, median),
                            "MedianDuplicatedReads" = tapply(dupReads, metadata$external_plate, median)) %>% 
    rownames_to_column("Plate")
  ### make DF of summarised medians data per well
  mappingStatsWell = data.frame("MedianTotalReads" = tapply(totalReads, metadata$well, median),
                            "MedianMappedReads" = tapply(mappedReads, metadata$well, median),
                            "MedianPercentageMapped" = tapply(pMapped, metadata$well, median),
                            "MedianDuplicatedReads" = tapply(dupReads, metadata$well, median)) %>% 
    rownames_to_column("Plate")
  # write to file 
  fwrite(mappingStatsPlate, glue("results/{ref}/qc/mappingStats_perPlate.tsv"), sep="\t", col.names = TRUE)
  fwrite(mappingStatsWell, glue("results/{ref}/qc/mappingStats_perWell.tsv"), sep="\t", col.names = TRUE)
  
  # plot histogram of totalReads
  pdf(glue("results/{ref}/qc/totalReads.pdf"))
  plt_tot = ggplot(data.frame(totalReads), aes(x=totalReads)) + 
    geom_density() + ggtitle(glue("{ref}_totalReads")) + 
    theme_light()
  print(plt_tot)
  null = dev.off()
  
  # plot duplicated reads
  pdf(glue("results/{ref}/qc/dupReads.pdf"))
  plt_dup = ggplot(data.frame(dupReads), aes(x=dupReads)) + 
    geom_density() +  ggtitle(glue("{ref}_dupReads")) + 
    theme_light()
  print(plt_dup)
  null = dev.off()
  
  # plot total mapped reads
  pdf(glue("results/{ref}/qc/mappedReads.pdf"))
  plt_mapped = ggplot(data.frame(mappedReads), aes(x=mappedReads)) + 
    geom_density() + ggtitle(glue("{ref}_mappedReads")) + 
    theme_light()
  print(plt_mapped)
  null = dev.off()
  
  # percentage mapped 
  pdf(glue("results/{ref}/qc/pMapped.pdf"))
  plt_pmapped = ggplot(data.frame(pMapped), aes(x=pMapped)) + 
    geom_density() + ggtitle(glue("{ref}_pMapped")) +  
    theme_light()
  print(plt_pmapped)
  null = dev.off()
  
  ### plot edge effects 
  edge = data.frame(edge) %>% rownames_to_column("Well") %>% 
      separate(Well,sep="(?<=[A-Za-z])(?=[0-9])", 
             into = c("Row", "column"))
  edge$column = as.numeric(edge$column)
  
  # edge effects 
  pdf(glue("results/{ref}/qc/edge_effects.pdf"))
  print(ggplot(edge, aes(x=column, y=reorder(Row, dplyr::desc(Row)), fill=edge)) +
    ggtitle(glue("Visualising edge effects - Median mapped Reads - {ref} reference")) + 
    geom_tile() + 
    scale_x_continuous(breaks = seq(1,12)) +
    theme_light() +
    labs(fill="Median mapped reads") + 
    ylab("Row") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 
  null = dev.off()

}

# mapped reads
metadata$totalReads = totalReads
metadata$pMapped = pMapped
metadata$MappedReads = mappedReads 
metadata$duplicated = dupReads



sessionInfo()

# 
# length(unique(metadata$external_plate))
# length(unique(metadata$plate))
