# Bam statistics
# mapped to reduced amplicon only reference

library(data.table)
library(tidyverse)
library(glue)


metadata = fread("config/samples.tsv", sep="\t")
present = fread("resources/present.samples.txt", sep="\t", header = FALSE)
metadata = metadata[metadata$sampleID %in% present$V1]


totalReads = c()
dupReads = c()
mappedReads = c()
pMapped = c()

for (sample in metadata$sampleID){
  flagstat = fread(glue("resources/wholegenome/alignments/bamStats/{sample}.flagstat"), fill=TRUE)
  
  totalReads = c(totalReads, flagstat$V1[1])
  dupReads = c(dupReads, flagstat$V1[4])
  mappedReads = c(mappedReads, flagstat$V1[5])
  pMapped = c(pMapped, ((flagstat$V1[5]/flagstat$V1[1]) * 100))
  
}


ggplot(data.frame(totalReads), aes(x=totalReads)) + 
  geom_density() + 
  theme_light()


ggplot(data.frame(dupReads), aes(x=dupReads)) + 
  geom_density() + 
  theme_light()

ggplot(data.frame(mappedReads), aes(x=mappedReads)) + 
  geom_density() + 
  theme_light()

ggplot(data.frame(pMapped), aes(x=pMapped)) + 
  geom_density() + 
  theme_light()


metadata$pMapped = pMapped


tapply(metadata$pMapped, metadata$plate, mean)

tapply(mappedReads, metadata$plate, mean)
tapply(mappedReads, metadata$well, mean)

exclude = metadata[metadata$pMapped < 50,]


############## Where is the coverage being aligned to? ##########

cov = fread(glue("results/wholegenome/coverage/{sample}.regions.bed.gz"))


cover = list()

#read in mosdepth coverage data and store coverage per 300bp
for (sample in metadata$sampleID){
  cov = fread(glue("results/coverage/{sample}.regions.bed.gz"))
  
  cover[[sample]] = cov$V4
}

#get mean of all positions across samples 
mediancov = apply(simplify2array(cover), 1, median, na.rm=T)


df$median = apply(df, 1, median, na.rm = T).

#### coverage per chrom 
chroms = c('2L', '2R', '3L', '3R', 'X')

for (chrom in chroms){
  
  bool_ = cov$V1 == chrom
  pos = cov$V2[bool_]
  
  coverage = meancov[bool_]
  
  df = data.frame("pos" = pos, "cov" = coverage)
  
  plt = ggplot(df, aes(x=pos, y=cov)) + geom_line() + ggtitle(chrom) + ylim(0, 200)
  
  print(plt)
  print(paste(chrom, sum(df$cov)))
}

targets[targets$chrom == '2L']

targets[]


sessionInfo()

