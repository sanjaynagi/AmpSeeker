library(data.table)
library(tidyverse)
library(glue)


## whole genome coverage 
## Read metadata 
metadata = fread("config/samples.tsv", sep="\t")
present = fread("resources/present.samples.txt", sep="\t", header = FALSE)
metadata = metadata[metadata$sampleID %in% present$V1]
#### targets
snps_targets = fread("resources/AgamDao_info.tsv", sep="\t")
chroms = c("2L", "2R", "3R", "3L", "X")

# Read in coverage file and susbet to chroms of choice toi reduce size
cov = fread("results/wholegenome/coverage/windowed/AgamDaoLSTM1_0022.regions.bed.gz")
bool_ = cov$V1 %in% chroms
cov = cov[bool_,]

cover = list()
#read in mosdepth coverage data and store coverage per 300bp
start.time = Sys.time()
for (sample in metadata$sampleID){
  cover[[sample]] = fread(glue("results/wholegenome/coverage/windowed/{sample}.regions.bed.gz"), select = 4)[bool_,] %>% deframe()
}
end.time = Sys.time()
time.taken = end.time - start.time
time.taken


df = apply(simplify2array(cover), 1, median, na.rm=T)
head(df)
max(df)


str(cover)
# get mean of all positions across samples 
mediancov = apply(simplify2array(cover), 1, median, na.rm=T)
max(mediancov)
max(cov$V4)


pdf(glue("results/AgamDao_wholegenome_coverage.pdf"))
#### coverage per chrom 
for (chr in chroms){
  
  bool_ = cov$V1 == chr
  pos = cov$V2[bool_]
  
  coverage = mediancov[bool_]
  
  df = data.frame("pos" = pos, "cov" = coverage)
  targets = snps_targets %>% filter(chrom == chr)
  
  plt = ggplot(df, aes(x=pos, y=cov)) + 
    geom_line() +
    ggtitle(glue("AgamDao Coverage - {chr}")) + 
    ylim(0, 100) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme_light() + 
    geom_vline(xintercept = targets$snploc, colour="grey", alpha=0.5)
  
  print(plt)
}

dev.off()


# Find target regions using foverlaps from data table 
cov$V4 = mediancov
colnames(cov) = c("chrom", "start", "end", "medianCov")

# make bed like data frame 
snps = snps_targets %>% 
  select(chrom, snploc) %>%
  mutate("end" = snploc+1) %>% rename(start=snploc)

head(snps)
head(cov)
#set Key for  foverlaps
setkey(cov)
#run foverlaps - merge based on SNP + amplicon ranges
onTargetdf = foverlaps(snps, 
                         cov, 
                         by.x = c("chrom", "start", "end"), 
                         by.y = c("chrom", "start", "end")) %>% select(chrom, start,end, medianCov) %>% distinct()

colnames(onTargetdf)


onTargetTotal = sum(onTargetdf$medianCov)
TotalCoverage = sum(mediancov)

# get % of on target reads 
(onTargetTotal/TotalCoverage)*100


