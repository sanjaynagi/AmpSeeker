# Bam statistics
# mapped to reduced amplicon only reference and whole genome 
library(data.table)
library(tidyverse)
library(glue)

############## Where are the reads being aligned to? ##########
metadata = fread("AmpSeq2023/config/metadata.tsv")





# Need to edit column name 
cov = fread("AmpSeq2023/results_amplicon/coverage/sample97.regions.bed.gz")[,1:3]
# read in mosdepth coverage data and store coverage per 300bp
for (sample in metadata$sampleID){
  covdf = fread(glue("AmpSeq2023/results_amplicon/coverage/{sample}.regions.bed.gz")) %>% dplyr::rename(!!sample := V4) %>% select(-c(V1,V2,V3))
  cov = cbind(cov, covdf)
}

## summarise each 
coverage_df = cov %>% mutate("targetID" = paste0(V1,"_", V2)) %>% 
  select(-c(V1,V2,V3)) %>% 
  column_to_rownames("targetID")


### Rebalancing ####
## cal total per sample 
total_per_sample = apply(coverage_df, FUN=sum, MARGIN = 2)

# Target read fractions 
targetReadFractionsdf = as.data.frame(t((t(coverage_df)) / total_per_sample))
targetReadFractionsdf = targetReadFractionsdf[,colSums(is.na(targetReadFractionsdf)) == 0]

# Median read fraction per target 
median_per_target = apply(targetReadFractionsdf, FUN = median, MARGIN = 1)

# Sum of Medians 
medianTotal = sum(median_per_target)

# Scale the read fractions to 1 
scaled_median_per_target = median_per_target * (1/medianTotal)

# to the power baby 0.561 oh yeah 
primer_volumes = scaled_median_per_target^-0.561
primer_volumes

# scale by the lowest weighting so that minimum weighting is 1
# for pipetting reasons and max is 10ul 
primer_volumes = primer_volumes/min(primer_volumes)
primer_volumes[primer_volumes > 10] = 10

# total 
total = sum(primer_volumes)

#get IQ mean 
iq = quantile(primer_volumes)
iq
# central 
central_primer_conc = iq[3] /total * 250000

#calc dilution factor
dilution_fac = central_primer_conc /40

# what are the primer vols
primer_volumes






sessionInfo()
