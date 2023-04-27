# AmpSeeker

Welcome to [AmpSeeker](https://github.com/sanjaynagi/AmpSeeker/), an Amplicon Sequencing data analysis pipeline written in Snakemake.

AmpSeeker is a compuational pipeline designed to analyse high-throughput Illumina MiSeq whole-genome sequencing data of custom amplicon panels. These panels target multiple loci of interest in a single assay. AmpSeeker is capable of analysing amplicon sequence data for potentially any organism, provided an appropriate reference sequence, SNP targets, and target genomic regions. Currently, the workflow allows users to perform genome alignment, sequence coverage calculation, and variant calling. In addition, the workflow produces an interactive IGV-Notebook, which facilitates read alignment investigation, and an HTML report that includes a diagrammatic job execution graph along with computation times of all software tools. Workflow tools are defined in Conda, ensuring reproducible computation that is platform-agnostic.


```{tableofcontents}
```