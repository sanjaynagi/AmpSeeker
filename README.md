<div align="center">

[<img src="https://raw.githubusercontent.com/sanjaynagi/AmpSeeker/main/docs/ampseeker-docs/logo.png" width="400"/>](https://raw.githubusercontent.com/sanjaynagi/AmpSeeker/main/docs/ampseeker-docs/logo.png)   


[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub release](https://img.shields.io/github/release/sanjaynagi/AmpSeeker?include_prereleases=&sort=semver&color=blue)](https://github.com/sanjaynagi/AmpSeeker/releases/)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)

</div>

**Documentation**: https://sanjaynagi.github.io/AmpSeeker/ 

AmpSeeker is a snakemake workflow for Amplicon Sequencing data analysis for both Illumina and Nanopore amplicon sequencing. The pipeline is generic and can work on any data, but has extra modules tailored towards insecticide resistance monitoring. It implements:

- BCL to Fastq conversion
- Genome alignment
- Variant calling
- Quality control
- Coverage
- Visualisation of reads in IGV
- VCF to DataFrame/.xlsx 
- Allele frequency calculation
- Population structure
- Geographic sample maps
- Genetic diversity

- Kdr origin analysis (Ag-vampIR panel)
- Species assignment (Ag-vampIR panel)

The workflow uses a combination of papermill and jupyter book, so that users can visually explore the results in a local webpage for convenience.

## Usage

Please see the [documentation](https://sanjaynagi.github.io/AmpSeeker/) for more information on running the workflow.

## Citation 

**Targeted genomic surveillance of insecticide resistance in African malaria vectors**  
Nagi, *et al*., 2025. *bioRxiv*. doi: https://doi.org/10.1101/2025.02.14.637727

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [GitHub Actions](https://github.com/features/actions).

## Contributing to AmpSeeker

1. Fork the repository to your own GitHub user account
2. Clone your fork
3. Create a branch to implement the changes or features you would like `git checkout -b my_new_feature-24-03-23`
4. Implement the changes
5. Use `git add FILES`, `git commit -m COMMENT`, and `git push` to push your changes back to the branch
6. Open a Pull request to the main repository 
7. Once the pull request is merged, either delete your fork, or switch back to the main branch `git checkout main` and use `git pull upstream main` to incorporate the changes back in your local repo. Prior to `git pull upstream main`, you may need to set sanjaynagi/AmpSeeker as the upstream remote url, with `git remote set-url upstream git@github.com:sanjaynagi/AmpSeeker.git`. 
8. At this stage, your local repo should be up to date with the main Ampseeker branch and you are ready to start from #3 if you have more contributions!
