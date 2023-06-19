# AmpSeeker

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.11.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![AmpSeq](https://github.com/sanjaynagi/AmpSeeker/workflows/AmpSeq/badge.svg)](https://github.com/sanjaynagi/AmpSeeker/actions?query=workflow:"AmpSeq")
[![GitHub release](https://img.shields.io/github/release/sanjaynagi/AmpSeeker?include_prereleases=&sort=semver&color=blue)](https://github.com/sanjaynagi/AmpSeeker/releases/)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)

**Documentation**: https://sanjaynagi.github.io/AmpSeeker/ 

AmpSeeker is a snakemake workflow for Amplicon Sequencing data analysis. The pipeline is a work in progress, however, it currently implements:

- BCL to FastQ conversion
- Genome alignment
- Coverage calculation
- Visualisation of reads in IGV-notebook
- Variant calling
- Principal component analysis
- Geographic sample maps

The workflow uses a combination of papermill and jupyter book, so that users can visually explore the results in a local webpage for convenience.

## Authors

* Sanjay C Nagi (@sanjaynagi)
* Trevor Mugoya (@ChabbyTMD)
* Edward Lukyamuzi (@eddUG) 

## Usage

Please see the [wiki](https://github.com/sanjaynagi/AmpSeeker/wiki) for more information on running the workflow. If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).

## Contributing to AmpSeeker

1. Fork the repository to your own github user account
2. Clone your fork
3. Create a branch to implement the changes or features you would like `git checkout -b my_new_feature-24-03-23`
4. Implement the changes
5. Use `git add FILES`, `git commit -m COMMENT`, and `git push` to push your changes back to the branch
6. Open a Pull request to the main repository 
7. Once the pull request is merged, either delete your fork, or switch back to the main branch `git checkout main` and use `git pull upstream main` to incorporate the changes back in your local repo. Prior to `git pull upstream main`, you may need to set sanjaynagi/AmpSeeker as the upstream remote url, with `git remote set-url upstream git@github.com:sanjaynagi/AmpSeeker.git`. 
8. At this stage, your local repo should be up to date with the main Ampseeker branch and you are ready to start from #3 if you have more contributions!
