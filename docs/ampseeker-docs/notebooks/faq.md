# Frequently Asked Questions (FAQ)

## General

### What does AmpSeeker do?
AmpSeeker is a Snakemake workflow for amplicon sequencing analysis. It supports Illumina and Nanopore data, performs QC, alignment, variant calling, annotation, and produces a Jupyter Book report.

### Who is AmpSeeker intended for?
AmpSeeker is designed for users who want a reproducible, automated pipeline for amplicon sequencing datasets, including cohort-level analyses and visual reporting.

### Do I need to install lots of tools manually?
No. The intended usage is with `--use-conda`, which creates per-rule environments automatically.

### Where should I run the workflow from?
Run Snakemake from the project root (the directory containing `workflow/`, `config/`, and `resources/`).

## Inputs and Configuration

### Which config file is used by default?
By default, AmpSeeker reads `config/config.yaml`. You can override this with:

```bash
snakemake --configfile path/to/your-config.yaml
```

### What is the minimum metadata requirement?
At minimum, metadata must contain `sample_id`.

For Illumina/Nanopore direct FASTQ mode (`from-bcl: False`):
- Illumina expects either:
  - `fq1` and `fq2` columns, or
  - automatic FASTQ naming in `resources/reads/` if `fq1/fq2` are absent.
- Nanopore expects `fq1`.

### What are `cohort-columns` used for?
They must exist in the metadata and are used to group/colour samples in notebooks and downstream analyses.

### What should the BED file look like?
The `targets` BED file should define target loci for your panel and match the reference coordinates used for mapping/calling.

## Illumina Modes

### Can AmpSeeker run from BCL folders?
Yes. Set:

```yaml
platform: illumina
from-bcl: True
```

Then provide `illumina-dir`.

### Can `illumina-dir` be a single path or multiple paths?
Yes. It supports both:

Single run:

```yaml
illumina-dir: resources/250110_M05658_0028_000000000-LTBV4
```

Multiple runs:

```yaml
illumina-dir:
  - resources/run_001
  - resources/run_002
```

### How does multi-run Illumina processing work?
AmpSeeker converts each run independently, then merges demultiplexed FASTQs per sample into standard outputs:
- `resources/reads/{sample}_1.fastq.gz`
- `resources/reads/{sample}_2.fastq.gz`

Run stats are also aggregated for downstream reporting.

### What if a sample appears in more than one run?
That is supported when metadata/sample sheet values are consistent for that sample. FASTQs are merged across runs.

### What if duplicate sample IDs conflict across runs?
AmpSeeker raises an error if duplicated `sample_id` entries across SampleSheets disagree on metadata values.

### Does this change direct FASTQ workflows?
No. `from-bcl: False` behavior is unchanged.

## Nanopore

### Is `from-bcl` valid for Nanopore?
No. Set `from-bcl: False` for Nanopore.

### How are Nanopore FASTQs provided?
Using metadata `fq1` paths.

## Running the Workflow

### What command should I use?
Typical run:

```bash
snakemake --cores 4 --use-conda
```

Dry run first:

```bash
snakemake -n --cores 4 --use-conda
```

### How do I run in `.test/`?
Example:

```bash
snakemake --cores 1 --directory .test --use-conda --configfile .test/config/config_fastqauto.yaml -n
```

### What if Snakemake says files are missing?
Check:
- config paths are correct relative to your run directory.
- `metadata`, `targets`, and reference files exist.
- FASTQ naming/paths match your chosen mode.

## Outputs

### Where are results written?
Primary outputs are under `results/`.

### Where is the report book?
The results book is built under `results/ampseeker-results/_build/html/` when `build-jupyter-book: True`.

### Why are some notebooks missing from outputs?
Notebook generation is conditional on config options (platform, QC toggles, analysis toggles, and panel-specific modules).

## Troubleshooting

### Error: unsupported platform
Set `platform` to either `illumina` or `nanopore`.

### Error: metadata file does not exist
Verify the `metadata` path is correct for your current working directory and config file.

### Error: `illumina-dir` must be provided
This occurs when `platform: illumina` and `from-bcl: True` but `illumina-dir` is empty.

### Error: cohort columns missing
Make sure every value listed in `cohort-columns` exists as a column in metadata (or SampleSheet-derived metadata in BCL mode).

### Why did index read QC fail?
Check whether index reads were produced in BCL conversion and whether combined index FASTQs exist in `resources/reads/`.

### Workflow reruns unexpectedly after config changes
Changes in config, inputs, or software environments can trigger reruns in Snakemake. Use dry runs to inspect what will execute.

## Reproducibility and Best Practices

### Should I edit files under `results/` manually?
No. Treat `results/` as generated output.

### Should I version-control my config and metadata?
Yes. Track your run config and metadata for reproducibility and provenance.

### How should I name datasets?
Use stable, descriptive `dataset` names (e.g., include project and run context) to keep outputs interpretable.

### What is the safest way to start on a new dataset?
1. Validate metadata format.
2. Dry-run with your final config.
3. Run with limited cores first.
4. Inspect QC notebooks before downstream interpretation.
