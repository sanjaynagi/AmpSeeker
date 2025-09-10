# Contributing to AmpSeeker

Welcome to the AmpSeeker contributing guide! This document provides detailed information about how the workflow works under the hood and how to contribute to or extend the pipeline for new panels and use cases.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Workflow Components](#workflow-components)
3. [Adding Support for New Panels](#adding-support-for-new-panels)
4. [Configuration System](#configuration-system)
5. [Development Guidelines](#development-guidelines)
6. [Testing](#testing)
7. [Troubleshooting](#troubleshooting)

## Architecture Overview

AmpSeeker is built using Snakemake, a workflow management system that enables scalable and reproducible data analysis. The pipeline follows a modular architecture with clear separation between:

- **Core functionality**: Generic amplicon sequencing analysis
- **Platform-specific processing**: Illumina vs. Nanopore handling
- **Panel-specific analyses**: Custom analyses for specific panels (e.g., Ag-vampIR)

### Key Design Principles

1. **Reproducibility**: All analyses are parameterized and version-controlled
2. **Modularity**: Components can be used independently and new modules easily added
3. **Interactive exploration**: Jupyter notebooks provide interactive result exploration
4. **Standardized outputs**: Consistent directory structure and file formats

### Directory Structure

```
├── config/                 # Configuration files and metadata
│   ├── config.yaml        # Main configuration file
│   ├── [panel].bed        # Target regions for specific panel
│   └── [panel]-metadata.tsv # Sample metadata
├── resources/             # Reference files and panel-specific resources
│   ├── reference/         # Genome references, GFF files
│   └── [panel]/          # Panel-specific resources
├── workflow/              # Core pipeline components
│   ├── ampseekertools.py # Utility functions
│   ├── envs/             # Conda environment specifications
│   ├── notebooks/        # Analysis notebooks
│   ├── rules/            # Snakemake rule files
│   └── Snakefile         # Main workflow entry point
└── results/              # Generated outputs (created at runtime)
```

## Workflow Components

### 1. Entry Point (`workflow/Snakefile`)

The main Snakefile orchestrates the entire workflow:

- Validates configuration parameters
- Loads metadata and determines sample processing strategy
- Conditionally includes platform-specific and panel-specific rules
- Defines the final output requirements through the `all` rule (see `AmpSeekerOutputs()` function within `rules/common.smk`)

**Key variables set in Snakefile:**
- `dataset`: Dataset name from config
- `panel`: Panel type from config
- `samples`: List of sample IDs from metadata
- `fastq_auto`: Boolean indicating if FASTQ paths are auto-generated
- `large_sample_size`: Boolean for >1000 samples (requires split merging)
- `plate_info`: Boolean indicating if plate layout information is available

### 2. Rule Files (`workflow/rules/`)

#### Core Rules
- **`common.smk`**: Shared functions, color mappings, output definitions
- **`metadata.smk`**: Metadata loading and validation
- **`qc.smk`**: Quality control with FastP, MultiQC
- **`qc-notebooks.smk`**: QC report generation notebooks
- **`from-bcl.smk`**: Illumina BCL to FASTQ conversion (if applicable)
- **`map-call-illumina.smk`** / **`map-call-nanopore.smk`**: Platform-specific alignment and variant calling
- **`snpeff.smk`**: Variant annotation
- **`analysis.smk`**: Population genetic analyses
- **`jupyter-book.smk`**: Report generation

#### Platform-Specific Processing

**Illumina Pipeline:**
```
Raw reads → FastP QC → BWA-MEM alignment → BCFtools variant calling → SnpEff annotation
```

**Nanopore Pipeline:**
```
Raw reads → FastP long QC → Minimap2 alignment → Clair3 variant calling → SnpEff annotation
```

#### Panel-Specific Rules
- **`ag-vampir.smk`**: Species identification and Kdr origin analysis for mosquito panels

### 3. Analysis Notebooks (`workflow/notebooks/`)

Jupyter notebooks implement all analyses using papermill for parameterization.
Papermill is a tool that allows you to execute Jupyter notebooks whilst passing in parameters, like a python script.  
This allows us to debug and develop the notebooks interactively, but then run them in an automated fashion as part of the Snakemake workflow.  

#### Quality Control Notebooks
- `read-quality.ipynb`: Sequencing quality metrics
- `coverage.ipynb`: Coverage analysis across targets
- `sample-quality-control.ipynb`: Sample filtering and QC metrics
- `run-statistics.ipynb`: Illumina run statistics (BCL processing)
- `reads-per-well.ipynb`: Per-well read distribution (plate-based data)

#### Population Genetics Notebooks
- `population-structure.ipynb`: PCA, neighbor-joining trees
- `genetic-diversity.ipynb`: Nucleotide diversity, Tajima's D
- `allele-frequencies.ipynb`: Allele frequency calculations
- `sample-map.ipynb`: Geographic visualization

#### Utility Notebooks
- `snp-dataframe.ipynb`: Convert VCF to DataFrame/Excel
- `IGV-explore.ipynb`: IGV session file generation (not part of workflow)

### 4. Utility Functions (`workflow/ampseekertools.py`)

Core functions used across the pipeline:

- **`load_vcf()`**: VCF loading with quality filtering
- **`create_color_mapping()`**: Consistent color schemes for visualizations
- **`get_fastqs()`**: FASTQ file path resolution
- **`AmpSeekerOutputs()`**: Dynamic output requirement definition

--- 

## Adding Support for New Panels

Adding a new panel to AmpSeeker involves several steps. Use the Ag-vampIR panel as a reference implementation.

### Step 1: Panel-Specific Resources

If the panel requires specific resource data, create a new directory under `resources/[panel-name]/`:

```
resources/
└── your-panel/
    ├── marker_SNPs.csv      # Panel-specific variant definitions
    ├── reference_data.tsv   # Any reference datasets
    └── ...                  # Other panel resources
```

### Step 2: Target Regions File

Create a BED file defining your SNP target region:

```
# config/your-panel.bed
chr1    1000    1001    amplicon_1
chr2    5000    5001    amplicon_2
...
```

### Step 3: Panel-Specific Analysis Notebooks

Create panel-specific notebooks in `workflow/notebooks/[panel-name]/`.
Copy how the Ag-vampIR notebooks are structured. 


```python
# Example notebook structure
import pandas as pd
import numpy as np
import allel

# Parameters (injected by papermill)
dataset = "your-dataset"
metadata_path = "path/to/metadata"
vcf_path = "path/to/vcf"
panel_resources = "resources/your-panel/"

# Load data
metadata = pd.read_csv(metadata_path, sep='\t')
vcf = allel.read_vcf(vcf_path)

# Panel-specific analysis logic
def your_panel_analysis():
    # Implement your specific analysis

# Generate outputs and visualizations
```

The second cell should contain the parameters that will be injected by papermill - it must be given the `parameters` tag. When notebooks  are executed, these parameters will be overwritten with those from the workflow.

To prevent code from being shown in the jupyter book, you will need to add the `remove-input` tag to the cells you want to hide.

### Step 4: Panel-Specific Rules

Create `workflow/rules/[panel-name].smk`:

```python
rule your_panel_analysis:
    input:
        nb=f"{workflow.basedir}/notebooks/{panel}/your-analysis.ipynb",
        vcf="results/vcfs/targets/{dataset}.annot.vcf",
        metadata="results/config/metadata.qcpass.tsv",
        panel_data=f"resources/{panel}/marker_SNPs.csv"
    output:
        nb="results/notebooks/{panel}/your-analysis.ipynb",
        docs_nb="docs/ampseeker-results/notebooks/{panel}/your-analysis.ipynb",
        results="results/{panel}/analysis_results.tsv"
    conda:
        "../envs/AmpSeeker-python.yaml"
    log:
        "logs/notebooks/{panel}/your-analysis.log"
    params:
        dataset=dataset,
        wkdir=wkdir
    shell:
        """
        papermill {input.nb} {output.nb} \
            -k AmpSeq_python \
            -p dataset {params.dataset} \
            -p metadata_path {input.metadata} \
            -p vcf_path {input.vcf} \
            -p panel_data_path {input.panel_data} 2> {log}
        
        cp {output.nb} {output.docs_nb}
        """
```

### Step 5: Integration with Main Workflow

Add panel detection logic to `workflow/Snakefile`. We only want to include this code if we are using this panel.

```python
# Add after existing panel includes
if panel == "your-panel":
    include: "rules/your-panel.smk"
```

Update `AmpSeekerOutputs()` function in `workflow/rules/common.smk`:

```python
# Add to AmpSeekerOutputs function
if config["panel"] == "your-panel":
    inputs.extend([
        "results/notebooks/your-panel/your-analysis.ipynb",
        "docs/ampseeker-results/notebooks/your-panel/your-analysis.ipynb",
    ])
```

Update the `workflow/notebooks/process-notebooks.ipynb` notebook to include processing of your new panel notebooks.
The processing ensures that injected papermill parameter cells are hidden in the final jupyter book.

```python
import glob

input_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/*ipynb")
if panel == 'ag-vampir':
    other_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/ag-vampir/*ipynb")
    input_nbs.extend(other_nbs)

# add this  
if panel == 'your-panel':
    your_panel_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/your-panel/*ipynb")
    input_nbs.extend(your_panel_nbs)
```


### Step 6: Configuration Template

Create an example configuration in `config/`. 
Do not add panel specific config options to the main `config/config.yaml` file, if additional config options are required, create a second config  (your-panel.config.yaml) and read this within the panel-specific analyses.

```yaml
# config/your-panel-example-config.yaml
dataset: your-panel-study
panel: your-panel
targets: config/your-panel.bed
metadata: config/your-panel-metadata.tsv

reference-fasta: resources/reference/genome.fasta
reference-gff3: resources/reference/annotations.gff
```

### Metadata Format

The metadata file must contain specific columns depending on the data source:

**BCL processing** (`from-bcl: True`):
Uses Illumina SampleSheet.csv format automatically.

**FASTQ processing** (`from-bcl: False`):
```tsv
sample_id    fq1                    fq2                    cohort_column1    ...
sample001    reads/sample001_1.fq   reads/sample001_2.fq   group_A          ...
```

**Nanopore**:
```tsv
sample_id    fq1                    cohort_column1    ...
sample001    reads/sample001.fq     group_A          ...
```

**Optional columns**:
- `plate`, `well_letter`, `well_number`: For plate layout visualization
- `latitude`, `longitude`: For geographic mapping
- Any number of cohort/grouping columns for analysis

## Development Guidelines

### Code Style

1. **Python**: Follow PEP 8 standards
2. **Plotting**: Use plotly for all plots and use consistent color schemes
3. **Snakemake**: Use consistent indentation and rule naming
4. **Notebooks**: Include markdown documentation for each analysis step
5. **Comments**: Document complex logic and parameter choices

### Adding New Analysis Modules

1. Create notebook in `workflow/notebooks/`
2. Implement as Snakemake rule in appropriate rule file
3. Add outputs to `AmpSeekerOutputs()` function
4. Test with example data
5. Update documentation

### Version Control

1. Use semantic versioning
2. Tag releases appropriately
3. Test thoroughly before merging to main

### Environment Management

All dependencies are managed through conda environments in `workflow/envs/`:

- `AmpSeeker-python.yaml`: Core Python analysis environment
- `AmpSeeker-cli.yaml`: Command-line tools (BWA, BCFtools, etc.)
- `AmpSeeker-qc.yaml`: Quality control tools
- `AmpSeeker-jupyterbook.yaml`: Report generation
- `AmpSeeker-snpeff.yaml`: Variant annotation

## Testing

### Test Structure

Tests are located in `.test/` directory with configurations for different scenarios:

```
.test/
├── config/
│   ├── config_agvampir.yaml      # Ag-vampIR panel test
│   ├── config_fastqauto.yaml     # Auto FASTQ detection test
│   └── config_fastqmetadata.yaml # Manual FASTQ paths test
└── resources/                    # Test data
```

### Running Tests

```bash
# Run specific test
snakemake --configfile .test/config/config_agvampir.yaml --use-conda -j4

# Run all tests (via GitHub Actions)
# Tests are automatically executed on pull requests
```

### Adding New Tests

1. Create test configuration in `.test/config/`
2. Download test data from external resource within the .github/workflow 
3. Ensure test completes successfully
4. Update GitHub Actions workflow if needed

## Contributing Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Implement changes following these guidelines
4. Add tests for new functionality
5. Update documentation as needed
6. Submit pull request with clear description
7. Address review comments
8. Merge after approval and testing

Thank you for contributing to AmpSeeker!