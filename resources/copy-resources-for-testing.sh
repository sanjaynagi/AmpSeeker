#!/bin/bash

# Create necessary directories
mkdir -p tests/snakemake/resources/reference
mkdir -p tests/snakemake/resources/ag-vampir
mkdir -p tests/snakemake/config
mkdir -p tests/snakemake/workflow/lib
mkdir -p tests/snakemake/docs

# Download and process reference genome with progress bar
# The progress bar will be shown on stderr while the content goes to stdout for the pipe
curl -L --progress-bar https://vectorbase.org/common/downloads/release-66/AgambiaePEST/fasta/data/VectorBase-66_AgambiaePEST_Genome.fasta | sed 's/AgamP4_//g' > tests/snakemake/resources/reference/AgamP4.fa

# Copy configuration files
cp config/ag-vampir.bed tests/snakemake/config/ag-vampir.bed

# Copy resource files
cp resources/ag-vampir/* tests/snakemake/resources/ag-vampir/.

# Copy documentation
cp -r docs/ tests/snakemake/.

# Copy workflow tools
cp -r workflow/lib tests/snakemake/workflow/

echo "Local test environment setup complete!"
