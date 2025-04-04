#!/bin/bash

# Create necessary directories
mkdir -p .test/resources/reference
mkdir -p .test/config
mkdir -p .test/workflow
mkdir -p .test/docs

# Download and process reference genome with progress bar
# The progress bar will be shown on stderr while the content goes to stdout for the pipe
curl -L --progress-bar https://vectorbase.org/common/downloads/release-66/AgambiaePEST/fasta/data/VectorBase-66_AgambiaePEST_Genome.fasta | sed 's/AgamP4_//g' > .test/resources/reference/AgamP4.fa

# Copy configuration files
cp config/ag-vampir.bed .test/config/ag-vampir.bed

# Copy resource files
cp resources/*npy resources/*csv resources/*pickle .test/resources/.

# Copy documentation
cp -r docs/ .test/docs/

# Copy workflow tools
cp workflow/ampseekertools.py .test/workflow/

echo "Local test environment setup complete!"
