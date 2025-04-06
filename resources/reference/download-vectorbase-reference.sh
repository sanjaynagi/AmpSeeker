#!/bin/bash

# Default values
GENOME_NAME="AgambiaePEST"
VERSION="66"
STRIP_STRING="AgamP4_"
OUTPUT_DIR="./resources/reference"

# Function to display usage information
usage() {
    echo "Usage: $0 [GENOME_NAME] [VERSION] [STRIP_STRING]"
    echo "  GENOME_NAME   : Name of the genome (default: AgambiaePEST)"
    echo "  VERSION       : Version number (default: 66)"
    echo "  STRIP_STRING  : String to strip from sequence headers (default: AgamP4_)"
    echo ""
    echo "Example: $0 AgambiaePEST 66 AgamP4_"
    exit 1
}

# Parse arguments
[ $# -ge 1 ] && GENOME_NAME=$1
[ $# -ge 2 ] && VERSION=$2
[ $# -ge 3 ] && STRIP_STRING=$3

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Extract the short name from the STRIP_STRING (removing trailing underscore if present)
SHORT_NAME=${STRIP_STRING%_}
[ -z "$SHORT_NAME" ] && SHORT_NAME=$GENOME_NAME

echo "Downloading genome and annotation files for $GENOME_NAME (version $VERSION)"
echo "Will strip '$STRIP_STRING' from sequence headers"

# Download and process the FASTA file
echo "Downloading FASTA file..."
FASTA_URL="https://vectorbase.org/common/downloads/release-$VERSION/$GENOME_NAME/fasta/data/VectorBase-$VERSION"_"$GENOME_NAME"_"Genome.fasta"
curl -L --progress-bar "$FASTA_URL" | sed "s/$STRIP_STRING//g" > "$OUTPUT_DIR/$SHORT_NAME.fa"

# Download and process the GFF file
echo "Downloading GFF file..."
GFF_URL="https://vectorbase.org/common/downloads/release-$VERSION/$GENOME_NAME/gff/data/VectorBase-$VERSION"_"$GENOME_NAME.gff"
curl -L --progress-bar "$GFF_URL" | sed "s/$STRIP_STRING//g" > "$OUTPUT_DIR/$SHORT_NAME.gff"

echo "Download complete. Files saved to $OUTPUT_DIR/"
echo "  - $OUTPUT_DIR/$SHORT_NAME.fa"
echo "  - $OUTPUT_DIR/$SHORT_NAME.gff"
