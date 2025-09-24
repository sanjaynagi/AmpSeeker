#!/bin/bash
set -euo pipefail

# Define target directory
TARGET_DIR="resources/models"

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

# Download the clair3 models archive
wget --quiet http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz

# Extract into the target directory
tar -zxvf clair3_models.tar.gz -C "$TARGET_DIR"

# Remove the tarball to save space
rm -f clair3_models.tar.gz

echo "Clair3 models have been installed in $TARGET_DIR"