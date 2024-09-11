#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR LIB"
    exit 1
fi

# Get input directory, output directory, and species from arguments
INPUT_DIR="$1"
OUTPUT_DIR="$2"
LIB="$3"

# Loop through all .fasta or .fna files in the input directory
for fasta_file in ${INPUT_DIR}*.{fasta,fna}; do
  # Run RepeatMasker for each file with species specified
  RepeatMasker -lib "$LIB" -dir "$OUTPUT_DIR" "$fasta_file"
done

# Delete all files in OUTPUT_DIR that do not end in .out
echo "Cleaning up files in $OUTPUT_DIR that don't end in .out"
find "$OUTPUT_DIR" -type f ! -name "*.out" -exec rm {} +

echo "Cleanup completed."

# Example using -lib:
# RepeatMasker -lib  /opt/RepeatMasker/Libraries/famdb/dfam38_full.8.h5 /opt/RepeatMasker/Libraries/famdb/sequences/c_elegans/c_elegans_chromosome_I.fasta