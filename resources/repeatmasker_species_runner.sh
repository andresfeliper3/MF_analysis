#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 INPUT_DIR OUTPUT_DIR SPECIES"
    exit 1
fi

# Get input directory, output directory, and species from arguments
INPUT_DIR="$1"
OUTPUT_DIR="$2"
SPECIES="$3"

# Loop through all .fasta or .fna files in the input directory
for fasta_file in ${INPUT_DIR}*.{fasta,fna}; do
  # Run RepeatMasker for each file with species specified
  RepeatMasker -species "$SPECIES" -dir "$OUTPUT_DIR" "$fasta_file"
  # Find and remove temporary directories
  find "$OUTPUT_DIR" -type d -name "RM_*" -exec rm -r {} +
  # Delete all files in OUTPUT_DIR that do not end in .out
  echo "Cleaning up files in $OUTPUT_DIR that don't end in .out"
  find "$OUTPUT_DIR" -type f ! -name "*.out" -exec rm {} +
  echo "Cleanup completed."
done


# Example using -species:
# RepeatMasker -species nematode -dir /working/results/c_elegans/CX /working/c_elegans/c_elegans_chromosome_X.fasta

# Example execution of the file
# ./repeatmasker_species_runner.sh /working/dna_sequences/Acropora_millepora/ /working/RM_resources/acropora_millepora/ cnidaria