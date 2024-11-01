import os
import argparse
from Bio import SeqIO


def split_fasta_with_biopython(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for record in SeqIO.parse(input_file, "fasta"):
        output_file = os.path.join(output_dir, f"{record.id}.fna")

        with open(output_file, "w") as out_file:
            SeqIO.write(record, out_file, "fasta")

        print(f"Created file: {output_file}")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Split a multi-FASTA file into individual chromosome files.")
    parser.add_argument("input_file", help="Path to the input FASTA file containing multiple sequences.")
    parser.add_argument("output_dir", help="Directory to save the split FASTA files.")

    # Parse arguments
    args = parser.parse_args()

    # Call the function with command-line arguments
    split_fasta_with_biopython(args.input_file, args.output_dir)
