import os
import sys

project_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..'))  # Get the absolute path of the parent directory
sys.path.append(project_root)  # Add the project root to the Python path

# from local_scripts.load import load_organism, whole_MFA
import argparse
import yaml

# Load configurations from the YAML file


sequences_file_path = os.path.join(project_root, 'src/sequences.yaml')

with open(sequences_file_path, 'r') as sequences_file:
    config = yaml.safe_load(sequences_file)


# Function to read fasta sequence
def read_fasta_sequence(file_path):
    sequence = ""

    with open(file_path, "r") as file:
        lines = file.readlines()

        # Skip header lines (lines starting with '>')
        sequence_lines = [line.strip() for line in lines if not line.startswith(">")]

        # Concatenate the lines to form the sequence
        sequence = "".join(sequence_lines)

    return sequence


def create_sequence_data_dict(path):
    if not os.path.exists(path):
        os.makedirs(path)

    files = os.listdir(path)
    sorted_files = sorted(files,
                          key=lambda x: int(x.rstrip('.fna')[3:]) if x.rstrip('.fna')[3:].isdigit() else float('inf'))

    return [
        {"path": os.path.join(path, file), "name": file.split(".")[0]}
        for file in sorted_files
    ]


organism = ""


def analyze_command(args):
    global organism
    if args.id:
        organism = args.jd

    elif args.name:
        organism = args.name
    else:
        print("Please provide either -id or -name.")


def graph_command(args):
    if args.id:
        print("Graphing by ID:", args.id)
    elif args.name:
        print("Graphing by name:", args.name)
    else:
        print("Please provide either -id or -name.")


def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_parser = subparsers.add_parser('analyze', help='Analyze command')
    analyze_parser.add_argument('-id', help='ID for analysis')
    analyze_parser.add_argument('-name', help='Name for analysis')

    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-id', help='ID for graphing')
    graph_parser.add_argument('-name', help='Name for graphing')

    args = parser.parse_args()

    if args.command == 'analyze':
        analyze_command(args)
    elif args.command == 'graph':
        graph_command(args)


if __name__ == "__main__":
    main()

sequences_folder = os.path.join(os.path.dirname(__file__), config["sequences_folder"])
ORGANISM_NAME = config[organism]['organism_name']
GCF = config[organism]['GCF']
REGIONS_NUMBER = config[organism]['regions_number']
ORGANISM_FOLDER = ORGANISM_NAME.replace(" ", "_")
DOWNLOAD_URL = config[organism]['download_url']

organism_path = os.path.abspath(os.path.join(sequences_folder, ORGANISM_FOLDER))

data = create_sequence_data_dict(organism_path)
AMOUNT_CHROMOSOMES = len(data)

# change sequences.yaml to include all sequences
# make them be selectable by name or GCF
# separation of responsibilities. load.py includes loading and commanding.
# after executing command, the load should be made and then analyze
