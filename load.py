import os
import yaml
from utils.logger import logger
from Bio import SeqIO



class Loader:
    def __init__(self):
        self.organism = None
        self.config = None
        self.sequences_file_path = 'resources/sequences.yaml'
        self.load_yaml_file()

    def load_yaml_file(self):
        with open(self.sequences_file_path, 'r') as sequences_file:
            self.config = yaml.safe_load(sequences_file)

    # Function to read fasta sequence
    def read_fasta_sequence(self, file_path):
        sequence = ""

        with open(file_path, "r") as file:
            lines = file.readlines()

            # Skip header lines (lines starting with '>')
            sequence_lines = [line.strip() for line in lines if not line.startswith(">")]

            # Concatenate the lines to form the sequence
            sequence = "".join(sequence_lines)

        return sequence

    def extract_refseq_accession_number(self, file_path) -> str:
        with open(file_path, 'r') as fna_file:
            for record in SeqIO.parse(fna_file, 'fasta'):
                # Split the description by space and check the first part
                parts = record.description.split()
                first_part = parts[0]

                if first_part.startswith("gi|"):
                    # If it starts with 'gi|', extract the third part (refseq accession number)
                    refseq_accession_number = first_part.split('|')[3]
                elif first_part.startswith("NC_") or first_part.startswith("XM_") or first_part.startswith("XP_"):
                    # If it starts directly with the accession number (e.g., NC_), use it as is
                    refseq_accession_number = first_part
                else:
                    # Handle any unexpected format if needed
                    logger.error(f"Unexpected format in sequence description: {record.description}")
                    continue

                return refseq_accession_number

    def extract_file_name(self, file_path) -> str:
        return os.path.basename(file_path).split(".")[0]

    def create_sequence_data_dict(self, path) -> list:
        if not os.path.exists(path):
            os.makedirs(path)

        files = os.listdir(path)
        sorted_files = sorted(files,
                              key=lambda x: int(x.rstrip('.fna')[3:]) if x.rstrip('.fna')[3:].isdigit() else float(
                                  'inf'))

        return [
            {"path": os.path.join(path, file), "name": file.split(".")[0], "organism_name": self.organism}
            for file in sorted_files
        ]

    def set_organism(self, organism: str):
        self.organism = organism

    def get_organism(self) -> str:
        return self.organism

    def get_sequences_folder(self) -> str:
        return self.config["sequences_folder"]

    def get_organism_name(self) -> str:
        try:
            return self.config[self.organism]['organism_name']
        except Exception as e:
            logger.error("Invalid organism name: %s", e)
            return None

    def get_gcf(self) -> str:
        try:
            return self.config[self.organism]['GCF']
        except KeyError as e:
            logger.error("Invalid GCF: %s", e)
            return ""

    def get_regions_number(self) -> int:
        try:
            return self.config[self.organism]['regions_number']
        except KeyError as e:
            logger.error("Invalid regions number: %s", e)
            return -1

    def get_organism_folder(self) -> str | None:
        try:
            return self.get_organism_name().replace(" ", "_")
        except Exception as e:
            logger.error("Invalid data entered")
            return None

    def get_organism_gtf_subfolder(self) -> str:
        return f"{self.get_organism_folder()}/gtf"

    def get_download_url(self) -> str:
        try:
            return self.config[self.organism]['download_url']
        except KeyError as e:
            logger.error("Invalid download url %s", e)
            return ""

    def get_download_gff_url(self) -> str:
        try:
            return self.config[self.organism]['download_files_url']
        except KeyError as e:
            logger.error("Invalid download url %s", e)
            return ""

    def get_data(self) -> list:
        organism_path = f"{self.get_sequences_folder()}/{self.get_organism_folder()}"
        data = self.create_sequence_data_dict(organism_path)
        return data

    def get_amount_chromosomes(self) -> int:
        return len(self.get_data())


loader = Loader()

