import os
import yaml
from utils.logger import logger


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

    def create_sequence_data_dict(self, path) -> list:
        if not os.path.exists(path):
            os.makedirs(path)

        files = os.listdir(path)
        sorted_files = sorted(files,
                              key=lambda x: int(x.rstrip('.fna')[3:]) if x.rstrip('.fna')[3:].isdigit() else float(
                                  'inf'))

        return [
            {"path": os.path.join(path, file), "name": file.split(".")[0]}
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

    def get_organism_folder(self) -> str:
        try:
            return self.get_organism_name().replace(" ", "_")
        except Exception as e:
            logger.error("Invalid data entered")
            return None

    def get_download_url(self) -> str:
        try:
            return self.config[self.organism]['download_url']
        except KeyError as e:
            logger.error("Invalid download url %s", e)
            return ""

    def get_data(self) -> dict:
        organism_path = f"{self.get_sequences_folder()}/{self.get_organism_folder()}"
        data = self.create_sequence_data_dict(organism_path)
        return data

    def get_amount_chromosomes(self) -> int:
        return len(self.get_data())


loader = Loader()

