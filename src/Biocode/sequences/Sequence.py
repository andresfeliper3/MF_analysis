from Bio import SeqIO


class Sequence:
    def __init__(self, sequence: str = None, sequence_data: dict = None, name: str = None,
                 refseq_accession_number: str = None, organism_name: str = None):
        if sequence:
            self.sequence = sequence
            self.name = name
            self.refseq_accession_number = refseq_accession_number
        elif sequence_data:
            self.sequence = str(SeqIO.read(sequence_data['path'], "fasta").seq)
            self.name = sequence_data['name']
            self.refseq_accession_number = sequence_data['refseq_accession_number']
        self.size = len(self.sequence)
        self.cover = []
        self.cover_percentage = 0.0
        self.organism_name = organism_name
        self.region_number = 1
        self.regions_total = 0

    def get_sequence(self) -> str:
        return self.sequence

    def get_size(self) -> int:
        return self.size

    def set_name(self, name: str):
        self.name = name

    def get_name(self) -> str:
        return self.name

    def set_cover(self, cover: list[int]):
        self.cover = cover

    def get_cover(self) -> list[int]:
        return self.cover

    def set_cover_percentage(self, cover_percentage: float):
        self.cover_percentage = cover_percentage

    def get_cover_percentage(self) -> float:
        return self.cover_percentage

    def get_refseq_accession_number(self) -> str:
        return self.refseq_accession_number

    def get_organism_name(self) -> str:
        return self.organism_name

    def set_organism_name(self, organism_name):
        self.organism_name = organism_name

    def set_region_number(self, region_number):
        self.region_number = region_number

    def set_regions_total(self, regions_total):
        self.regions_total = regions_total

    def get_region_number(self):
        return self.region_number

    def get_regions_total(self):
        return self.regions_total