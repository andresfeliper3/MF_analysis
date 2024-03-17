from Bio import SeqIO


class Sequence:
    def __init__(self, sequence: str = None, sequence_data: dict = None, name: str = None):
        if sequence:
            self.sequence = sequence
            self.name = name
        elif sequence_data:
            self.sequence = str(SeqIO.read(sequence_data['path'], "fasta").seq)
            self.name = sequence_data['name']
        self.size = len(self.sequence)
        self.cover = []
        self.cover_percentage = 0.0

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
