import random
from src.Biocode.sequences.Sequence import Sequence


class RandomSequence(Sequence):

    def __init__(self, size: int):
        super().__init__("")
        self.nucleotides = ['A', 'T', 'C', 'G']
        self.size = size
        self.refresh()

    def _generate_random_nucleotide(self):
        random_number = random.randint(0, 3)
        random_nucleotide = self.nucleotides[random_number]
        return random_nucleotide

    def refresh(self):
        seq = ""
        for i in range(self.size):
            seq += self._generate_random_nucleotide()
        self.sequence = seq
