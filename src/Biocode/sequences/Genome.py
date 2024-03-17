from Bio import SeqIO
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.sequences.RegionSequence import RegionSequence


class Genome:
    def __init__(self, chromosomes: list[Sequence] = None, chromosomes_data: list[dict] = None,
                 regions_number: int = 0):
        if regions_number < 0:
            raise Exception("Enter a valid regions_number for the Genome")
        elif regions_number == 0:
            if chromosomes:
                self.chromosomes = chromosomes
            elif chromosomes_data:
                self.chromosomes = []
                for data in chromosomes_data:
                    self.chromosomes.append(Sequence(str(SeqIO.read(data['path'], "fasta").seq), name=data['name']))
        else:  # with regions number
            if chromosomes:
                self.chromosomes = [RegionSequence(sequence=chromosome.get_sequence(), regions_number=regions_number)
                                    for chromosome in chromosomes]
            elif chromosomes_data:
                self.chromosomes = []
                for data in chromosomes_data:
                    self.chromosomes.append(
                        RegionSequence(sequence=str(SeqIO.read(data['path'], "fasta").seq), name=data['name'],
                                       regions_number=regions_number))

        self.number_of_chromosomes = len(self.chromosomes)
        self.chromosomes_names = [chromosome.get_name() for chromosome in self.chromosomes]

    def add_chromosome(self, chromosome):
        self.chromosomes.append(chromosome)

    def get_number_of_chromosomes(self) -> int:
        return self.number_of_chromosomes

    def get_chromosomes(self) -> list[Sequence]:
        return self.chromosomes

    def _calculate_chromosomes_sizes(self):
        self.chromosomes_sizes = [chromosome.get_size() for chromosome in self.chromosomes]

    def get_chromosomes_sizes(self) -> list[int]:
        self._calculate_chromosomes_sizes()
        return self.chromosomes_sizes

    def get_genome_size(self) -> int:
        self._calculate_chromosomes_sizes()
        return sum(self.chromosomes_sizes)

    def get_chromosomes_names(self) -> list[str]:
        return self.chromosomes_names
