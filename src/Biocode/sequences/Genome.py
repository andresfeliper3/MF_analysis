from Bio import SeqIO
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.sequences.RegionSequence import RegionSequence

from Loader import Loader
from utils.logger import logger

from utils.decorators import Inject


@Inject(loader = Loader)
class Genome:
    def __init__(self, chromosomes: list[Sequence] = None, chromosomes_data: list[dict] = None,
                 regions_number: int = 0, window_length: int = 0, loader: Loader = None):
        self.loader = loader
        if (regions_number is None or regions_number < 0) and (
                window_length is None or window_length < 0):
            raise Exception('Enter a valid regions_1 number or window length')
        elif regions_number == 0 and window_length == 0:
            if chromosomes:
                self.chromosomes = chromosomes
            elif chromosomes_data:
                self.chromosomes = []
                for data in chromosomes_data:
                    self.chromosomes.append(Sequence(str(SeqIO.read(data['path'], "fasta").seq), name=data['name'],
                                        organism_name=data['organism_name'],
                                        refseq_accession_number=self.loader.extract_refseq_accession_number(data['path'])))
        else:  # with regions_1 number or window length
            if chromosomes:
                self.chromosomes = [RegionSequence(sequence=chromosome.get_sequence(), regions_number=regions_number,
                                                   refseq_accession_number=chromosome.get_refseq_accession_number(),
                                                   organism_name=chromosome.get_organism_name(), window_length=window_length)
                                    for chromosome in chromosomes]
            elif chromosomes_data:
                self.chromosomes = []
                for data in chromosomes_data:
                    self.chromosomes.append(
                        RegionSequence(sequence=str(SeqIO.read(data['path'], "fasta").seq), name=data['name'],
                                       regions_number=regions_number,
                                       refseq_accession_number=self.loader.extract_refseq_accession_number(data['path']),
                                       organism_name=data['organism_name'], window_length=window_length))

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
