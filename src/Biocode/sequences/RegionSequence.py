import numpy as np
from src.Biocode.sequences.Sequence import Sequence

from utils.logger import logger


class RegionSequence(Sequence):
    def __init__(self, sequence: str = None, sequence_data: dict = None, regions_number: int = 0, window_length: int = 0,
                 name: str = None, refseq_accession_number: str = None, organism_name: str = None):
        super().__init__(sequence=sequence, sequence_data=sequence_data, name=name,
                         refseq_accession_number=refseq_accession_number, organism_name= organism_name)
        self.region_common_size = None
        self.regions_number = regions_number #or len(sequence) // window_length
        self.window_length = window_length
        self.regions = []
        self.regions_names = []
        self.calculate_regions()

    def set_regions_number(self, regions_number):
        self.regions_number = regions_number
        self.calculate_regions()

    def calculate_regions(self):
        if (self.regions_number is None or self.regions_number < 0) and (self.window_length is None or self.window_length < 0):
            raise Exception('Enter a valid regions number or window length')
        elif self.regions_number is not None and self.regions_number > 0:
            self.region_common_size = self.get_size() // self.regions_number
            for i in range(self.regions_number):
                start = i * self.region_common_size
                end = (i + 1) * self.region_common_size if i < self.regions_number - 1 else None
                region = self.sequence[start:end]
                region_name = f"{self.name}_region_{i + 1}_of_{self.regions_number}"
                refseq_accession_number = f"{self.refseq_accession_number}_region_{i + 1}"
                region_sequence = Sequence(sequence=region, name=region_name,
                                             refseq_accession_number=refseq_accession_number,
                                             organism_name=self.organism_name)
                region_sequence.set_region_number(i + 1)
                region_sequence.set_regions_total(self.regions_number)
                self.regions.append(region_sequence)
                self.regions_names.append(region_name)
        elif self.window_length is not None and self.window_length > 0:
            regions_sequences = self._split_sequence_per_window_length(self.sequence, self.window_length)
            self.regions_number = len(regions_sequences)
            for i, region_sequence in enumerate(regions_sequences):
                region_name = f"{self.name}_region_{i + 1}_of_{self.regions_number}"
                refseq_accession_number = f"{self.refseq_accession_number}_region_{i + 1}"
                self.regions.append(Sequence(sequence=region_sequence, name=region_name,
                                             refseq_accession_number=refseq_accession_number,
                                             organism_name=self.organism_name))
                self.regions_names.append(region_name)

    def _split_sequence_per_window_length(self, sequence, window_length):
        return [sequence[i:i + window_length] for i in range(0, len(sequence), window_length)]

    def get_regions(self) -> list[Sequence]:
        return self.regions

    def get_regions_number(self) -> int:
        return self.regions_number

    def get_regions_sizes(self) -> list[int]:
        return np.array([len(region.get_sequence()) for region in self.regions])

    def get_regions_names(self) -> list[str]:
        return self.regions_names

    def get_window_length(self) -> int:
        return self.window_length
