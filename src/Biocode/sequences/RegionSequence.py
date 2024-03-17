import numpy as np
from src.Biocode.sequences.Sequence import Sequence


class RegionSequence(Sequence):
    def __init__(self, sequence: str = None, sequence_data: dict = None, regions_number: int = 0, name: str = None):
        super().__init__(sequence=sequence, sequence_data=sequence_data, name=name)
        self.region_common_size = None
        self.regions_number = regions_number
        self.regions = []
        self.regions_names = []
        self.calculate_regions()

    def set_regions_number(self, regions_number):
        self.regions_number = regions_number
        self.calculate_regions()

    def calculate_regions(self):
        if self.regions_number is None or self.regions_number <= 0:
            raise Exception('Enter a valid regions number')
        else:
            self.region_common_size = self.get_size() // self.regions_number
            for i in range(self.regions_number):
                start = i * self.region_common_size
                end = (i + 1) * self.region_common_size if i < self.regions_number - 1 else None
                region = self.sequence[start:end]
                region_name = f"{self.name}_region_{i + 1}"
                self.regions.append(Sequence(sequence=region, name=region_name))
                self.regions_names.append(region_name)

    def get_regions(self) -> list[Sequence]:
        return self.regions

    def get_regions_number(self) -> int:
        return self.regions_number

    def get_regions_sizes(self) -> list[int]:
        return np.array([len(region.get_sequence()) for region in self.regions])

    def get_regions_names(self) -> list[str]:
        return self.regions_names
