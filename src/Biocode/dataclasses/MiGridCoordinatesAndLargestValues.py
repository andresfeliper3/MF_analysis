from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class MiGridCoordinatesAndLargestValues:
    k: int
    largest_values: List[int]
    coordinates: List[Tuple[int, int]]
    nucleotides_strings: List[str]

    def get_k(self) -> int:
        return self.k

    def get_largest_values(self) -> List[int]:
        return self.largest_values

    def get_coordinates(self) -> List[Tuple[int, int]]:
        return self.coordinates

    def set_nucleotides_strings(self, nucleotides_strings: List[str]):
        self.nucleotides_strings = nucleotides_strings

    def __str__(self):
        return f"MiGridCoordinatesAndLargestValues(k={self.k}, largest_values={self.largest_values}, coordinates={self.coordinates}, nucleotides_strings={self.nucleotides_strings})"