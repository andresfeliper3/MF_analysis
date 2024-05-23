from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class MiGridCoordinatesAndLargestValues:
    k: int
    largest_values: List[int]
    coordinates: List[Tuple[int, int]]

    def get_k(self) -> int:
        return self.k

    def get_largest_values(self) -> List[int]:
        return self.largest_values

    def get_coordinates(self) -> List[Tuple[int, int]]:
        return self.coordinates