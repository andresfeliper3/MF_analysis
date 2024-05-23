from dataclasses import dataclass

@dataclass
class MiGridCoordinatesAndLargestValues:
    k: int
    largest_values: list[int]
    coordinates: list[tuple]

