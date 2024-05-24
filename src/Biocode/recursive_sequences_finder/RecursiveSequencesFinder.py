from typing import Tuple, Any, List
import numpy as np

from src.Biocode.dataclasses.MiGridCoordinatesValuesAndNucleotides import MiGridCoordinatesValuesAndNucleotides


class RecursiveSequencesFinder:


    def __init__(self, grid_exponents: List[int], cgrs_mi_grids: List[List[List[int]]]):
        self.largest_mi_grid_values_for_different_k = None
        self.GRID_EXPONENTS = grid_exponents
        self.cgrs_mi_grids = cgrs_mi_grids



    def get_largest_mi_grid_values_given_k_n_and_coordinates(self, k: int, n: int) -> Tuple[Any, List]:
        """
        get_largest_mi_grid_value_given_k_n
        :param k: exponent of the grid size. E.g. if largest values with grid size 1024 is desired, k must be 10.
        :param n: amount of largest values desired. If the 10 max values are desired, n must be 10.
        :return: array of size n with the largest values using the grid size 2^k
        """
        if k not in self.GRID_EXPONENTS:
            raise Exception("Invalid k value for mfa.get_largest_mi_grid_value_given_k_n")

        exponent_index = np.where(self.GRID_EXPONENTS == k)[0][0] # to access numeric value
        selected_mi_grid = self.cgrs_mi_grids[exponent_index]
        selected_mi_grid_flattened = selected_mi_grid.flatten()

        # Get the indices of the n largest values in the flattened array
        indices_of_largest_n_values = np.argpartition(-selected_mi_grid_flattened, n)[:n]
        # Get the n largest values
        largest_n_values = selected_mi_grid_flattened[indices_of_largest_n_values]

        coordinates = self._find_coordinates_of_values_in_matrix_flattened(indices_of_largest_n_values, (2**k, 2**k))
        return largest_n_values, coordinates

    def _find_coordinates_of_values_in_matrix_flattened(self, indices: List, matrix_shape: Tuple) -> List:
        return [np.unravel_index(index, matrix_shape) for index in indices]


    def get_largest_mi_grid_values_strings_for_different_k(self, k1:int, k2: int, k_step: int, amount_sequences: int) -> List[MiGridCoordinatesValuesAndNucleotides]:
        self.largest_mi_grid_values_for_different_k = []
        for k in range(k1, k2, k_step):
            largest_values, coordinates =  self.get_largest_mi_grid_values_given_k_n_and_coordinates(k=k, n=amount_sequences)

            self.largest_mi_grid_values_for_different_k.append(MiGridCoordinatesValuesAndNucleotides(
                k=k, largest_values=largest_values, coordinates=coordinates, nucleotides_strings=[]))

        self._find_nucleotide_strings_from_coordinates()
        self.largest_mi_grid_values_for_different_k = np.array(self.largest_mi_grid_values_for_different_k)
        return self.largest_mi_grid_values_for_different_k

    def _find_nucleotide_strings_from_coordinates(self, coordinates_list: List[MiGridCoordinatesValuesAndNucleotides]=None):
        if self.largest_mi_grid_values_for_different_k is None:
            self.largest_mi_grid_values_for_different_k = coordinates_list


        for values_obj in self.largest_mi_grid_values_for_different_k:
            grid_size = 2 ** values_obj.get_k()
            n_largest_nucleotides_strings_list = []
            for coordinates in values_obj.get_coordinates():
                nucleotide_string_reversed = self._find_nucleotide_string_from_coordinate_reversed(
                    coordinates=coordinates, grid_size=grid_size, nucleotide_string="")

                nucleotide_string = nucleotide_string_reversed[::-1]
                n_largest_nucleotides_strings_list.append(nucleotide_string)

            values_obj.set_nucleotides_strings(n_largest_nucleotides_strings_list)




    def _find_nucleotide_string_from_coordinate_reversed(self, coordinates: Tuple, grid_size: int,
                                                         nucleotide_string: str) -> str:
        row = coordinates[0]
        col = coordinates[1]
        half_grid_index = (grid_size / 2) - 1

        is_row_lower_than_half = row <= half_grid_index
        is_col_lower_than_half = col <= half_grid_index

        if grid_size == 1:
            return nucleotide_string
        elif is_row_lower_than_half and is_col_lower_than_half:
            return self._find_nucleotide_string_from_coordinate_reversed(coordinates, grid_size // 2, nucleotide_string + "A")

        elif is_row_lower_than_half and not is_col_lower_than_half:
            return self._find_nucleotide_string_from_coordinate_reversed((row, col - half_grid_index - 1), grid_size // 2,
                                                                nucleotide_string + "C")
        elif not is_row_lower_than_half and not is_col_lower_than_half:
            return self._find_nucleotide_string_from_coordinate_reversed((row - half_grid_index - 1, col - half_grid_index - 1),
                                                                grid_size // 2, nucleotide_string + "G")
        elif not is_row_lower_than_half and is_col_lower_than_half:
            return self._find_nucleotide_string_from_coordinate_reversed((row - half_grid_index - 1, col), grid_size // 2,
                                                                nucleotide_string + "T")


