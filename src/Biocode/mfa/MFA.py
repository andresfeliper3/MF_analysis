import numpy as np
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.mfa.CGR import CGR
from utils.logger import logger

from src.Biocode.graphs.Graphs import Graphs


class MFA:
    def __init__(self, sequence: Sequence):
        self._10_largest_mi_grid_values_for_k_from_10_to_4 = None
        self.DDq = None
        self.Dqmin = None
        self.Dqmax = None
        self.fq = None
        self.total_fractal_points = None
        self.sequence = sequence
        self.q_min = -20
        self.q_max = 20
        self.q_values = np.arange(self.q_min, self.q_max + 1)
        self.Dq_values = np.zeros(len(self.q_values))
        self.tau_q_values = np.zeros(len(self.q_values))
        self.FIX_1_ERROR_VALUE = 0.99999

        self.result = {}

        self.GRID_EXPONENTS = np.array([10, 9, 8, 7, 6, 5, 4, 3, 2])
        self.sizes = np.power(2, self.GRID_EXPONENTS)
        self.epsilons = 1 / self.sizes

    def _generate_cgr_mi_grids_thoroughly(self):  # NOT BEING USED
        self.cgr_gen = CGR(self.sequence)
        self.cgrs_mi_grids = []

        for index, epsilon in enumerate(self.epsilons):
            cgr_mi_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=False, epsilon=epsilon)
            self.cgrs_mi_grids.append(cgr_mi_grid)
            logger.info(f"Ready mi_grid with size {self.sizes[index]}, and epsilon {epsilon}")

    def _generate_cgr_mi_grids_quickly(self):
        self.cgr_gen = CGR(self.sequence)
        self.cgrs_mi_grids = [0] * len(self.sizes)

        cgr_largest_mi_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=False, epsilon=self.epsilons[0])
        self.cgrs_mi_grids[0] = cgr_largest_mi_grid

        logger.info(f"Ready mi_grid with size {self.sizes[0]} and epsilon {self.epsilons[0]}")

        for i in range(1, len(self.sizes)):
            self.cgrs_mi_grids[i] = self.__resize_matrix(original_matrix=self.cgrs_mi_grids[i - 1],
                                                         target_size=self.sizes[i])
            logger.info(f"Ready mi_grid with size {self.sizes[i]} and epsilon {self.epsilons[i]}")

    def __resize_matrix(self, original_matrix, target_size):
        result_matrix = [[0 for _ in range(target_size)] for _ in range(target_size)]

        for i in range(target_size):
            for j in range(target_size):
                result_matrix[i][j] = (
                        original_matrix[2 * i][2 * j] +
                        original_matrix[2 * i + 1][2 * j] +
                        original_matrix[2 * i][2 * j + 1] +
                        original_matrix[2 * i + 1][2 * j + 1]
                )
        return np.array(result_matrix)

    def multifractal_discrimination_analysis(self):
        self._generate_cgr_mi_grids_quickly()
        self.total_fractal_points = np.sum(self.cgrs_mi_grids[0])
        assert np.sum(self.cgrs_mi_grids[1]) == self.total_fractal_points, "total fractal points should be the same"
        self.fq = []

        for q_index, q in enumerate(self.q_values):
            if q == 1:
                q = self.FIX_1_ERROR_VALUE
            self.fq.append({'q': q, 'fq': np.zeros(len(self.cgrs_mi_grids))})
            for index, mi_grid in enumerate(self.cgrs_mi_grids):
                cgr_mi_grid_flattened = mi_grid.reshape(-1)
                no_zeros = cgr_mi_grid_flattened[np.where(cgr_mi_grid_flattened != 0)]
                division = no_zeros / self.total_fractal_points
                powered = np.power(division, q)
                sum_term = np.sum(powered, dtype=np.float64)
                numerator = np.log(sum_term)
                denominator = (q - 1)
                self.fq[-1]['fq'][index] = numerator / denominator

            logger.debug(self.fq)
            linear_coefficients = np.polyfit(np.log(self.epsilons), self.fq[-1]['fq'], 1)

            # CHECK THIS OUT
            #Graphs.graph_linear_regression(self.fq[-1], self.epsilons, self.sequence.get_name(), "Linear_regression_test")


            Dq = linear_coefficients[0]  # slope of the linear function
            self.Dq_values[q_index] = Dq

            tau_q = (q - 1) * Dq
            self.tau_q_values[q_index] = tau_q


        # Find the maximum and minimum Dq values
        self.Dqmax = np.max(self.Dq_values)
        self.Dqmin = np.min(self.Dq_values)

        # Calculate the degree of multifractality (DDq)
        self.DDq = self.Dqmax - self.Dqmin

        # plt.show()
        print("finished chromosome:", self.sequence.get_name())
        print("*************************************")
        self.result = {
            'q_values': self.q_values,
            'Dq_values': self.Dq_values,
            'tau_q_values': self.tau_q_values,
            'DDq': self.DDq,
            'sequence_name': self.sequence.get_name()
        }
        return self.result

    def print_tab_values(self):
        # Print or visualize the calculated Dq and τ(q) values
        print("q\tDq\tτ(q)")
        for q, Dq, tau_q in zip(self.q_values, self.Dq_values, self.tau_q_values):
            print(f"{q}\t{Dq}\t{tau_q}")

    def print_degree_of_multifractality(self):
        # Print the DDq value
        print(f"Degree of Multifractality (DDq): {self.DDq}")

    def graph_multifractal_spectrum(self):
        # Dq vs q
        Graphs.graph_one(x_array=self.q_values, y_array=self.Dq_values, x_label='q', y_label='Dq',
                         title='Multifractal spectrum (Dq vs q)', name=self.sequence.get_name())

    def graph_correlation_exponent(self):
        # t(q) vs q
        Graphs.graph_one(x_array=self.q_values, y_array=self.tau_q_values, x_label='q', y_label='t(q)',
                         title='Correlation exponent (t(q) vs q)', name=self.sequence.get_name())

    def get_q_values(self):
        return self.q_values

    def get_Dq_values(self):
        return self.Dq_values

    def get_tau_q_values(self):
        return self.tau_q_values

    def get_DDq(self) -> float:
        return self.DDq

    def get_Dqmax(self) -> float:
        return self.Dqmax

    def get_Dqmin(self) -> float:
        return self.Dqmin

    def get_result(self) -> dict:
        return self.result

    def get_total_fractal_points(self) -> float:
        return self.total_fractal_points

    def get_cgrs_mi_grids(self) -> list:
        return self.cgrs_mi_grids

    def get_epsilons(self) -> list[float]:
        return self.epsilons

    def get_sizes(self) -> list:
        return self.sizes

    def get_q_min(self) -> int:
        return self.q_min

    def get_q_max(self) -> int:
        return self.q_max

    def get_cgr_gen(self) -> CGR:
        return self.cgr_gen

    def get_fq(self) -> list[dict]:
        return self.fq

    def get_largest_mi_grid_values_given_k_n(self, k: int, n: int) -> list[int]:
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
        selected_mi_grid_flattened_sorted_desc = np.sort(selected_mi_grid_flattened)[::-1]
        largest_n_values = selected_mi_grid_flattened_sorted_desc[:n]
        return largest_n_values

    def get_10_largest_mi_grid_values_for_k_from_10_to_4(self):
        self._10_largest_mi_grid_values_for_k_from_10_to_4 = []
        for k in range(10, 3, -1):
            self._10_largest_mi_grid_values_for_k_from_10_to_4.append(
                self.get_largest_mi_grid_values_given_k_n(k=k, n=10))

        self._10_largest_mi_grid_values_for_k_from_10_to_4 = np.array(self._10_largest_mi_grid_values_for_k_from_10_to_4)
        return self._10_largest_mi_grid_values_for_k_from_10_to_4

