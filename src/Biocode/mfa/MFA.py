import numpy as np
import pickle
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.services.WholeMiGridsService import WholeMiGridsService

from utils.decorators import Timer, DBConnection

from src.Biocode.mfa.CGR import CGR
from utils.logger import logger

from src.Biocode.graphs.Graphs import Graphs

import gc

class MFA:
    def __init__(self, sequence: Sequence):
        self.cgr_gen = None
        self.cgr_initial_grid = None
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
        self.grid_sizes = np.power(2, self.GRID_EXPONENTS)
        self.epsilons = 1 / self.grid_sizes

    def _generate_cgr_mi_grids_thoroughly(self):  # NOT BEING USED
        self.cgr_gen = CGR(self.sequence)
        cgrs_mi_grids = []

        for index, epsilon in enumerate(self.epsilons):
            cgr_mi_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=False, epsilon=epsilon)
            cgrs_mi_grids.append(cgr_mi_grid)
            logger.info(f"Ready mi_grid with size {self.grid_sizes[index]}, and epsilon {epsilon}")
        return cgrs_mi_grids

    @DBConnection
    @Timer
    def _generate_cgr_mi_grids_quickly_using_database(self):
        self.cgr_gen = CGR(self.sequence)
        cgr_largest_mi_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=False, epsilon=self.epsilons[0])

        cgr_largest_mi_grid_np = np.array(cgr_largest_mi_grid, dtype=np.float64)
        del cgr_largest_mi_grid

        # Option 1: Using pickle for serialization
        binary_data = pickle.dumps(cgr_largest_mi_grid_np)
        # Option 2: Using NumPy's tobytes() method
        binary_data_numpy = cgr_largest_mi_grid_np.tobytes()
        del cgr_largest_mi_grid_np

        mi_grids_service = WholeMiGridsService()
        mi_grids_service.insert(record=(binary_data_numpy, 1, self.epsilons[0]))
        del binary_data_numpy
        ##
        gc.collect()
        retrieved_data = mi_grids_service.extract_mi_grid_by_chromosome_id(1)
        # Deserialize the data - use the same method used for storing
        # If using pickle:
        #cgr_largest_mi_grid = pickle.loads(retrieved_data)
        # If using NumPy's tobytes():
        cgr_largest_mi_grid = np.frombuffer(retrieved_data, dtype=np.float64)

        # Reshape the matrix
        cgr_largest_mi_grid = cgr_largest_mi_grid.reshape((self.grid_sizes[0], self.grid_sizes[0]))
        return self.generate_cgr_mi_grids_from_initial_grid(cgr_largest_mi_grid)

    def generate_initial_grid(self, name, title):
        self.cgr_gen = CGR(self.sequence)
        self.cgr_initial_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=True, epsilon=self.epsilons[0],
                                                                              name=name, title=title)
        return self.cgr_initial_grid

    def _generate_cgr_mi_grids_quickly(self) -> list:
        self.cgr_gen = CGR(self.sequence)
        cgr_largest_mi_grid = self.cgr_gen.generate_cgr_counting_grid_cells(graph=False, epsilon=self.epsilons[0])
        return self.generate_cgr_mi_grids_from_initial_grid(cgr_largest_mi_grid)

    def generate_cgr_mi_grids_from_initial_grid(self, cgr_largest_mi_grid: list) -> list:
        cgrs_mi_grids = [0] * len(self.grid_sizes)
        cgrs_mi_grids[0] = cgr_largest_mi_grid
        logger.info(f"Ready largest_mi_grid with size {self.grid_sizes[0]} and epsilon {self.epsilons[0]}")
        for i in range(1, len(self.grid_sizes)):
            cgrs_mi_grids[i] = self.__resize_matrix(original_matrix=cgrs_mi_grids[i - 1],
                                                    target_size=self.grid_sizes[i])
            logger.info(f"Ready mi_grid with size {self.grid_sizes[i]} and epsilon {self.epsilons[i]}")
        return cgrs_mi_grids

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

    def multifractal_discrimination_analysis(self, cgrs_mi_grids):
        self.total_fractal_points = np.sum(cgrs_mi_grids[0])
        assert np.sum(cgrs_mi_grids[1]) == self.total_fractal_points, "total fractal points should be the same"
        self.fq = []

        for q_index, q in enumerate(self.q_values):
            if q == 1:
                q = self.FIX_1_ERROR_VALUE
            self.fq.append({'q': q, 'fq': np.zeros(len(cgrs_mi_grids))})
            for index, mi_grid in enumerate(cgrs_mi_grids):
                cgr_mi_grid_flattened = mi_grid.reshape(-1)
                no_zeros = cgr_mi_grid_flattened[np.where(cgr_mi_grid_flattened != 0)]
                division = no_zeros / self.total_fractal_points
                powered = np.power(division, q)
                sum_term = np.sum(powered, dtype=np.float64)
                numerator = np.log(sum_term)
                denominator = (q - 1)
                self.fq[-1]['fq'][index] = numerator / denominator

            #logger.info("fq - ", self.fq)
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
        logger.info(f"Finished chromosome: {self.sequence.get_name()}")
        logger.info("*************************************")
        self.result = {
            'q_values': self.q_values,
            'Dq_values': self.Dq_values,
            'tau_q_values': self.tau_q_values,
            'DDq': self.DDq,
            'sequence_name': self.sequence.get_name(),
            'sequence_size': self.sequence.get_size(),
            'refseq_accession_number': self.sequence.get_refseq_accession_number()
        }
        return self.result

    def print_tab_values(self):
        # Print or visualize the calculated Dq and τ(q) values
        print("q\tDq\tτ(q)")
        for q, Dq, tau_q in zip(self.q_values, self.Dq_values, self.tau_q_values):
            print(f"{q}\t{Dq}\t{tau_q}")

    def print_degree_of_multifractality(self):
        # Print the DDq value
        logger.info(f"Degree of Multifractality (DDq): {self.DDq}")

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

    def get_cgr_initial_grid(self):
        return self.cgr_initial_grid

    def set_cgr_initial_grid(self, cgr_initial_grid):
        self.cgr_initial_grid = cgr_initial_grid

    def get_epsilons(self) -> list[float]:
        return self.epsilons

    def get_sizes(self) -> list:
        return self.grid_sizes

    def get_q_min(self) -> int:
        return self.q_min

    def get_q_max(self) -> int:
        return self.q_max

    def get_cgr_gen(self) -> CGR:
        return self.cgr_gen

    def get_fq(self) -> list[dict]:
        return self.fq

    def get_grid_sizes(self) -> list[int]:
        return self.grid_sizes



