from src.Biocode.managers.SequenceManagerInterface import SequenceManagerInterface
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.mfa.MFA import MFA
from src.Biocode.graphs.Graphs import Graphs
from src.Biocode.recursive_sequences_finder.RecursiveSequencesFinder import RecursiveSequencesFinder
from src.Biocode.dataclasses.MiGridCoordinatesValuesAndNucleotides import MiGridCoordinatesValuesAndNucleotides

from src.Biocode.services.WholeResultsService import WholeResultsService
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.RecursiveRepeatsWholeChromosomesService import RecursiveRepeatsWholeChromosomesService
from src.Biocode.services.WholeMiGridsService import WholeMiGridsService

from utils.decorators import Inject

from utils.logger import logger
from typing import List
import gc
import numpy as np

@Inject(whole_results_service=WholeResultsService,
        organisms_service = OrganismsService,
        whole_chromosomes_service = WholeChromosomesService,
        repeats_service = RepeatsService,
        recursive_repeats_whole_chromosomes_service = RecursiveRepeatsWholeChromosomesService,
        whole_mi_grids_service = WholeMiGridsService
        )
class SequenceManager(SequenceManagerInterface):
    def __init__(self, sequence: Sequence = None, sequence_data: dict = None, sequence_name: str = None,
                 organism_name: str = None,
                 whole_results_service: WholeResultsService = None,
                 organisms_service: OrganismsService = None,
                 whole_chromosomes_service: WholeChromosomesService = None,
                 repeats_service: RepeatsService = None,
                 recursive_repeats_whole_chromosomes_service: RecursiveRepeatsWholeChromosomesService = None,
                 whole_mi_grids_service: WholeMiGridsService = None):

        self.chromosome_id = None
        self.organism_id = None
        self.whole_chromosomes_service = whole_chromosomes_service
        self.whole_results_service = whole_results_service
        self.organisms_service = organisms_service
        self.repeats_service = repeats_service
        self.recursive_repeats_whole_chromosomes_service = recursive_repeats_whole_chromosomes_service
        self.whole_mi_grids_service = whole_mi_grids_service


        self.n_largest_mi_grid_values_strings_for_k = None
        if sequence:
            self.sequence = sequence
            self.refseq_accession_number = sequence.get_refseq_accession_number()
            if sequence_name:
                self.sequence.set_name(sequence_name)
            if organism_name:
                self.sequence.set_organism_name(organism_name)

            self.sequence_name = self.sequence.get_name()
            self.organism_name = self.sequence.get_organism_name()
        elif sequence_data:
            self.sequence = Sequence(sequence_data=sequence_data)
            self.sequence_name = sequence_data['name']
            self.organism_name = sequence_data['organism_name'] or organism_name
            self.refseq_accession_number = self.sequence.get_refseq_accession_number()

        self.mfa_generator = MFA(sequence=self.sequence)
        self.mfa_results = None
        self.fq = None
        self.degree_of_multifractality = None
        self.cover = None
        self.cover_percentage = None
        self._10_largest_values_from_k_10_to_4 = None



    def generate_mfa(self, GCF):
        mi_grid_whole_chromosome_id = self.whole_mi_grids_service\
            .get_whole_chromosome_id_if_mi_grid_exists_by_refseq_accession_number(self.refseq_accession_number)
        if mi_grid_whole_chromosome_id is None:
            chromosome_id = self._insert_chromosome(GCF)
            logger.info(f"No previous results found in mi_grids table for {self.refseq_accession_number} - executing CGR algorithm")
            initial_cgr = self.mfa_generator.generate_initial_grid()
            self._save_initial_grid_to_database(initial_cgr, chromosome_id)
            cgr_results = self.mfa_generator.generate_cgr_mi_grids_from_initial_grid(initial_cgr)
        else:
            logger.info(f"Extracting {self.refseq_accession_number} mi_grid from DB")
            retrieved_data = self.whole_mi_grids_service.extract_mi_grid_by_chromosome_id(mi_grid_whole_chromosome_id)
            initial_cgr = np.frombuffer(retrieved_data, dtype=np.float64)
            initial_cgr = initial_cgr.reshape((self.mfa_generator.get_grid_sizes()[0],
                                               self.mfa_generator.get_grid_sizes()[0]))
            cgr_results = self.mfa_generator.generate_cgr_mi_grids_from_initial_grid(initial_cgr)

        self.mfa_results = self.mfa_generator.multifractal_discrimination_analysis(cgrs_mi_grids=cgr_results)
        logger.info(f"MFA results were generated for chromosome: {self.refseq_accession_number}")
        self.fq = self.mfa_generator.get_fq()

    def _save_initial_grid_to_database(self, cgr_largest_mi_grid, whole_chromosome_id):
        cgr_largest_mi_grid_np = np.array(cgr_largest_mi_grid, dtype=np.float64)

        binary_data_numpy = cgr_largest_mi_grid_np.tobytes()

        del cgr_largest_mi_grid_np
        self.whole_mi_grids_service.insert(record=(binary_data_numpy, whole_chromosome_id,
                                                   self.mfa_generator.get_epsilons()[0]))
        del binary_data_numpy
        gc.collect()


    def generate_degree_of_multifractality(self):
        self.degree_of_multifractality = self.mfa_generator.get_DDq()

    def graph_cgr(self):
        self.mfa_generator.get_cgr_gen().graph_cgr(title=f"CGR of {self.sequence_name}", name=f"{self.organism_name}/whole")

    def graph_3d_cgr(self, grid_size=512):
        epsilons = self.mfa_generator.get_epsilons()
        sizes = self.mfa_generator.get_sizes()
        for index, mi_grid in enumerate(self.mfa_generator.get_cgrs_mi_grids()):
            if sizes[index] == grid_size:
                Graphs.graph_3d_cgr(count_matrix=mi_grid, name=f"{self.organism_name}/whole",
                                    title=f"MMR for {self.sequence.get_name()} with grid size {sizes[index]}",
                                    epsilon=epsilons[index])

    def graph_linear_fit(self):
        Graphs.graph_linear_fit(fq_values=self.fq, epsilons=self.mfa_generator.get_epsilons(),
                                sequence_name=self.sequence_name, name=f"{self.organism_name}/whole")

    def graph_multifractal_spectrum(self):
        Graphs.graph_many(results_array=[self.mfa_results], X='q_values', Y='Dq_values', x_label='q', y_label='Dq',
                          title=f'Dq vs q for {self.sequence_name}', name=f"{self.organism_name}/whole",
                          labels_array=[self.sequence.get_name()])

    def graph_correlation_exponent(self):
        Graphs.graph_many(results_array=[self.mfa_results], X='q_values', Y='tau_q_values', x_label='q', y_label='t(q)',
                          title=f't(q) vs q for chromosomes {self.sequence_name}',
                          name=f"{self.organism_name}/whole",
                          labels_array=[self.sequence.get_name()], markersize=3)

    def _attach_cover_data(self):
        self.cover_percentage = self.sequence.get_cover_percentage()
        self.cover = self.sequence.get_cover()

    def graph_multifractal_analysis(self, _3d_cgr=False, linear_fit=True, degrees_of_multifractality=True,
                                    multifractal_spectrum=True, correlation_exponent=True, top_labels=True):
        super().graph_multifractal_analysis(_3d_cgr, linear_fit, degrees_of_multifractality,
                                            multifractal_spectrum, correlation_exponent)
        if degrees_of_multifractality:
            print(f"The degree of multifractality of {self.sequence_name} is {self.degree_of_multifractality}")

    def graph_coverage(self, subfolder="whole"):
        Graphs.graph_coverage(values=self.cover, sequence_name=self.sequence_name, name=f"{self.organism_name}/{subfolder}")


    def find_nucleotides_strings_recursively(self, k1: int, k2: int, k_step:int, amount_sequences: int) -> List[MiGridCoordinatesValuesAndNucleotides]:
        recursive_sequences_finder = RecursiveSequencesFinder(grid_exponents=self.mfa_generator.GRID_EXPONENTS,
                                                              cgrs_mi_grids=self.mfa_generator.get_cgrs_mi_grids())
        self.n_largest_mi_grid_values_strings_for_k = recursive_sequences_finder.get_largest_mi_grid_values_strings_for_different_k(
            k1=k1, k2=k2, k_step=k_step, amount_sequences=amount_sequences
        )
        return self.n_largest_mi_grid_values_strings_for_k

    def set_organism_name(self, organism_name):
        self.organism_name = organism_name

    def set_cover(self, cover):
        self.cover = cover

    def set_cover_percentage(self, cover_percentage):
        self.cover_percentage = cover_percentage

    def get_sequence_name(self) -> str:
        return self.sequence_name

    def get_mfa_generator(self) -> MFA:
        return self.mfa_generator

    def get_degree_of_multifractality(self) -> float:
        return self.degree_of_multifractality

    def get_mfa_results(self) -> dict:
        return self.mfa_results

    def get_fq(self) -> list[dict]:
        return self.fq

    def get_cover(self) -> list[int]:
        return self.cover

    def get_cover_percentage(self) -> float:
        return self.cover_percentage

    def get_organism_id(self) -> int:
        return self.organism_id

    def get_chromosome_id(self) -> int:
        return self.chromosome_id

    def get_refseq_accession_number(self) -> str:
        return self.refseq_accession_number


    def save_results_to_db_during_execution(self, GCF):
        self._insert_chromosome(GCF)
        self.whole_results_service.insert(record=(self.chromosome_id, self.mfa_results['Dq_values'].tolist(),
                                             self.mfa_results['tau_q_values'].tolist(),
                                             self.mfa_results['DDq']))
        del self.mfa_results
        logger.info(f"************* Saved to DB {self.sequence_name} *************")


    def find_only_kmers_recursively(self) -> List[MiGridCoordinatesValuesAndNucleotides]:
        kmers_list = self.find_nucleotides_strings_recursively(k1=10, k2=4, k_step=-1, amount_sequences=10)
        logger.info("kmers - " + str(kmers_list))
        return kmers_list


    def save_repeats_found_recursively_to_db(self, kmers_list: List[MiGridCoordinatesValuesAndNucleotides], GCF: str,
                                             method_to_find_it: str = "Recursively"):

        self.organism_id = self.organisms_service.extract_by_GCF(GCF=GCF)

        self.chromosome_id = self._insert_chromosome()

        for kmers in kmers_list:
            nucleotides_strings = kmers.get_nucleotides_strings()
            largest_values = kmers.get_largest_values()
            coordinates = kmers.get_coordinates()
            for index, string in enumerate(nucleotides_strings):
                repeats_service_id = self.repeats_service.insert(record=(string, "", method_to_find_it))
                self.recursive_repeats_whole_chromosomes_service.insert(
                    record=(repeats_service_id, self.chromosome_id , len(string), int(largest_values[index]),
                            str(coordinates[index])))


    def _insert_chromosome(self, GCF) -> int:
        self.organism_id = self.organisms_service.extract_by_GCF(GCF=GCF)
        whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(self.refseq_accession_number)
        if whole_chromosome_id is None:
            return self.whole_chromosomes_service.insert(
                record=(self.sequence.get_name(), self.refseq_accession_number, self.organism_id,
                        self.cover_percentage,
                        self.cover,
                        self.sequence.get_size()))
        else:
            return self.whole_chromosomes_service.update(
                pk_value=whole_chromosome_id,
                record=(self.sequence.get_name(), self.refseq_accession_number, self.organism_id,
                        self.cover_percentage,
                        self.cover,
                        self.sequence.get_size()))


