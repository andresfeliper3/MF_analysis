from src.Biocode.managers.GenomeManagerInterface import GenomeManagerInterface
from src.Biocode.sequences.Genome import Genome
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.graphs.Graphs import Graphs

from src.Biocode.services.WholeResultsService import WholeResultsService
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.ChromosomesService import ChromosomesService

from src.Biocode.utils.utils import list_to_str, str_to_list

import os


class GenomeManager(GenomeManagerInterface):
    def __init__(self, genome: Genome = None, genome_data: list[dict] = None, chromosomes: list[Sequence] = None,
                 organism_name: str = None):
        super().__init__(genome, genome_data, chromosomes, organism_name, 0)

    def generate_degrees_of_multifractality(self):
        for manager in self.managers:
            manager.generate_degree_of_multifractality()
            self.degrees_of_multifractality.append(manager.get_degree_of_multifractality())

    def graph_degrees_of_multifractality(self, y_range=None, top_labels=False):
        # Check if the lengths of x_array and y_array match
        y_array = [mfa_result['DDq'] for mfa_result in self.mfa_results] if len(
            self.degrees_of_multifractality) == 0 else self.degrees_of_multifractality

        if len(self.genome.get_chromosomes_names()) != len(y_array):
            raise ValueError(
                f"Number of chromosome names {len(self.genome.get_chromosomes_names())} ({self.genome.get_chromosomes_names()}) \
        does not match the number of y_array {len(y_array)} ({y_array}).")

        Graphs.graph_bars(x_array=self.genome.get_chromosomes_names(), y_array=y_array,
                          title=f"Degree of multifractality by chromosomes of {self.organism_name}",
                          name=f"{self.organism_name}/whole",
                          y_label="Degree of multifractality", y_range=y_range, top_labels=top_labels)

    def graph_multifractal_spectrum(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='Dq_values', x_label='q', y_label='Dq',
                          title=f"Dq vs q by chromosomes of {self.organism_name}", name=f"{self.organism_name}/whole",
                          labels_array=self.genome.get_chromosomes_names())

    def graph_correlation_exponent(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='tau_q_values', x_label='q', y_label='t(q)',
                          title=f"t(q) vs q by chromosomes of {self.organism_name}", name=f"{self.organism_name}/whole",
                          labels_array=self.genome.get_chromosomes_names(), markersize=3)

    def graph_multifractal_analysis_merged(self, y_range_degrees_of_multifractality=None,
                                           degrees_of_multifractality=True,
                                           multifractal_spectrum=True, correlation_exponent=True, top_labels=True):
        if multifractal_spectrum:
            self.graph_multifractal_spectrum()
        if correlation_exponent:
            self.graph_correlation_exponent()
        if degrees_of_multifractality:
            self.graph_degrees_of_multifractality(y_range=y_range_degrees_of_multifractality, top_labels=top_labels)

    def calculate_and_graph(self):
        self.calculate_multifractal_analysis_values()
        self.graph_multifractal_analysis()

    def calculate_and_graph_plus_merged(self):
        self.calculate_multifractal_analysis_values()
        self.graph_multifractal_analysis()
        self.graph_multifractal_analysis_merged()

    def calculate_and_graph_only_merged(self):
        self.calculate_multifractal_analysis_values()
        self.graph_multifractal_analysis_merged()

    def generate_df_results(self, selected_columns: list = None):
        # Extract sequence names and use them as row labels
        row_labels = self.genome.get_chromosomes_names()
        # Extract sequence names and use them as row labels
        q_min = self.managers[0].get_mfa_generator().get_q_min()
        q_max = self.managers[0].get_mfa_generator().get_q_max()
        return super().generate_df_results(self.mfa_results, row_labels, q_min, q_max, "Whole Genome",
                                           selected_columns)

    def save_to_db(self, GCF):
        whole_results_service = WholeResultsService()
        organisms_service = OrganismsService()
        chromosomes_service = ChromosomesService()
        """
        [(val1, val2), (val1, val2)]
        ["chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        [{"q_values", "Dq_values", "tau_q_values", "DDq"}]
        """
        organism_id = int(organisms_service.extract_by_GCF(GCF=GCF).loc[0, 'id'])

        for index, result in enumerate(self.mfa_results):
            chromosome_id = chromosomes_service.insert(record=(result['sequence_name'], organism_id,
                                                               self.cover_percentage[index],
                                                               list_to_str(self.cover[index])))
            whole_results_service.insert(record=(chromosome_id, list_to_str(result['Dq_values'].tolist()),
                                                 list_to_str(result['tau_q_values'].tolist()),
                                                 list_to_str(result['DDq'])))

    def set_cover(self, cover: list):
        self.cover = cover
        for index, manager in enumerate(self.managers):
            manager.set_cover(self.cover[index])

    def set_cover_percentage(self, cover_percentage: list):
        self.cover_percentage = cover_percentage
        for index, manager in enumerate(self.managers):
            manager.set_cover_percentage(self.cover_percentage[index])

    def get_organism_name(self):
        return self.organism_name

    def get_cgr_generators(self):
        return self.cgr_generators

    def get_genome(self):
        return self.genome

    def get_mfa_results(self):
        return self.mfa_results

    def get_degrees_of_multifractality(self):
        return self.degrees_of_multifractality

    def get_managers(self) -> list:
        return self.managers

    def get_df_results(self):
        return self.df_results

    def get_cover(self) -> list[list[int]]:
        return self.cover

    def get_cover_percentage(self) -> list[float]:
        return self.cover_percentage
