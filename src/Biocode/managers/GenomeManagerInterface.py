import pandas as pd
import os

from src.Biocode.sequences.Genome import Genome
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager

from src.Biocode.dataclasses.MiGridCoordinatesValuesAndNucleotides import MiGridCoordinatesValuesAndNucleotides


from utils.logger import logger
from typing import List
import gc


class GenomeManagerInterface:

    def __init__(self, genome: Genome = None, genome_data: list[dict] = None, chromosomes: list[Sequence] = None,
                 organism_name: str = None, window_length: int = 0, regions_number: int = 0):
        self.n_largest_mi_grid_values_strings = None
        self.df_results = None

        self.genome = None
        if genome:
            self.genome = genome
        elif genome_data:
            self.genome = Genome(chromosomes_data=genome_data, regions_number=regions_number, window_length=window_length)
        elif chromosomes:
            self.genome = Genome(chromosomes=chromosomes, regions_number=regions_number, window_length=window_length)

        # Managers
        if not self.genome:
            logger.warning(f"There is not self.genome object in GenomeManagerInterface")
        else:
            self.managers = []
            if (regions_number is None or regions_number < 0) and (
                    window_length is None or window_length < 0):
                raise Exception('Enter a valid regions_1 number or window length for the GenomeManagerInterface')
            elif regions_number == 0 and window_length == 0:
                for chromosome in self.genome.get_chromosomes():
                    self.managers.append(SequenceManager(sequence=chromosome, sequence_name=chromosome.get_name(),
                                                         organism_name=organism_name))
            else:  # > 0
                for chromosome in self.genome.get_chromosomes():
                    self.managers.append(RegionSequenceManager(sequence=chromosome, sequence_name=chromosome.get_name(),
                                                               regions_number=regions_number, window_length=window_length,
                                                               organism_name=organism_name))
                # regions_1 names
                self.regions_names = []
                self._attach_regions_names()

        #self.regions_number = regions_number or self.genome.get_chromosomes()[0].get_regions_number()
        self.window_length = window_length

        # name of organism
        self.organism_name = organism_name
        # mfa resultsc
        self.mfa_results = []
        # degrees of multifractality
        self.degrees_of_multifractality = []
        # coverage
        self.cover_percentage = []
        self.cover = []

    def graph_cgr(self):
        """Graph CGR for all the chromosomes"""
        for manager in self.managers:
            manager.graph_cgr()

    def graph_3d_cgr(self, all=True, index=0, grid_size=512):
        """Graph CGR density in a 3D chart"""
        if not all:
            self.managers[index].graph_3d_cgr(grid_size)
        else:
            for manager in self.managers:
                manager.graph_3d_cgr(grid_size)

    def graph_linear_fit(self):
        """Graph linear fit for fq vs ln(epsilon)"""
        for manager in self.managers:
            manager.graph_linear_fit()

    def generate_degrees_of_multifractality(self):
        """Generate degrees of multifractality"""
        pass

    def graph_multifractal_spectrum(self):
        """Graph multifractal spectrum of Dq vs q"""
        pass

    def graph_correlation_exponent(self):
        """Graph t(q) vs q"""
        pass

    def calculate_multifractal_analysis_values(self, GCF: str, save_to_db: bool, name: str, only_cgr: bool):
        """Generate mfa generators, generate mfa values, the cover and cover
        percentage and send them to the DB"""
        for manager in self.managers:
            logger.info(f"Starting chromosome: {manager.get_sequence_name()}")
            manager.calculate_multifractal_analysis_values(GCF, name, only_cgr)
            if save_to_db:
                manager.save_results_to_db_during_execution(GCF=GCF)
                logger.info(f"Saved results for GCF: {GCF}")
            del manager
            gc.collect()

        """
        self.n_largest_mi_grid_values_strings = self._find_nucleotides_strings_recursively(
            manager=manager, k1=10, k2=4, k_step=-1, amount_sequences=10
        )
        logger.critical(self.n_largest_mi_grid_values_strings)
        
        """


    def find_only_kmers_recursively_and_calculate_multifractal_analysis_values(self, GCF: str, save_to_db: bool,
                                                                               method_to_find_it: str, k_range: tuple):
        for manager in self.managers:
            manager.calculate_multifractal_analysis_values(GCF)
            kmers_list = self._find_nucleotides_strings_recursively(manager, k1=k_range[1], k2=k_range[0], k_step=-1,
                                                                    amount_sequences=10)
            if save_to_db:
                manager.save_repeats_found_recursively_to_db(kmers_list=kmers_list, GCF=GCF,
                                                             method_to_find_it=method_to_find_it)
            del manager


    def _find_nucleotides_strings_recursively(self, manager: SequenceManager, k1: int, k2: int, k_step: int,
                                              amount_sequences: int) -> List[MiGridCoordinatesValuesAndNucleotides]:
        return manager.find_nucleotides_strings_recursively(k1, k2, k_step, amount_sequences)

    def graph_multifractal_analysis(self, _3d_cgr=True, degrees_of_multifractality=True,
                                    multifractal_spectrum=True, correlation_exponent=True, top_labels=True):
        """Graph the multifractal analysis values using the chart of 3D density of points, the linear fit fq vs q, 
        the multifractal spectrum Dq vs q, and the correlation exponent t(q) vs q"""
        for manager in self.managers:
            manager.graph_multifractal_analysis(_3d_cgr, degrees_of_multifractality, multifractal_spectrum,
                                                correlation_exponent, top_labels=top_labels)

    def graph_multifractal_analysis_merged(self):
        """Graph mutifractal analysis charts with merged values"""
        pass

    def graph_coverage(self):
        """Graph the barplot representing the coverage of the DNA sequence (the nucleotides that have been 
        identified)"""
        for manager in self.managers:
            manager.graph_coverage()

    def calculate_and_graph(self):
        """Generate MFA values and graph them along with the coverage"""
        pass

    def generate_df_results(self, mfa_results, row_labels, q_min, q_max, sheet, selected_columns: list = None):
        # Create an empty DataFrame with the correct number of rows
        self.df_results = pd.DataFrame(index=row_labels, columns=["D%d" % i for i in range(q_min, q_max + 1)])

        # Code to populate the DataFrame
        for i, data_dict in enumerate(mfa_results):
            self.df_results.loc[row_labels[i]] = data_dict['Dq_values']

        self.df_results['DDq'] = [data_dict['DDq'] for data_dict in mfa_results]
        self.df_results['t(q=20)'] = [data_dict['tau_q_values'][-1] for data_dict in mfa_results]

        # Set the index based on the type of row_labels
        if any('_region_' in label for label in row_labels):
            # Row labels correspond to chromosome regions_1
            self.df_results.index.name = 'Region'
        else:
            # Row labels correspond to whole chromosomes
            self.df_results.index.name = 'Chromosome'

        selected_df_results = self.df_results[selected_columns] if selected_columns else self.df_results
        self._save_df_results(selected_df_results, sheet)
        return selected_df_results

    def _save_df_results(self, df, sheet):
        directory = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
                                 "out/results")
        os.makedirs(directory, exist_ok=True)

        # File path for the Excel file
        output_file = f'{directory}/{self.organism_name}.xlsx'

        # Check if the file already exists
        if os.path.exists(output_file):
            # Load the existing Excel file
            with pd.ExcelFile(output_file) as xls:
                # Read the existing sheets into a dictionary of DataFrames
                sheets = {sheet_name: xls.parse(sheet_name, index_col=0) for sheet_name in xls.sheet_names}

            # Add the new sheet to the dictionary
            sheets[sheet] = df
            # Save all sheets back to the Excel file
            with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
                for sheet_name, sheet_df in sheets.items():
                    sheet_df.to_excel(writer, sheet_name=sheet_name, index=True)

        else:
            # If the file does not exist, create a new one
            with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
                df.to_excel(writer, sheet_name=sheet, index=True)

    def save_to_db_after_execution(self, GCF):
        pass

    def set_organism_name(self, organism_name):
        self.organism_name = organism_name

    def set_mfa_results(self, mfa_results: list):
        self.mfa_results = mfa_results

    def set_flattened_mfa_results(self, flattened_mfa_results: list):
        self.flattened_mfa_results = flattened_mfa_results

    def set_degrees_of_multifractality(self, degrees_of_multifractality: list):
        self.degrees_of_multifractality = degrees_of_multifractality

    def set_cover(self, cover: list):
        pass

    def set_cover_percentage(self, cover_percentage: list):
        pass