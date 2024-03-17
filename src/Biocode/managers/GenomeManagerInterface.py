import pandas as pd
import os
import sqlite3

from src.Biocode.sequences.Genome import Genome
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager


from src.Biocode.utils.utils import list_to_str, str_to_list
from utils.logger import logger


class GenomeManagerInterface:

    def __init__(self, genome: Genome = None, genome_data: list[dict] = None, chromosomes: list[Sequence] = None,
                 organism_name: str = None,
                 regions_number: int = 0):
        self.df_results = None
        self.regions_number = regions_number
        if genome:
            self.genome = genome
        elif genome_data:
            self.genome = Genome(chromosomes_data=genome_data, regions_number=regions_number)
        elif chromosomes:
            self.genome = Genome(chromosomes=chromosomes, regions_number=regions_number)

        # Managers
        if self.genome:
            self.managers = []
            if regions_number < 0:
                raise Exception("Not a valid regions_number for the GenomeManager constructor")
            elif regions_number == 0:
                for chromosome in self.genome.get_chromosomes():
                    self.managers.append(SequenceManager(sequence=chromosome, sequence_name=chromosome.get_name(),
                                                         organism_name=organism_name))
            else:  # > 0
                for chromosome in self.genome.get_chromosomes():
                    self.managers.append(RegionSequenceManager(sequence=chromosome, sequence_name=chromosome.get_name(),
                                                               regions_number=regions_number,
                                                               organism_name=organism_name))
                # regions names
                self.regions_names = []
                self._attach_regions_names()


        # name of organism
        self.organism_name = organism_name
        # mfa results
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

    def graph_linear_fit(self, name):
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

    def calculate_multifractal_analysis_values(self):
        """Generate mfa generators, generate mfa values, attach the degrees of multifractality, the cover and cover 
        percentage"""
        for manager in self.managers:
            logger.info(f"Starting chromosome: {manager.get_sequence_name()}")
            manager.calculate_multifractal_analysis_values()
            self.cover.append(manager.get_cover())
            self.cover_percentage.append(manager.get_cover_percentage())
            self.mfa_results.append(manager.get_mfa_results())

        self.generate_degrees_of_multifractality()

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
            # Row labels correspond to chromosome regions
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

    def save_to_db(self, GCF):
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