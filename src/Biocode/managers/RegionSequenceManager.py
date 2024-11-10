import pandas as pd
import os

from src.Biocode.managers.SequenceManagerInterface import SequenceManagerInterface
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.sequences.RegionSequence import RegionSequence
from src.Biocode.graphs.Graphs import Graphs

from src.Biocode.services.RegionResultsService import RegionResultsService
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.RegionChromosomesService import RegionChromosomesService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.RegionMiGridService import RegionMiGridsService

from utils.logger import logger

from utils.decorators import Inject


@Inject(region_results_service=RegionResultsService,
        organisms_service=OrganismsService,
        region_chromosomes_service=RegionChromosomesService,
        whole_chromosomes_service=WholeChromosomesService,
        region_mi_grids_service=RegionMiGridsService)
class RegionSequenceManager(SequenceManagerInterface):
    def __init__(self, sequence: Sequence = None, sequence_data: dict = None, regions: list[Sequence] = None,
                 sequence_name: str = None, organism_name: str = None, regions_number: int = 0, window_length: int = 0,
                 region_results_service: RegionResultsService = None,
                 organisms_service: OrganismsService = None,
                 region_chromosomes_service: RegionChromosomesService = None,
                 whole_chromosomes_service: WholeChromosomesService = None,
                 region_mi_grids_service: RegionMiGridsService = None):

        self.region_results_service = region_results_service
        self.organisms_service = organisms_service
        self.region_chromosomes_service = region_chromosomes_service
        self.whole_chromosomes_service = whole_chromosomes_service
        self.region_mi_grids_service = region_mi_grids_service
        self.regions_total = regions_number

        if sequence:
            if type(sequence) == RegionSequence:
                self.sequence = sequence
                if sequence_name:
                    self.sequence.set_name(sequence_name)
                if organism_name:
                    self.sequence.set_organism_name(organism_name)
            elif type(sequence) == Sequence:
                self.sequence = RegionSequence(sequence=sequence.get_sequence(), regions_number=regions_number,
                                               name=sequence.get_name(), window_length=window_length,
                                               refseq_accession_number=sequence.get_refseq_accession_number())
            self.sequence_name = self.sequence.get_name()
            self.organism_name = self.sequence.get_organism_name()

        elif sequence_data:
            self.sequence = RegionSequence(sequence_data=sequence_data, regions_number=regions_number,
                                           name=sequence_data['name'], window_length=window_length)
            self.sequence_name = sequence_data['name']
            self.organism_name = sequence_data['organism_name'] or organism_name
        elif regions:
            temp = ''.join(region.get_sequence() for region in regions)
            self.sequence = RegionSequence(sequence=temp, regions_number=regions_number, window_length=window_length,
                                           name=sequence_name)
            self.sequence_name = sequence_name
            self.organism_name = organism_name

        self.regions_total = self.sequence.get_regions_total()
        self.window_length = self.sequence.get_window_length()
        # list of regions (Sequences)
        self.regions = self.sequence.get_regions()
        # managers
        self.managers = []
        self._attach_managers()

        # mfa results
        self.mfa_results = []
        # fq
        self.fq = []
        # degree of multifractality
        self.degrees_of_multifractality = []

        # coverage
        self.cover_percentage = []
        self.cover = []

    def _attach_managers(self):
        for region in self.regions:
            self.managers.append(SequenceManager(sequence=region, sequence_name=region.get_name()))

    def calculate_multifractal_analysis_values(self, GCF):
        self.generate_mfa(GCF, self.region_mi_grids_service, self.region_chromosomes_service)
        self.generate_degrees_of_multifractality()
        self._attach_cover_data()

    def generate_mfa(self, GCF, mi_grids_service, chromosomes_service):
        for manager in self.managers:
            manager.generate_mfa(GCF, mi_grids_service, chromosomes_service)
            self.mfa_results.append(manager.get_mfa_results())
            self.fq.append(manager.get_fq())

    def generate_degrees_of_multifractality(self):
        for manager in self.managers:
            manager.generate_degrees_of_multifractality()
            self.degrees_of_multifractality.append(manager.get_degrees_of_multifractality())

    def graph_cgr(self):
        for manager in self.managers:
            manager.graph_cgr(name=f"regions/{self.organism_name}")

    def graph_3d_cgr(self, grid_size=512):
        for manager in self.managers:
            manager.graph_3d_cgr(grid_size, name=f"regions/{self.organism_name}")

    def graph_linear_fit(self):
        for index, regional_fq in enumerate(self.fq):
            Graphs.graph_linear_fit(fq_values=regional_fq,
                                    epsilons=self.managers[index].get_mfa_generator().get_epsilons(),
                                    sequence_name=self.sequence.get_regions_names()[index],
                                    name=f"regions/{self.organism_name}")

    def graph_degrees_of_multifractality(self, y_range=None, top_labels=False):
        # Check if the lengths of x_array and y_array match
        if len(self.sequence.get_regions_names()) != len(self.degrees_of_multifractality):
            raise ValueError(
                f"Number of chromosome names {len(self.sequence.get_regions_names())} ({self.sequence.get_regions_names()}) \
        does not match the number of degrees of multifractality {len(self.degrees_of_multifractality)} ({self.degrees_of_multifractality}).")

        Graphs.graph_bars(x_array=self.sequence.get_regions_names(), y_array=self.degrees_of_multifractality,
                          title=f"Degree of multifractality by regions for {self.sequence_name}",
                          name=f"regions/{self.organism_name}",
                          y_label="Degree of multifractality", y_range=y_range, top_labels=top_labels)

    def graph_multifractal_spectrum(self, color_by='region'):
        Graphs.graph_many_grouped(results_array=self.flattened_mfa_results, X='q_values', Y='Dq_values', x_label='q',
                                  y_label='Dq',
                                  title=f"Multifractal spectrum of chromosome {self.sequence_name} by regions of {self.organism_name}",
                                  name=f"regions/{self.organism_name}",
                                  regions_number=self.sequence.get_regions_number(), labels_array=self.sequence.get_regions_names(),
                                  color_by=color_by)
    def graph_multifractal_spectrum_old(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='Dq_values', x_label='q', y_label='Dq',
                          title=f'Dq vs q by regions for {self.sequence_name}', name=f"regions/{self.organism_name}",
                          labels_array=self.sequence.get_regions_names())

    def graph_correlation_exponent(self, color_by='region', markersize=3):
        Graphs.graph_many_grouped(results_array=self.flattened_mfa_results, X='q_values', Y='tau_q_values', x_label='q',
                                  y_label='t(q)',
                                  title=f"Correlation exponent of chromosome {self.sequence_name} by regions of {self.organism_name}",
                                  name=f"regions/{self.organism_name}",
                                  regions_number=self.sequence.get_regions_number(), labels_array=self.sequence.get_regions_names(),
                                  markersize=markersize, color_by=color_by)

    def graph_correlation_exponent_old(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='tau_q_values', x_label='q', y_label='t(q)',
                          title=f't(q) vs q by regions for {self.sequence_name}', name=f"regions/{self.organism_name}",
                          labels_array=self.sequence.get_regions_names(), markersize=2)

    def _attach_cover_data(self):
        for sequence in self.regions:
            self.cover_percentage.append(sequence.get_cover_percentage())
            self.cover.append(sequence.get_cover())

    def graph_multifractal_analysis(self, y_range_degrees_of_multifractality=None, _3d_cgr=False, linear_fit=True,
                                    degrees_of_multifractality=True,
                                    multifractal_spectrum=True, correlation_exponent=True, top_labels=True):
        super().graph_multifractal_analysis(_3d_cgr, linear_fit, degrees_of_multifractality, multifractal_spectrum,
                                            correlation_exponent)
        if degrees_of_multifractality:
            self.graph_degrees_of_multifractality(y_range=y_range_degrees_of_multifractality, top_labels=top_labels)

    def graph_multifractal_analysis_merged(self,  multifractal_spectrum=True, correlation_exponent=True,
                                           regions_labels=None,
                                           color_by='region'):
        if multifractal_spectrum:
            self.graph_multifractal_spectrum(color_by=color_by)
        if correlation_exponent:
            self.graph_correlation_exponent(color_by=color_by, markersize=3)

    def graph_coverage(self):
        for manager in self.managers:
            manager.set_organism_name(self.organism_name)
            manager.graph_coverage(subfolder="regions")

    def generate_df_results(self, selected_columns=None):
        sheet = "With regions"
        # Extract sequence names and use them as row labels
        row_labels = self.sequence.get_regions_names()
        # Extract sequence names and use them as row labels
        q_min = self.get_managers()[0].get_mfa_generator().get_q_min()
        q_max = self.get_managers()[0].get_mfa_generator().get_q_max()
        # Create an empty DataFrame with the correct number of rows
        self.df_results = pd.DataFrame(index=row_labels, columns=["D%d" % i for i in range(q_min, q_max + 1)])

        # Code to populate the DataFrame
        for i, data_dict in enumerate(self.flattened_mfa_results):
            self.df_results.loc[row_labels[i]] = data_dict['Dq_values']

        self.df_results['DDq'] = [data_dict['DDq'] for data_dict in self.flattened_mfa_results]
        self.df_results['t(q=20)'] = [data_dict['tau_q_values'][-1] for data_dict in self.flattened_mfa_results]

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

    def set_sequence_name(self, sequence_name):
        self.sequence_name = sequence_name

    def set_mfa_results(self, mfa_results: list):
        self.mfa_results = mfa_results

    def set_degrees_of_multifractality(self, degrees_of_multifractality: list):
        self.degrees_of_multifractality = degrees_of_multifractality

    def set_flattened_mfa_results(self, flattened_mfa_results: list):
        self.flattened_mfa_results = flattened_mfa_results

    def set_cover(self, cover):
        for index, manager in enumerate(self.managers):
            manager.set_cover(cover[index])

    def set_cover_percentage(self, cover_percentage):
        for index, manager in enumerate(self.managers):
            manager.set_cover_percentage(cover_percentage[index])

    def get_sequence_name(self):
        return self.sequence_name

    def get_sequence(self):
        return self.sequence

    def get_mfa_results(self):
        return self.mfa_results

    def get_regions_names(self):
        return self.sequence.get_regions_names()

    def get_regions(self):
        return self.regions

    def get_degrees_of_multifractality(self):
        return self.degrees_of_multifractality

    def get_cover(self) -> list[list[int]]:
        return self.cover

    def get_cover_percentage(self) -> list[float]:
        return self.cover_percentage

    def get_managers(self) -> list[SequenceManager]:
        return self.managers

    def get_whole_refseq_accession_number(self) -> str:
        return self.sequence.get_refseq_accession_number()

    def get_regions_refseq_accessions_numbers(self) -> list[str]:
        return [region.get_refseq_accession_number() for region in self.regions]

    def save_results_to_db_during_execution(self, GCF):
        """
        [(val1, val2), (val1, val2)]
        ["chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        [{"q_values", "Dq_values", "tau_q_values", "DDq"}]
        """
        #try:
        self.organism_id = self.organisms_service.extract_by_GCF(GCF=GCF)
        whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
            self.sequence.get_refseq_accession_number())
        for index, result in enumerate(self.mfa_results):
            chromosome_id = self.region_chromosomes_service.extract_id_by_refseq_accession_number(
                self.regions[index].get_refseq_accession_number())
            self.region_chromosomes_service.update_when_null(pk_value=chromosome_id,
                                                             record=(self.regions[index].get_name(),
                                                                     self.regions[index].get_refseq_accession_number(),
                                                                     self.organism_id,
                                                                     self.cover_percentage[index],
                                                                     self.cover[index],
                                                                     self.sequence.get_regions_number(),
                                                                     index + 1,
                                                                     result['sequence_size'],
                                                                     whole_chromosome_id))
            self.region_results_service.insert(
                record=(chromosome_id, result['Dq_values'].tolist(),
                        result['tau_q_values'].tolist(),
                        result['DDq']))
        del self.mfa_results
        # except AttributeError as e:
        #    logger.error(f"Was not able to extract ids from organism and whole_chromosome to insert "
        #                 f"a new region sequence {self.sequence_name} with error name: {e.name}")
        logger.info(f"************* Saved to DB {self.sequence_name} *************")
