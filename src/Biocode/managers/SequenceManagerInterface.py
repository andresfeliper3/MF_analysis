from typing import List
from src.Biocode.dataclasses.MiGridCoordinatesValuesAndNucleotides import MiGridCoordinatesValuesAndNucleotides
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.RegionChromosomesService import RegionChromosomesService

from utils.logger import logger
from src.Biocode.utils.utils import remove_region_part


class SequenceManagerInterface:
    def generate_mfa_generator(self):
        """Generate the MFA class"""
        pass

    def generate_mfa(self, GCF, mi_grids_service, chromosomes_service, name, only_cgr):
        """Generate mfa_results with the MFA values"""
        pass

    def generate_degrees_of_multifractality(self):
        """Generate degrees of multifractality"""
        pass

    def graph_cgr(self):
        "Generate CGR without the grid"
        pass

    def graph_3d_cgr(self, grid_size):
        """Graph CGR density in a 3D chart"""
        pass

    def graph_linear_fit(self):
        """Graph linear fit for fq vs ln(epsilon)"""
        pass

    def graph_multifractal_spectrum(self):
        """Graph multifractal spectrum of Dq vs q"""
        pass

    def graph_correlation_exponent(self):
        """Graph t(q) vs q"""
        pass

    def _attach_cover_data(self):
        """Attach to the fields of the class the cover and cover_percentage of the sequence"""
        pass

    def calculate_multifractal_analysis_values(self, GCF, name, only_cgr):
        """Generate mfa generators, generate mfa values, attach the degrees of multifractality, the cover and cover
        percentage"""
        pass

    def graph_multifractal_analysis(self, _3d_cgr=False, linear_fit=True, degrees_of_multifractality=None,
                                    multifractal_spectrum=True, correlation_exponent=True):
        """Graph the multifractal analysis values using the chart of 3D density of points, the linear fit fq vs q,
        the multifractal spectrum Dq vs q, and the correlation exponent t(q) vs q"""
        if _3d_cgr:
            self.graph_3d_cgr()
        if linear_fit:
            self.graph_linear_fit()
        if multifractal_spectrum:
            self.graph_multifractal_spectrum()
        if correlation_exponent:
            self.graph_correlation_exponent()

    def graph_coverage(self):
        """Graph the barplot representing the coverage of the DNA sequence (the nucleotides that have been
        identified)"""
        pass

    def calculate_and_graph(self, GCF:str):
        """Generate MFA values and graph them along with the coverage"""
        self.calculate_multifractal_analysis_values(GCF)
        self.graph_multifractal_analysis()
        self.graph_coverage()

    def find_nucleotides_strings_recursively(self, k1: int, k2: int, k_step:int, amount_sequences: int) -> List[MiGridCoordinatesValuesAndNucleotides]:
        """Find nucleotides strings recursively using CGR"""
        pass

    def set_cover(self, cover):
        pass

    def set_cover_percentage(self, cover_percentage):
        pass

    def save_results_to_db_during_execution(self, GCF: str):
       """ Save to DB after each sequence is executed"""
       pass

    def _insert_chromosome(self, GCF, chromosomes_service) -> int:
        self.organism_id = self.organisms_service.extract_by_GCF(GCF=GCF)
        chromosome_id = chromosomes_service.extract_id_by_refseq_accession_number(self.refseq_accession_number)

        if isinstance(chromosomes_service, WholeChromosomesService):
            record = (self.sequence.get_name(), self.refseq_accession_number, self.organism_id, self.cover_percentage,
                        self.cover, self.sequence.get_size())
            if chromosome_id is None:
                return chromosomes_service.insert(record=record)
            else:
                result = chromosomes_service.update_when_null(pk_value=chromosome_id, record=record)
                return result

        elif isinstance(chromosomes_service, RegionChromosomesService):
            whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
                remove_region_part(self.sequence.get_refseq_accession_number()))

            record = (self.sequence_name, self.refseq_accession_number, self.organism_id,
                            self.cover_percentage, self.cover, self.sequence.get_regions_total(),
                            self.sequence.get_region_number(), self.sequence.get_size(), whole_chromosome_id)
            if chromosome_id is None:
                return chromosomes_service.insert(record=record)
            else:
                return chromosomes_service.update_when_null(pk_value=chromosome_id, record=record)