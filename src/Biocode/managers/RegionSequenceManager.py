from src.Biocode.managers.SequenceManagerInterface import SequenceManagerInterface
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.sequences.RegionSequence import RegionSequence
from src.Biocode.graphs.Graphs import Graphs

from src.Biocode.services.RegionResultsService import RegionResultsService
from src.Biocode.services.OrganismsService import OrganismsService
from src. Biocode.services.RegionChromosomesService import RegionChromosomesService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService

from utils.logger import logger

from utils.decorators import Inject

@Inject(region_results_service = RegionResultsService,
        organisms_service = OrganismsService,
        region_chromosomes_service = RegionChromosomesService,
        whole_chromosomes_service = WholeChromosomesService)
class RegionSequenceManager(SequenceManagerInterface):
    def __init__(self, sequence: Sequence = None, sequence_data: dict = None, regions: list[Sequence] = None,
                 sequence_name: str = None, organism_name: str = None, regions_number: int = 0,
                 region_results_service: RegionResultsService = None,
                 organisms_service: OrganismsService = None,
                 region_chromosomes_service: RegionChromosomesService = None,
                 whole_chromosomes_service: WholeChromosomesService = None):

        self.region_results_service = region_results_service
        self.organisms_service = organisms_service
        self.region_chromosomes_service = region_chromosomes_service
        self.whole_chromosomes_service = whole_chromosomes_service

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
                                               name=sequence.get_name(),
                                               refseq_accession_number=sequence.get_refseq_accession_number())
            self.sequence_name = self.sequence.get_name()
            self.organism_name = self.sequence.get_organism_name()

        elif sequence_data:
            self.sequence = RegionSequence(sequence_data=sequence_data, regions_number=regions_number,
                                           name=sequence_data['name'])
            self.sequence_name = sequence_data['name']
            self.organism_name = sequence_data['organism_name'] or organism_name
        elif regions:
            temp = ''.join(region.get_sequence() for region in regions)
            self.sequence = RegionSequence(sequence=temp, regions_number=regions_number, name=sequence_name)
            self.sequence_name = sequence_name
            self.organism_name = organism_name

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
        self.degree_of_multifractality = []

        # coverage
        self.cover_percentage = []
        self.cover = []

    def _attach_managers(self):
        for region in self.regions:
            self.managers.append(SequenceManager(sequence=region, sequence_name=region.get_name()))

    def generate_mfa(self):
        for manager in self.managers:
            manager.generate_mfa()
            self.mfa_results.append(manager.get_mfa_results())
            self.fq.append(manager.get_fq())

    def generate_degree_of_multifractality(self):
        for manager in self.managers:
            manager.generate_degree_of_multifractality()
            self.degree_of_multifractality.append(manager.get_degree_of_multifractality())

    def graph_cgr(self):
        for manager in self.managers:
            manager.graph_cgr(name=f"{self.organism_name}/regions")

    def graph_3d_cgr(self, grid_size=512):
        for manager in self.managers:
            manager.graph_3d_cgr(grid_size, name=f"{self.organism_name}/regions")

    def graph_linear_fit(self):
        for index, regional_fq in enumerate(self.fq):
            Graphs.graph_linear_fit(fq_values=regional_fq,
                                    epsilons=self.managers[index].get_mfa_generator().get_epsilons(),
                                    sequence_name=self.sequence.get_regions_names()[index],
                                    name=f"{self.organism_name}/regions")

    def graph_degree_of_multifractality(self, y_range=None, top_labels=True):
        # Check if the lengths of x_array and y_array match
        if len(self.sequence.get_regions_names()) != len(self.degree_of_multifractality):
            raise ValueError(
                f"Number of chromosome names {len(self.sequence.get_regions_names())} ({self.sequence.get_regions_names()}) \
        does not match the number of degrees of multifractality {len(self.degree_of_multifractality)} ({self.degree_of_multifractality}).")

        Graphs.graph_bars(x_array=self.sequence.get_regions_names(), y_array=self.degree_of_multifractality,
                          title=f"Degree of multifractality by regions for {self.sequence_name}",
                          name=f"{self.organism_name}/regions",
                          y_label="Degree of multifractality", y_range=y_range, top_labels=top_labels)

    def graph_multifractal_spectrum(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='Dq_values', x_label='q', y_label='Dq',
                          title=f'Dq vs q by regions for {self.sequence_name}', name=f"{self.organism_name}/regions",
                          labels_array=self.sequence.get_regions_names())

    def graph_correlation_exponent(self):
        Graphs.graph_many(results_array=self.mfa_results, X='q_values', Y='tau_q_values', x_label='q', y_label='t(q)',
                          title=f't(q) vs q by regions for {self.sequence_name}', name=f"{self.organism_name}/regions",
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
            self.graph_degree_of_multifractality(y_range=y_range_degrees_of_multifractality, top_labels=top_labels)

    def graph_coverage(self):
        for manager in self.managers:
            manager.set_organism_name(self.organism_name)
            manager.graph_coverage(subfolder="regions")


    def set_sequence_name(self, sequence_name):
        self.sequence_name = sequence_name

    def set_cover(self, cover):
        for index, manager in enumerate(self.managers):
            manager.set_cover(cover[index])

    def set_cover_percentage(self, cover_percentage):
        for index, manager in enumerate(self.managers):
            manager.set_cover_percentage(cover_percentage[index])

    def get_sequence_name(self):
        return self.sequence_name

    def get_mfa_results(self):
        return self.mfa_results

    def get_regions_names(self):
        return self.sequence.get_regions_names()

    def get_regions(self):
        return self.regions

    def get_degree_of_multifractality(self):
        return self.degree_of_multifractality

    def get_cover(self) -> list[list[int]]:
        return self.cover

    def get_cover_percentage(self) -> list[float]:
        return self.cover_percentage

    def get_managers(self) -> list[SequenceManager]:
        return self.managers

    def save_to_db_during_execution(self, GCF):
        """
        [(val1, val2), (val1, val2)]
        ["chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        [{"q_values", "Dq_values", "tau_q_values", "DDq"}]
        """
        #try:

        organism_id = self.organisms_service.extract_by_GCF(GCF=GCF)
        whole_chromosome_id = int(self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
            self.sequence.get_refseq_accession_number()))

        for index, result in enumerate(self.mfa_results):
            chromosome_id = self.region_chromosomes_service.insert(record=(result['sequence_name'],
                                                                           self.sequence.get_refseq_accession_number(),
                                                                           organism_id,
                                                                           self.cover_percentage[index],
                                                                           self.cover[index],
                                                                           self.regions_total,
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