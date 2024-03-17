from src.Biocode.managers.SequenceManagerInterface import SequenceManagerInterface
from src.Biocode.sequences.Sequence import Sequence

from src.Biocode.mfa.MFA import MFA
from src.Biocode.graphs.Graphs import Graphs


class SequenceManager(SequenceManagerInterface):
    def __init__(self, sequence: Sequence = None, sequence_data: dict = None, sequence_name: str = None,
                 organism_name: str = None):
        self.organism_name = organism_name
        if sequence:
            self.sequence = sequence
            self.sequence_name = sequence_name
            if sequence_name:
                self.sequence.set_name(sequence_name)
        elif sequence_data:
            self.sequence = Sequence(sequence_data=sequence_data)
            self.sequence_name = sequence_data['name']

        self.mfa_generator = MFA(sequence=self.sequence)
        self.mfa_results = None
        self.fq = None
        self.degree_of_multifractality = None
        self.cover = None
        self.cover_percentage = None

    def generate_mfa(self):
        self.mfa_results = self.mfa_generator.multifractal_discrimination_analysis()
        self.fq = self.mfa_generator.get_fq()

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