from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class WholeResultsService(AbstractService):
    def __init__(self):
        self.table_name = "chr_whole_results"
        self.columns = ["whole_chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        self.pk_column = "id"

    def extract_results(self, GCF):
        query = f"SELECT chr_whole_results.id as results_id, o.name as organism_name, whole_chromosomes.id as whole_chromosome_id \
                        , whole_chromosomes.name as sequence_name, DDq, o.GCF AS GCF, Dq_values, tau_q_values, cover, " \
                f"cover_percentage FROM chr_whole_results  JOIN whole_chromosomes ON " \
                f"chr_whole_results.whole_chromosome_id = whole_chromosomes.id JOIN organisms o on whole_chromosomes.organism_id = o.id " \
                f"WHERE o.GCF = '{GCF}';"
        return self.extract_with_custom_query(query)
