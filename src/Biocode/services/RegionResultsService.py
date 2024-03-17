from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.managers.DBConnectionManager import DBConnectionManager


class RegionResultsService(AbstractService):
    def __init__(self):
        self.table_name = "chr_region_results"
        self.columns = ["regions_number", "chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        self.pk_column = "id"

    def extract_results(self, GCF):
        query = f"SELECT chr_region_results.id as results_id, o.name as organism_name, chromosomes.id as chromosome_id \
                    , chromosomes.name as sequence_name, DDq, Dq_values, tau_q_values, cover, cover_percentage FROM chr_region_results  JOIN chromosomes ON " \
                f"chr_region_results.chromosome_id = chromosomes.id JOIN organisms o on chromosomes.organism_id = o.id " \
                f"WHERE o.GCF = '{GCF}';"
        return DBConnectionManager.extract_with_custom_query(query)

