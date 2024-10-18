from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.services_context.service_decorator import Service


@Service
class RegionResultsService(AbstractService):
    def __init__(self):
        self.table_name = "chr_region_results"
        self.columns = ["region_chromosome_id", "Dq_values", "tau_q_values", "DDq"]
        self.pk_column = "id"

    def extract_results(self, GCF):
        query = f"SELECT chr_region_results.id as results_id, o.name as organism_name, rc.id as chromosome_id \
                    , rc.name as sequence_name, DDq, Dq_values, tau_q_values, cover, cover_percentage FROM chr_region_results  JOIN region_chromosomes AS rc ON " \
                f"chr_region_results.region_chromosome_id = rc.id JOIN organisms o on rc.organism_id = o.id " \
                f"WHERE o.GCF = '{GCF}';"
        return self.extract_with_custom_query(query)

    def extract_ddq_by_refseq_accession_number(self, refseq_accession_number: str):
        query = f"SELECT rc.id, crr.DDq, rc.refseq_accession_number AS region_ran, " \
                f"wc.refseq_accession_number AS whole_ran, rc.region_number, rc.regions_total, rc.SIZE " \
                f"FROM chr_region_results crr JOIN region_chromosomes rc ON crr.region_chromosome_id = rc.id " \
                f"JOIN whole_chromosomes wc ON rc.whole_chromosome_id = wc.id " \
                f"WHERE whole_ran='{refseq_accession_number}' " \
                f"ORDER BY rc.region_number;"
        return self.extract_with_custom_query(query)
