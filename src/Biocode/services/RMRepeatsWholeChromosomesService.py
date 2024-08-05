from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class RMRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "RM_repeats_whole_chromosomes"
        self.columns = [
            "repeats_id", "whole_chromosomes_id", "sw_score", "percentage_divergence", "percentage_deletions",
            "percentage_insertions", "query_begin", "query_end", "repeat_length", "query_left", "strand",
            "repeat_begin", "repeat_end", "repeat_left"]

        self.pk_column = "id"

    def extract_info_by_chromosome(self, refseq_accession_number: str):
        query = f"SELECT sw_score, percentage_divergence, percentage_deletions, percentage_insertions, " \
                f"refseq_accession_number, query_begin, query_end, repeat_length, query_left, strand, repeats.name, " \
                f"repeats.class_family, repeat_begin, repeat_end, repeat_left FROM repeats " \
                f"JOIN RM_repeats_whole_chromosomes Rrwc on repeats.id = Rrwc.repeats_id JOIN whole_chromosomes wc " \
                f"on Rrwc.whole_chromosomes_id = wc.id WHERE refseq_accession_number = '{refseq_accession_number}';"
        return DBConnectionManager.extract_with_custom_query(query)